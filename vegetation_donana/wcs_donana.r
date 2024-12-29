
library(ows4R)
library(sf)
library(terra)
library(dplyr)
library(ggplot2)
library(readr)

# Set your working directory
setwd("/Users/maria/Dropbox/collaborations/EEFI/workshop")

## The script has the following sections:

# 1. Load the spatial data. I load the coordinates of the monitoring plots (and convert them into spatial points) and a series placeholder flooding satellite-based rasters using a WCS client
# 2. Process the spatial data: I get the value of the raster at a given point to be used as a covariate to predict abundacne change 
    # NOTE: Because flooding values are 0, I created a placeholder covariate: density of interspecific neighbors in a 5x5m plot
# 3. Model abundance change. Here, I just use a simple autoregerssive model done in JAGS
# 4. Do the forecast and quantify forecast skill
############### 
#1. LOAD SPATIAL DATA

### Load location of the study plots with abundace monitoring since 2007:
### Reference system is WGS84
coords_long=read_csv(url("https://raw.githubusercontent.com/MariaPaniw/workshops_EFFI/refs/heads/main/vegetation_donana/coords_plot_since2007.csv?token=GHSAT0AAAAAAC2TAOO4TOA3HUOY7SYFFRE4Z3RUREQ"))

crdref <- "+proj=longlat +datum=WGS84"
pts_long <- vect(cbind(coords_long$Long,coords_long$Lat), atts=coords_long[,c("ID","Elevation")],crs=crdref)
pts_long

newcrs <- "+proj=utm +zone=29 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

rob <- project(pts_long, newcrs)

### Load the hydro period maps
WCS <- WCSClient$new("https://icts-donana.es/geoserver/hidroperiodov1/wms", "2.0.1")
caps <- WCS$getCapabilities()


############### 
#2. GET SPATIAL COVARIATE DATA AT POINT COORDINATES

### go through maps (37 years): 
df.hydro=NULL 

for(i in 1:37){
 
  print(paste("Running year ",i))
  name=caps$getCoverageSummaries()[[i]]$CoverageId # ID of a map
  
  hydro_map <- caps$findCoverageSummaryById(name, exact = T)
  hydro_map_dims <- hydro_map$getDimensions()
  hydro_map_des <- hydro_map$getDescription()
  
  
  cov_data <- hydro_map$getCoverage(
    # bbox = OWSUtils$toBBOX(700000,408000,740000, 4130000),
    # time = hydro_map_dims[[1]]$uom[[1]],  
    # filename = paste("/Users/maria/Dropbox/collaborations/EEFI/workshop/maps/hydro/hydro_",i,".tiff",sep="")
    
  )
  
  #Get hydrology value for point coordinates
  df.hydro=rbind(df.hydro,data.frame(plot=rob$ID,hydroperiod=extract(cov_data, rob, sp = T)[,2],year=i))
}

unique(df.hydro$hydroperiod)
unique(df.hydro$year)

df.hydro[is.na(df.hydro$hydroperiod),]

df.hydro$year=factor(df.hydro$year)
levels(df.hydro$year)=1985:2022

plot(cov_data,xlim = c(719000,722100),
     ylim = c(4097200, 4099400))

points(rob)

plot(pts_long)

### Add shrub abundance data: 
num=read_csv(url("https://raw.githubusercontent.com/MariaPaniw/workshops_EFFI/refs/heads/main/vegetation_donana/shrub_number.csv?token=GHSAT0AAAAAAC2TAOO5VDUV3XKSXUPCSRRUZ3RURUA"))

num$year=factor(num$year)

# Add satellite derived measure to abundance data 

num$cov=left_join(num,df.hydro,by=c("plot","year"))$hydroperiod

#### In the absence of the correct satellite data, we just fill the covariate data with interspecific abundances
num$cov=NA

years=unique(num$year)
plots=unique(num$plot)
species=unique(num$species)

colnames(num)[5:7]=c("adults","saplings","seedlings")

for(a in years){
  
  for(b in plots){
    
    for(c in species){
      
      if(length(num$cov[num$species%in%c&num$plot%in%b&num$year%in%a])>0){
        
        num$cov[num$species%in%c&num$plot%in%b&num$year%in%a]=sum(num$adults[!num$species%in%c&num$plot%in%b&num$year%in%a],na.rm=T)
        
      }
    }
  }
}

# we focus on two very abundant shrubs with lots of data

sub=num[num$species%in%c("Halimium halimifolium","Lavandula stoechas"),]


############### 
#3. STATISTICAL MODEL

# We could use the new SpAbudance package https://doi.org/10.1111/2041-210X.14332
# To do: make it a proper spatial model 
# Here I just use JAGS

#Create C matrix 

n.plots=18
n.years=length(2007:2022)

### Halimium
C.hal=array(NA,c(n.plots,n.years))
cov.hal=array(NA,c(n.plots,n.years))

hal=sub[sub$species%in%"Halimium halimifolium",]

hal$saplings[is.na(hal$saplings)]=0
hal$adults[is.na(hal$adults)]=0

for(i in 1:n.plots){ # site loop
  for(k in 1:n.years){
    
    year=as.character(2007:2022)[k]
    plot=unique(sub$plot)[i]
  
    sum=as.numeric(hal[hal$plot%in%plot&hal$year%in%year,c("adults")])
    cov= hal$cov[hal$plot%in%plot&hal$year%in%year]
    if(length(sum)>0) C.hal[i,k]=sum
    
    if(length(cov)>0)  cov.hal[i,k]=cov
    
  }
  
}

C.hal[5,16]=0

# We have to interpolate the neighbors

cov.hal[,3]=round(rowMeans(cov.hal[,c(2,4)]))
cov.hal[,5:6]=round(rowMeans(cov.hal[,c(4,7)]))
cov.hal[,8:10]=round(rowMeans(cov.hal[,c(7,11)]))
cov.hal[,12]=round(rowMeans(cov.hal[,c(11,13)]))

cov.hal[5,16]=10

## Lavandula
C.lav=array(NA,c(n.plots,n.years))
cov.lav=array(NA,c(n.plots,n.years))
lav=sub[sub$species%in%"Lavandula stoechas",]

lav$saplings[is.na(lav$saplings)]=0
lav$adults[is.na(lav$adults)]=0

for(i in 1:n.plots){ # site loop
  for(k in 1:n.years){
    
    year=as.character(2007:2022)[k]
    plot=unique(sub$plot)[i]
    
    sum=as.numeric(lav[lav$plot%in%plot&lav$year%in%year,c("adults")])
    cov= lav$cov[lav$plot%in%plot&lav$year%in%year]
    
    if(length(sum)>0) C.lav[i,k]=sum
    
    if(length(cov)>0)  cov.lav[i,k]=cov
    
  }
  
}

C.lav[5,16]=0

# We have to interpolate the neighbors

cov.lav[,3]=round(rowMeans(cov.lav[,c(2,4)]))
cov.lav[,5:6]=round(rowMeans(cov.lav[,c(4,7)]))
cov.lav[,8:10]=round(rowMeans(cov.lav[,c(7,11)]))
cov.lav[,12]=round(rowMeans(cov.lav[,c(11,13)]))

cov.lav[5,16]=10

#### Bundle data for JAGS model

bdata <- list(n.hal = C.hal, C.hal=C.hal,
              n.lav = C.lav, C.lav=C.lav,
              nsites = dim(C.hal)[1],
              cov.hal=cov.hal,
              cov.lav=cov.lav,
              nyears = dim(C.hal)[2]
)

# Specify model in BUGS language
cat(file = "abundance change shrubs.txt","
    model {
    # Priors

    a0.h ~ dnorm( 0 , 1.0E-05 )
    a1.h ~ dnorm( 0 , 1.0E-05 ) 
    a2.h ~ dnorm( 0 , 1.0E-05 )
  
    a0.l ~ dnorm( 0 , 1.0E-05 )
    a1.l ~ dnorm( 0 , 1.0E-05 ) 
    a2.l ~ dnorm( 0 , 1.0E-05 )

    for(i in 1:nsites) { # Loop over sites
    
    #Initial abundance
    
    N.h[i,1] <- C.hal[i,1]
    
    N.l[i,1] <- C.lav[i,1]
  

    #Specify the model for years 2 through nYears
    
    for(t in 1:(nyears-1)) {
    
    
# Halimium halimifolium

    n.hal[i,t+1] ~ dpois(N.h[i,t+1])


    log(N.h[i,t+1]) <- a0.h + a1.h *cov.hal[i,t]  + a2.h * log(N.h[i,t]+0.001) 
 
# Lavandula stoechas 

    n.lav[i,t+1] ~ dpois(N.l[i,t+1])
 
     log(N.l[i,t+1]) <- a0.l + a1.l *cov.lav[i,t]  + a2.l * log(N.l[i,t]+0.001) 
  

    # # Goodness of fit

    n.pred[i,t+1] ~ dpois(N.h[i,t+1]) # Part 2 of HM
    e1[i,t+1] <- N.h[i,t+1]
    resid1[i,t+1] <- pow(pow(n.hal[i,t+1], 0.5) - pow(e1[i,t+1], 0.5), 2)
    resid1.pred[i,t+1] <- pow(pow(n.pred[i,t+1], 0.5) - pow(e1[i,t+1], 0.5), 2)

    }

    }
     
    fit1.data <- sum(resid1[,2:nyears]) # Fit statistic 
    fit1.pred <- sum(resid1.pred[,2:nyears])
    
    for( t in 2:nyears){
    Ntot.h[t] <- sum(N.h[,t])
   
    
    Ntot.l[t] <- sum(N.l[,t])
    
    
    }
    
}
    ")

inits <- function(){list(a0.h=rnorm(1,0,0.01),
                         a1.h=rnorm(1,0,0.01),
                         a2.h=rnorm(1,0,0.01),
                         
                         a0.l=rnorm(1,0,0.01),
                         a1.l=rnorm(1,0,0.01),
                         a2.l=rnorm(1,0,0.01))}

# What to monitor 

params <- c( "a0.h",
             "a1.h",
             "a2.h",
             "Ntot.h",
             "N.h",
             "a0.l",
             "a1.l",
             "a2.l",
             "Ntot.l",
             "N.l",
             "fit1.data",
             "fit1.pred")

library(jagsUI)

na <- 10000 ; ni <- 450000 ; nt <- 500 ; nb <- 100000 ; nc <- 3


out1 <- jags(bdata, inits, params, "abundance change shrubs.txt", n.adapt = na, n.chains = nc,
             n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T)

traceplot(out1) # the output is not too bad; chains look convereged

print(out1, 3) # Looks good

# Bayesian p-value

mean(out1$sims.list$fit1.pred > out1$sims.list$fit1.data)

# Predicted values underestimate observed (for Halimium)
plot(out1$sims.list$fit1.pred, out1$sims.list$fit1.data)
abline(1,1)

out1$mean

## Plot observed vs Predicted

# Total abundance: 

# HALIMIUM
bay.out=out1$sims.list$Ntot.h[sample(c(1:dim(out1$sims.list$Ntot.h)[1]),100),] # take a random sample of 100 posterior values

df.pred=data.frame(N=as.numeric(bay.out[,-1]),
                   year=rep(c(2008:2022),each=dim(bay.out)[1]),
                   sim=1:dim(bay.out)[1])

hal.tot=aggregate(cbind(adults,saplings)~year,data=hal,function(x) sum(x,na.rm=T))
hal.tot$sim=NA
hal.tot$tot=hal.tot$adults
hal.tot$year=as.numeric(as.character(hal.tot$year))

a=ggplot(df.pred,aes(year,N,group=sim))+
  geom_line(alpha=0.1)+
  geom_point(data=hal.tot, aes(year,tot),size=3,col="blue")+
  # scale_color_viridis(discrete = T,option="A")+
  xlab("Year")+ylab("Abundance")+theme_bw(base_size=20)+
  theme(panel.grid = element_blank())+
  ggtitle("Halimium halimifolium")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none",
        legend.title = element_blank())

a
# LAVANDULA
bay.out=out1$sims.list$Ntot.l[sample(c(1:dim(out1$sims.list$Ntot.l)[1]),100),] # take a random sample of 100 posterior values

df.pred=data.frame(N=as.numeric(bay.out[,-1]),
                   year=rep(c(2008:2022),each=dim(bay.out)[1]),
                   sim=1:dim(bay.out)[1])

lav.tot=aggregate(cbind(adults,saplings)~year,data=lav,function(x) sum(x,na.rm=T))
lav.tot$sim=NA
lav.tot$tot=lav.tot$adults
lav.tot$year=as.numeric(as.character(lav.tot$year))

b=ggplot(df.pred,aes(year,N,group=sim))+
  geom_line(alpha=0.1)+
  geom_point(data=lav.tot, aes(year,tot),size=3,col="blue")+
  # scale_color_viridis(discrete = T,option="A")+
  xlab("Year")+ylab("Abundance")+theme_bw(base_size=20)+
  theme(panel.grid = element_blank())+
  ggtitle("Lavandula stoechas")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        legend.position = "none",
        legend.title = element_blank())

b
############### 
#4. FORECAST NEXT YEAR & NEW SPACE

