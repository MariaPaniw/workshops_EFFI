
#### Script to cut large NDVI raster for entire Doñana Area to the size of the Doñana Biological Reserve
#### Author: Maria Paniw
### Date: 07 January 2025

library(sf)
library(terra)
library(raster)
library(dplyr)
library(ggplot2)
library(readr)

# Set your working directory
setwd("/Users/maria/Dropbox/collaborations/EEFI/workshop")

############### 
#1. LOAD AND CUT NDVI TO RESERVE POLYGON

ndvi.names <- list.files("NDVI", pattern="*.tif")

crdref <- "+proj=utm +zone=30 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

newcrs <- "+proj=utm +zone=29 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

for(i in 1:length(ndvi.names)){
  # Load each rater
  temp <- raster(paste0("NDVI/", ndvi.names[i]))
  
  #Cut exactly to the shape of the reserve (res)
  temp2 <- crop(temp, extent(res))
  temp3 <- mask(temp2, res)
  
  # Deine CRS
  crs(temp3) <- crdref
  
  # Project to 29 zone as this is what the other data are in
  temp_p=projectRaster(temp3,crs=newcrs)
  
  writeRaster(temp_p,paste0("NDVI_reserve/",ndvi.names[i]))
}

