# Functions ----
# Function to preprocess NDVI data with BISE correction and interpolation
process_ndvi_data <- function(raw_ndvi_df) {
  # Ensure data is sorted by date
  raw_ndvi_df <- raw_ndvi_df %>%
    arrange(date)
  
  # Initialize vectors for processed data
  ndvi_raw <- raw_ndvi_df$ndvi
  ndvi_processed <- ndvi_raw
  
  # 1. BISE correction
  diff_vals <- diff(ndvi_raw)
  decrease_idx <- which(diff_vals < 0) + 1  # points showing decrease
  
  # Using a slope threshold of 0.2 (20%)
  data_threshold <- c(NA, ndvi_raw[-1] - 0.2 * diff_vals)
  
  # Forward sliding period (3 values window)
  val_sliding_period <- running(ndvi_raw, fun=max, width=3)
  val_sliding_period <- c(val_sliding_period, NA, NA)
  
  # Identify points to reject
  law_check <- val_sliding_period[decrease_idx] - data_threshold[decrease_idx]
  reject_decrease <- decrease_idx[which(law_check > 0)]
  
  # Check for sudden increases
  increase_idx <- which(diff_vals > 0.2)
  reject_increase <- increase_idx[!increase_idx %in% decrease_idx]
  reject_points <- c(reject_increase, reject_decrease)
  
  # Apply corrections
  ndvi_processed[reject_points] <- NA
  
  # 2. Fill start and end gaps
  start_window <- ndvi_processed[7:11]
  ndvi_processed[1:6] <- start_window[!is.na(start_window)][1]
  
  end_window <- na.omit(ndvi_processed[(length(ndvi_processed)-4):(length(ndvi_processed)-1)])
  if(length(end_window) > 0) {
    ndvi_processed[length(ndvi_processed)] <- min(end_window)
  } else {
    ndvi_processed[length(ndvi_processed)] <- 0
  }
  
  # 3. Interpolate missing values and apply Savitzky-Golay filter
  ndvi_interpolated <- na.spline(ndvi_processed)
  ndvi_smoothed <- sgolayfilt(ndvi_interpolated, p=3, n=7, m=0)
  
  # Create daily interpolation
  dates_seq <- seq(min(raw_ndvi_df$date), max(raw_ndvi_df$date), by="days")
  daily_values <- spline(x=as.numeric(raw_ndvi_df$date), 
                         y=ndvi_smoothed, 
                         xout=as.numeric(dates_seq))$y
  
  # Return results
  result_df <- data.frame(
    date = dates_seq,
    int.NDVI = daily_values
  )
  
  return(result_df)
}

get_slope_significance <- function(data) {
  # Convert date to numeric for linear model
  data$time <- as.numeric(data$date)
  
  # Fit linear model
  model <- lm(int.NDVI ~ time, data = data)
  
  # Get model statistics
  stats <- tidy(model)
  
  # Return slope and p-value
  return(data.frame(
    slope = stats$estimate[2],
    p.value = stats$p.value[2]
  ))
}

# Function to calculate NDVI metrics
calculate_ndvi_metrics <- function(data) {
  
  # Add month column if not present
  if(!"month" %in% colnames(data)) {
    data$month <- month(data$date)
  }
  
  # Calculate metrics for each plot and year
  ndvi_metrics <- data %>%
    group_by(plot, year) %>%
    summarise(
      # Annual metrics
      integrated_ndvi = sum(int.NDVI, na.rm = TRUE),
      annual_mean = mean(int.NDVI, na.rm = TRUE),
      annual_median = median(int.NDVI, na.rm = TRUE),
      
      # Winter-Spring metrics (January-May)
      winter_spring_mean = mean(int.NDVI[month %in% 1:5], na.rm = TRUE),
      winter_spring_median = median(int.NDVI[month %in% 1:5], na.rm = TRUE),
      winter_spring_max = max(int.NDVI[month %in% 1:5], na.rm = TRUE),
      winter_spring_integrated = sum(int.NDVI[month %in% 1:5], na.rm = TRUE),
      
      # Summer metrics (July-September)
      summer_mean = mean(int.NDVI[month %in% 7:9], na.rm = TRUE),
      summer_median = median(int.NDVI[month %in% 7:9], na.rm = TRUE),
      summer_max = max(int.NDVI[month %in% 7:9], na.rm = TRUE),
      summer_integrated = sum(int.NDVI[month %in% 7:9], na.rm = TRUE),
      
      # Count number of observations in each period
      n_winter_spring = sum(month %in% 1:5),
      n_summer = sum(month %in% 7:9),
      n_total = n()
    ) %>%
    ungroup()
  
  return(ndvi_metrics)
}

# Function to calculate slope and significance
# Function to calculate slope and significance with proper percent change
calculate_trends <- function(data) {
  # Convert the data to long format for the three variables
  long_data <- data %>%
    select(plot, year, 
           integrated_ndvi, 
           summer_integrated, 
           winter_spring_integrated) %>%
    pivot_longer(cols = c(integrated_ndvi, summer_integrated, winter_spring_integrated),
                 names_to = "variable",
                 values_to = "value")
  
  # Calculate trends for each plot and variable
  trends <- long_data %>%
    group_by(plot, variable) %>%
    summarise(
      # Get first and last year values from the fitted model
      model = list(lm(value ~ year)),
      first_year = min(year),
      last_year = max(year),
      n_years = last_year - first_year,
      # Get fitted values for first and last year
      first_fitted = predict(first(model), newdata = data.frame(year = min(year))),
      last_fitted = predict(first(model), newdata = data.frame(year = max(year))),
      # Calculate total and annual percent change
      total_percent_change = ((last_fitted - first_fitted) / first_fitted) * 100,
      percent_change_per_year = total_percent_change / n_years,
      # Extract model statistics using broom
      stats = list(tidy(lm(value ~ year))),
      # Get R-squared
      rsq = summary(first(model))$r.squared,
      .groups = "drop"
    ) %>%
    # Extract slope and p-value
    mutate(
      slope = sapply(stats, function(x) x$estimate[2]),
      p_value = sapply(stats, function(x) x$p.value[2]),
      significant = p_value < 0.05
    ) %>%
    # Select final columns
    select(plot, variable, slope, p_value, significant, 
           rsq, total_percent_change, percent_change_per_year,
           first_year, last_year, first_fitted, last_fitted)
  
  return(trends)
}

