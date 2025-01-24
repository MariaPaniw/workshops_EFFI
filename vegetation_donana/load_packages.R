# Function to check, install, and load required packages
install_and_load_packages <- function(required_packages) {
  # Check which packages are not installed
  missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  
  # Install missing packages
  if(length(missing_packages) > 0) {
    message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages, dependencies = TRUE)
  }
  
  # Load all required packages
  for(package in required_packages) {
    library(package, character.only = TRUE)
  }
  
  message("All required packages have been installed and loaded.")
}

# Example usage:
required_packages <- required_packages <- c("tidyverse", 
                                            "signal",
                                            "zoo",
                                            "lubridate",
                                            "ggplot2",
                                            "here",
                                            "gtools",
                                            "broom",
                                            "FactoMineR",
                                            "factoextra",
                                            "ggrepel",
                                            "ggpubr")
install_and_load_packages(required_packages)