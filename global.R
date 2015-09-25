
# Guillaume Lobet - University of Liege


# ROOT-FIT
# Julkowska et al. (2014). 
# Capturing Arabidopsis root architecture dynamics with ROOT-FIT reveals 
# diversity in responses to salinity. 
# Plant Physiology. doi:10.1104/pp.114.248963



# Global libraries
  packages <- c("plyr", "ggplot2", "gridExtra", "multcomp")
  for(p in packages){
    if (!require(p,character.only = TRUE)){
      install.packages(p,dep=TRUE)
      if(!require(p,character.only = TRUE))stop("Package not found")
    }
  }
  
  fitting <<- "null"
  
  
  
  
  

  
  
  
  