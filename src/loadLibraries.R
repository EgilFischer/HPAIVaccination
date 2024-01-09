# loading required libraries ####
packages <- c("ggplot2","dplyr","tidyverse","readr","rootSolve","readxl","distributionsrd","dbscan","ggh4x")


## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

# User-defined generic functions ####
rm.exist <- function(x){
  if(exists(deparse(substitute(x))))rm(x); 
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
