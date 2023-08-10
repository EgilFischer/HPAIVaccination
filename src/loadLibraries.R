# loading required libraries ####
packages <- c("ggplot2","dplyr","tidyverse","readr","rootSolve","readxl")


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