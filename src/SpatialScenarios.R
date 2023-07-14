#########################################################
#                                                        
#                 Spatial Scenarios                            
#                  Description                           
#                                                        
#                  Author:                               
#                  Contact:                              
#                  Creation date                         
#########################################################
source("./src/loadLibraries.R")
source("./src/probMajorOutbreak.R")
#load location data ####
param.list.baseline.layer <- list(
  scenario = "baseline_Layer", #scenario
  runs = 10, #number of runs
  max.time = 17*30,#length of the run
  itypes = 2, #type
  N0 = 45000, #population size - agramatie website
  initial= 10 , #initially infected - choosen value
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = 2),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976 and Gemeraard et al 2023 #Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2/20, #Variance infectious period - Hobbelen et al uses shape parameter of 20 -> variance = m^2 / shape 
  transRate = matrix(c(0,0.0,0.0,0), nrow = 2), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.99,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005/7 #per capita death rate #Mortality events - based on performance reports and pers. com. mieke matthijs 0.5% per week +Gonzales & elbers figure 1 
)

source("./src/SpatialSellkeModel_init_simuls.R")
library("moments")
spatial.input%>%reframe(.by = TYPE, 
                        meanSize = mean(SIZE),
                        medianSize = median(SIZE),
                        skewnessSize = skewness(SIZE),
                        kurtosisSize = kurtosis(SIZE),
                        perc25 = quantile(SIZE,0.25),
                        perc75 = quantile(SIZE,0.75), 
                        perc95 = quantile(SIZE,0.95), 
                        seSize = sqrt(var(SIZE)/length(SIZE)),
                        sdSize = sqrt(var(SIZE)))


jarque.test(unlist(spatial.input%>%filter(TYPE == "LAYER")%>%select(SIZE)))
jarque.test(unlist(spatial.input%>%filter(TYPE == "BROILER")%>%select(SIZE)))





