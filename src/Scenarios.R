#########################################################
#                                                        
#                  Simulate scenarios                               
#                                                        
#                  Author:Egil Fischer                              
#                  Contact:           e.a.j.fischer@uu.nl                   
#                  Creation date: 2-5-2023                         
#########################################################

#load libraries
source("./src/loadLibraries.R") 
source("./src/postprocesSimulations.R")
source("./src/multitypetransitionsSIRsellke_simfunc.R")
#define number of types
itypes = 2;

# SCENARIOS #

#Layers####

#Baseline parameters for layer flock ####
param.list.baseline.layer <- list(
  scenario = "baseline_Layer", #scenario
  runs = 10, #number of runs
  max.time = 17*30,#length of the run
  itypes = itypes, #type
  N0 = 45000, #population size - agramatie website
  initial= 10 , #initially infected - choosen value
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976 and Gemeraard et al 2023 #Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2/20, #Variance infectious period - Hobbelen et al uses shape parameter of 20 -> variance = m^2 / shape 
  transRate = matrix(c(0,0.0,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.99,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005/7 #per capita death rate #Mortality events - based on performance reports and pers. com. mieke matthijs 0.5% per week +Gonzales & elbers figure 1 
)


#Clinical protection scenarios
# Clinical protection ####
scenario.list.clinic <-list()
for(i in c(1:3)){
  param.list <- param.list.baseline.layer;
  param.list$runs <- 10;
  param.list$transRate <-  0*param.list$transRate
  param.list$p.hightitre <- c(0,0.5,0)[i]
  param.list$pdie <- c(c(0.5,0.5,0.005)[i],0.005)
  param.list$scenario <- paste0("layerClinic",param.list$pdie[1]*100 ,"Phigh",param.list$p.hightitre*100);
  scenario.list.clinic[[i]]<- param.list
}

simulate.multitypeSIR(scenario.list.clinic[[3]])

for(i in c(1:length(scenario.list.clinic))){
  print(scenario.list.clinic[[i]]$scenario);
  print(Sys.time())
  simulate.multitypeSIR(scenario.list.clinic[[i]])
}



#Create other scenarios
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 1, 10, 100)  #initial high-titre = exp(-transRate*introtime)
waning <- c(0,0.012,0.038)

#create scenarios
scenario.list.clinic <-list()
for(i in c(1:3)){
  param.list <- param.list.baseline.layer;
  param.list$runs <- 10;
  param.list$transRate <-  0*param.list$transRate
  param.list$p.hightitre <- c(0,0.5,0)[i]
  param.list$pdie <- c(c(0.5,0.5,0.005)[i],0.005)
  param.list$scenario <- paste0("layerClinic",param.list$pdie[1]*100 ,"Phigh",param.list$p.hightitre*100);
  scenario.list.clinic[[i]]<- param.list
}



#Broilers ####
#Baseline parameters for broiler flock ####
param.list.baseline.broiler <- list(
  scenario = "baseline_Layer", #scenario
  runs = 10, #number of runs
  max.time = 46,#length of the run
  itypes = itypes, #type
  N0 = 75000, #population size - agramatie website
  initial= 10 , #initially infected - choosen value
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976 and Gemeraard et al 2023 #Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2/20, #Variance infectious period - Hobbelen et al uses shape parameter of 20 -> variance = m^2 / shape 
  transRate = matrix(c(0,0.0,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.99,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005/7 #per capita death rate #Mortality events - based on performance reports and pers. com. mieke matthijs 0.5% per week +Gonzales & elbers figure 1 
)


# Clinical protection ####
scenario.list.clinic <-list()
for(i in c(1:3)){
  param.list <- param.list.baseline.broiler;
  param.list$runs <- 10;
  param.list$transRate <-  0*param.list$transRate
  param.list$p.hightitre <- c(0,0.5,0)[i]
  param.list$pdie <- c(c(0.5,0.5,0.005)[i],0.005)
  param.list$scenario <- paste0("broilerClinic",param.list$pdie[1]*100 ,"Phigh",param.list$p.hightitre*100);
  scenario.list.clinic[[i]]<- param.list
}

for(i in c(1:length(scenario.list.clinic))){
  print(scenario.list.clinic[[i]]$scenario);
  print(Sys.time())
  simulate.multitypeSIR(scenario.list.clinic[[i]])
}


#Create other scenarios####
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 1, 10, 30)  #initial high-titre = exp(-transRate*introtime)
waning <- c(0,0.012,0.038)

#create scenarios
scenario.list <- list()
for(i in c(1:4)){
  for(j in c(1:3)){
    param.list <- param.list.baseline.broiler;
    param.list$max.time  <- param.list.baseline.broiler$max.time - intro.time[i];
    param.list$transRate[2,1]  <- waning[j];
    param.list$p.hightitre <- exp(-waning[j]*intro.time[i]);
    param.list$scenario <- paste0("broilerwaning",waning[j],"introtime",intro.time[i]);
    scenario.list[[(i-1)*length(waning) + j]]<- param.list
  }
  
}




#do simulations#
output.baseline.broiler<- simulate.multitypeSIR(param.list.baseline.broiler)
for(i in c(1:length(scenario.list))){
  print(scenario.list[[i]]$scenario);
  print(Sys.time())
  simulate.multitypeSIR(scenario.list[[i]])
}





#Broilers ####
#Baseline parameters for broiler flock ####
param.list.baseline.broiler <- list(
  scenario = "baseline_Layer", #scenario
  runs = 10, #number of runs
  max.time = 46,#length of the run
  itypes = itypes, #type
  N0 = 75000, #population size - agramatie website
  initial= 10 , #initially infected - choosen value
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976 and Gemeraard et al 2023 #Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2/20, #Variance infectious period - Hobbelen et al uses shape parameter of 20 -> variance = m^2 / shape 
  transRate = matrix(c(0,0.0,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.99,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005/7 #per capita death rate #Mortality events - based on performance reports and pers. com. mieke matthijs 0.5% per week +Gonzales & elbers figure 1 
)


#Create other scenarios
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 1, 10, 30)  #initial high-titre = exp(-transRate*introtime)
waning <- c(0,0.012,0.038)

#create scenarios
scenario.list <- list()
for(i in c(1:4)){
  for(j in c(1:3)){
    param.list <- param.list.baseline.broiler;
    param.list$max.time  <- param.list.baseline.broiler$max.time - intro.time[i];
    param.list$transRate[2,1]  <- waning[j];
    param.list$p.hightitre <- exp(-waning[j]*intro.time[i]);
    param.list$scenario <- paste0("broilerwaning",waning[j],"introtime",intro.time[i]);
    scenario.list[[(i-1)*length(waning) + j]]<- param.list
  }
  
}




#do simulations#
output.baseline.broiler<- simulate.multitypeSIR(param.list.baseline.broiler)
for(i in c(1:length(scenario.list))){
  print(scenario.list[[i]]$scenario);
  print(Sys.time())
  simulate.multitypeSIR(scenario.list[[i]])
}
