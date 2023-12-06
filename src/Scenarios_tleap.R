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
source("./src/multitypetransitionsSEIR_tleap.R")
if(!dir.exists("./output")){dir.create("./output")}
# SCENARIOS #
#Layers####
#Baseline parameters for layer flock ####
param.list.baseline.layer <- list(
  scenario = "Test", #scenario
  runs =5, #number of runs
  max.time = 17*30,#length of the run
  itypes = itypes, #types####
  deltat = 0.05, #time step#
  N0 = 32000, #population size####
  initial= 10 , #initially infected
  p.protect = 1, #1 - 6/26,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  latency.period = c(.1,.1),#choose very short latency period
  k.latency =2,##k parameter Erlang distribution
  infectious.period = c(3.0,4.0),#Duration infectious period 
  #variance.infectious.period = c(3.0,4.0)^2, #Variance infectious period
  k.infectious = 20, #k parameter Erlang distribution
  transRate = matrix(c(0,0.012,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  trans.mean = 514, #only considers transistion from type 2 to 1
  trans.var = 85^2, #only considers transistion from type 2 to 1
  pdie = c(1.0,0.01),#probability of dying at end of infectious period
  mortRate = 0,#0.0005 #per capita death rate #Mortality events
  intro.time = 0 #time at which the infection is introduced
)

#Size scenarios
scenario.list.size.vaccination <- list()
for(i in c(1:3)){
  for(j in c(1:6)){
  param.list <- param.list.baseline.layer;
  param.list$runs <- 50;
  param.list$N0 <-  c(15000,32000,64000)[i]
  param.list$p.hightitre <- c(0.0, 0.5,0.6,0.7,0.8,0.9)[j]
  param.list$scenario <- paste0("layerSize",param.list$N0,"Vac",param.list$p.hightitre*100);
  scenario.list.size.vaccination[[7*(i-1)+j]]<- param.list
  }
}

for(i in  c(1:length(scenario.list.size.vaccination))){
  print(scenario.list.size.vaccination[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.size.vaccination[[i]],{
    #adjust p.protect with waning
    p.protect.adjust <- p.protect*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE);
    
     return(list(
      L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected
      I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
      R0 = c(0,0), # number of initially recovered
      S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
    )
    })
  sim.out<- sim.multitypeSEIR_tleap_distwaning(scenario.list.size.vaccination[[i]],
                                               inits.gamma,
                                               gamma.waning.distribution)
  #record output and parameters
  op <- list(out = sim.out, pars = scenario.list.size.vaccination[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y%m%d"),"output",scenario.list.size.vaccination[[i]]$scenario,".RDS"))
  print(Sys.time()-start)
  }



#Clinical protection scenarios
# Clinical protection ####
scenario.list.clinic <-list()
for(i in c(1:3)){
  param.list <- param.list.baseline.layer;
  param.list$runs <- 10;
  param.list$transRate <-  0*param.list$transRate
  param.list$p.hightitre <- c(0,0.5,0)[i]
  param.list$pdie <- c(c(0.5,0.5,0.001)[i],0.001)
  param.list$scenario <- paste0("layerClinic",param.list$pdie[1]*100 ,"Phigh",param.list$p.hightitre*100);
  scenario.list.clinic[[i]]<- param.list
}


for(i in c(1:length(scenario.list.clinic))){
  print(scenario.list.clinic[[i]]$scenario);
  print(Sys.time())
  simulate.multitypeSIR(scenario.list.clinic[[i]])
}



#Waning same strain####
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 50, 100, 200,300,400,500)  #initial high-titre = exp(-transRate*introtime)
1-pgamma(intro.time, shape = 36,rate =0.07)
#Waning scenarios
scenario.list.waning <- list()
for(i in c(1:length(intro.time))){
    t0 <-intro.time[i];
    param.list <- param.list.baseline.layer;
    param.list$transShape <- matrix(c(Inf,30,Inf,Inf), nrow = 2) #value based on Rudolf et al 2010
    param.list$transRate <- matrix(c(0,0.4,0.0,0), nrow = 2) #value based on Rudolf et al 2010
    param.list$runs <- 10;
    param.list$p.hightitre <- 1-pgamma(intro.time[i], shape = 36,rate =0.07)
    param.list$max.time <-param.list$max.time;
    param.list$intro.time <- t0;
    param.list$scenario <- paste0("layerStartTime",t0);
    scenario.list.waning[[i]]<- param.list
}



for(i in c(1:length(scenario.list.waning))){
  print(scenario.list.waning[[i]]$scenario);
  print(Sys.time())
  simulate.multitypeSIR(scenario.list.waning[[i]])
}







#Waning different strain ####
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 50, 100, 200,300,400,500)  #initial high-titre = exp(-transRate*introtime)
0.7*(1-pgamma(intro.time, shape = 4,rate =0.1/7))
#Waning scenarios
scenario.list.waning <- list()
for(i in c(1:length(intro.time))){
  t0 <-intro.time[i];
  param.list <- param.list.baseline.layer;
  param.list$transShape <- matrix(c(Inf,4,Inf,Inf), nrow = 2) #value based on Rudolf et al 2010
  param.list$transRate <- matrix(c(0,0.1/7,0.0,0), nrow = 2) #value based on Rudolf et al 2010
  param.list$runs <- 10;
  param.list$p.hightitre <- .7*(1-pgamma(intro.time[i], shape = 4,rate =0.1/7))
  param.list$max.time <-param.list$max.time;
  param.list$intro.time <- t0;
  param.list$scenario <- paste0("layerStartTime",t0,"DifStrain");
  scenario.list.waning[[i]]<- param.list
}



for(i in c(1:length(scenario.list.waning))){
  print(scenario.list.waning[[i]]$scenario);
  print(Sys.time())
  simulate.multitypeSIR(scenario.list.waning[[i]])
}








#Broilers ####
#Baseline parameters for broiler flock ####
param.list.baseline.broiler <- list(
  scenario = "baseline_Broiler", #scenario
  runs = 10, #number of runs
  max.time = 46,#length of the run
  itypes = 2, #type
  N0 = 38000, #population size - agramatie website
  initial= 10 , #initially infected - choosen value
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = 2),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976 and Gemeraard et al 2023 #Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2/20, #Variance infectious period - Hobbelen et al uses shape parameter of 20 -> variance = m^2 / shape 
  transShape = matrix(c(0,Inf,0.0,0), nrow = 2), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  transRate = matrix(c(0,0,0.0,0), nrow = 2), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.99,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005/7 #per capita death rate #Mortality events - based on performance reports and pers. com. mieke matthijs 0.5% per week +Gonzales & elbers figure 1 
)

#Size scenarios
scenario.list.size.vaccination <- list()
for(i in c(1:3)){
  for(j in c(1:6)){
    param.list <- param.list.baseline.broiler;
    param.list$runs <- 10;
    param.list$N0 <-  c(20000, 38000, 73000)[i]
    param.list$p.hightitre <- c(0.0, 0.5,0.6,0.7,0.8,0.9)[j]
    param.list$scenario <- paste0("broilerSize",param.list$N0,"Vac",param.list$p.hightitre*100);
    scenario.list.size.vaccination[[6*(i-1)+j]]<- param.list
  }
}

for(i in  c(1:length(scenario.list.size.vaccination))){
  print(scenario.list.size.vaccination[[i]]$scenario);
  print(Sys.time())
  simulate.multitypeSIR(scenario.list.size.vaccination[[i]])
}


# Clinical protection ####
scenario.list.clinic <-list()
for(i in c(1:3)){
  param.list <- param.list.baseline.broiler;
  param.list$runs <- 10;
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
intro.time <- c(0, 1, 10, 50, 100,200,300,400)  #initial high-titre = exp(-transRate*introtime)

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
