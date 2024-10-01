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
source("./src/WithinFarm/postprocesSimulations.R")
source("./src/WithinFarm/multitypetransitionsSEIR_tleap.R")
if(!dir.exists("./output")){dir.create("./output")}
# SCENARIOS #
#Layers####
#Baseline parameters for layer flock ####
param.list.baseline.layer <- list(
  scenario = "Layer", #scenario
  runs =100, #number of runs
  max.time = 17*30,#length of the run
  itypes = itypes, #types#
  deltat = 0.05, #time step#
  N0 = 32000, #population size
  initial= 10 , #initially infected
  p.hightitre = 1, #1 - 6/26,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  latency.period = c(.1,.1),#choose very short latency period
  k.latency =2,##k parameter Erlang distribution
  infectious.period = c(3.0,4.0),#Duration infectious period 
  k.infectious = 20, #k parameter Erlang distribution
#  transRate = matrix(c(0,0.012,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  trans.mean = 514, #only considers transistion from type 2 to 1
  trans.var = 85^2, #only considers transistion from type 2 to 1
  pdie = c(1.0,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005, #per capita death rate 
  intro.time = 0 #time at which the infection is introduced
)

#Size scenarios
scenario.list.size.vaccination <- list()
for(i in c(1:3)){
  for(j in c(1:6)){
  param.list <- param.list.baseline.layer;
  param.list$runs <- 100;
  param.list$N0 <-  c(15000,32000,64000)[i]
  param.list$p.hightitre <- c(0.0, 0.5,0.6,0.7,0.8,0.9)[j]
  param.list$scenario <- paste0("layerSize",param.list$N0,"Vac",param.list$p.hightitre*100);
  scenario.list.size.vaccination[[6*(i-1)+j]]<- param.list
  }
}

for(i in  c(1:length(scenario.list.size.vaccination))){
  print(scenario.list.size.vaccination[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.size.vaccination[[i]],{
    #adjust p.protect with waning
    p.protect.adjust <- p.hightitre*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE);
    
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
  param.list$runs <- 100;
  param.list$trans.mean <-  1000 #no waning
  param.list$p.hightitre <- c(0,0.5,0)[i]
  param.list$pdie <- c(c(0.5,0.5,0.001)[i],0.001)
  param.list$scenario <- paste0("layerClinic",param.list$pdie[1]*100 ,"Phigh",param.list$p.hightitre*100);
  scenario.list.clinic[[i]]<- param.list
}




for(i in  c(1:length(scenario.list.clinic))){
  print(scenario.list.clinic[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.clinic[[i]],{
    #adjust p.protect with waning
    if(is.infinite(trans.mean)){p.protect.adjust = 0}else{
    p.protect.adjust <- p.hightitre*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE)};
    
    return(list(
      L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected
      I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
      R0 = c(0,0), # number of initially recovered
      S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
    )
  })
  sim.out<- sim.multitypeSEIR_tleap_distwaning(scenario.list.clinic[[i]],
                                               inits.gamma,
                                               gamma.waning.distribution)
  #record output and parameters
  op <- list(out = sim.out, pars = scenario.list.size.vaccination[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y%m%d"),"output",scenario.list.clinic[[i]]$scenario,".RDS"))
  print(Sys.time()-start)
}


#Waning same strain####
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 50, 100, 200,300,400,500)  #

#Waning scenarios
scenario.list.waning <- list()
for(i in c(1:length(intro.time))){
    t0 <-intro.time[i];
    param.list <- param.list.baseline.layer;
    param.list$runs <- 100;
    param.list$p.hightitre <- 1
    param.list$max.time <-param.list$max.time;
    param.list$intro.time <- t0;
    param.list$scenario <- paste0("layerStartTime",t0);
    scenario.list.waning[[i]]<- param.list
}



for(i in  c(1:length(scenario.list.waning))){
  print(scenario.list.waning[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.waning[[i]],{
    #adjust p.protect with waning
    if(is.infinite(trans.mean)){p.protect.adjust = 0}else{
      p.protect.adjust <- p.hightitre*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE)};
    
    return(list(
      L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected
      I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
      R0 = c(0,0), # number of initially recovered
      S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
    )
  })
  sim.out<- sim.multitypeSEIR_tleap_distwaning(scenario.list.waning[[i]],
                                               inits.gamma,
                                               gamma.waning.distribution)
  #record output and parameters
  op <- list(out = sim.out, pars = scenario.list.waning[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y%m%d"),"output",scenario.list.waning[[i]]$scenario,".RDS"))
  print(Sys.time()-start)
}







#Waning different strain ####
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 50, 100, 200,300,400,500)  #initial high-titre = exp(-transRate*introtime)
#Waning scenarios
scenario.list.waning.hetero <- list()
for(i in c(1:length(intro.time))){
  param.list <- param.list.baseline.layer;
  param.list$trans.mean = 280; #only considers transistion from type 2 to 1
  param.list$trans.var = 140^2; #only considers transistion from type 2 to 1
  param.list$max.time <-param.list$max.time;
  param.list$intro.time <- t0;
  param.list$scenario <- paste0("layerStartTime",t0,"DifStrain");
  scenario.list.waning.hetero[[i]]<- param.list
}




for(i in  c(1:length(scenario.list.waning.hetero))){
  print(scenario.list.waning.hetero[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.waning.hetero[[i]],{
    #adjust p.protect with waning
    if(is.infinite(trans.mean)){p.protect.adjust = 0}else{
      p.protect.adjust <- p.hightitre*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE)};
    
    return(list(
      L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected
      I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
      R0 = c(0,0), # number of initially recovered
      S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
    )
  })
  sim.out<- sim.multitypeSEIR_tleap_distwaning(scenario.list.waning.hetero[[i]],
                                               inits.gamma,
                                               gamma.waning.distribution)
  #record output and parameters
  op <- list(out = sim.out, pars = scenario.list.waning.hetero[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y%m%d"),"output",scenario.list.waning.hetero[[i]]$scenario,".RDS"))
  print(Sys.time()-start)
}










#Broilers ####
#Baseline parameters for broiler flock ####
param.list.baseline.broiler <- param.list.baseline.layer;
param.list.baseline.broiler$scenario <- "Broiler";
param.list.baseline.broiler$max.time<- 46;
param.list.baseline.broiler$N0<- 38000;




#Size scenarios
scenario.list.size.vaccination <- list()
for(i in c(1:3)){
  for(j in c(1:6)){
    param.list <- param.list.baseline.broiler;
    param.list$runs <- 100;
    param.list$N0 <-  c(20000,38000,73000)[i]
    param.list$p.hightitre <- c(0.0, 0.5,0.6,0.7,0.8,0.9)[j]
    param.list$scenario <- paste0("broilerSize",param.list$N0,"Vac",param.list$p.hightitre*100);
    scenario.list.size.vaccination[[6*(i-1)+j]]<- param.list
  }
}

for(i in  c(1:length(scenario.list.size.vaccination))){
  print(scenario.list.size.vaccination[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.size.vaccination[[i]],{
    #adjust p.protect with waning
    p.protect.adjust <- p.hightitre*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE);
    
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
  param.list$runs <- 100;
  param.list$trans.mean <-  1000 #no waning
  param.list$p.hightitre <- c(0,0.5,0)[i]
  param.list$pdie <- c(c(0.5,0.5,0.001)[i],0.001)
  param.list$scenario <- paste0("layerClinic",param.list$pdie[1]*100 ,"Phigh",param.list$p.hightitre*100);
  scenario.list.clinic[[i]]<- param.list
}




for(i in  c(1:length(scenario.list.clinic))){
  print(scenario.list.clinic[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.clinic[[i]],{
    #adjust p.protect with waning
    if(is.infinite(trans.mean)){p.protect.adjust = 0}else{
      p.protect.adjust <- p.hightitre*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE)};
    
    return(list(
      L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected
      I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
      R0 = c(0,0), # number of initially recovered
      S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
    )
  })
  sim.out<- sim.multitypeSEIR_tleap_distwaning(scenario.list.clinic[[i]],
                                               inits.gamma,
                                               gamma.waning.distribution)
  #record output and parameters
  op <- list(out = sim.out, pars = scenario.list.size.vaccination[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y%m%d"),"output",scenario.list.clinic[[i]]$scenario,".RDS"))
  print(Sys.time()-start)
}


#Waning same strain####
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 50, 100, 200,300,400,500)  #

#Waning scenarios
scenario.list.waning <- list()
for(i in c(1:length(intro.time))){
  t0 <-intro.time[i];
  param.list <- param.list.baseline.layer;
  param.list$runs <- 100;
  param.list$p.hightitre <- 1
  param.list$max.time <-param.list$max.time;
  param.list$intro.time <- t0;
  param.list$scenario <- paste0("layerStartTime",t0);
  scenario.list.waning[[i]]<- param.list
}



for(i in  c(1:length(scenario.list.waning))){
  print(scenario.list.waning[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.waning[[i]],{
    #adjust p.protect with waning
    if(is.infinite(trans.mean)){p.protect.adjust = 0}else{
      p.protect.adjust <- p.hightitre*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE)};
    
    return(list(
      L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected
      I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
      R0 = c(0,0), # number of initially recovered
      S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
    )
  })
  sim.out<- sim.multitypeSEIR_tleap_distwaning(scenario.list.waning[[i]],
                                               inits.gamma,
                                               gamma.waning.distribution)
  #record output and parameters
  op <- list(out = sim.out, pars = scenario.list.waning[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y%m%d"),"output",scenario.list.waning[[i]]$scenario,".RDS"))
  print(Sys.time()-start)
}







#Waning different strain ####
#introduction at different times since peak protection at t = 0
intro.time <- c(0, 50, 100, 200,300,400,500)  #initial high-titre = exp(-transRate*introtime)
#Waning scenarios
scenario.list.waning.hetero <- list()
for(i in c(1:length(intro.time))){
  t0 <-intro.time[i];
  param.list <- param.list.baseline.layer;
  param.list$trans.mean = 280; #only considers transistion from type 2 to 1
  param.list$trans.var = 140^2; #only considers transistion from type 2 to 1
  param.list$max.time <-param.list$max.time;
  param.list$intro.time <- t0;
  param.list$scenario <- paste0("layerStartTime",t0,"DifStrain");
  scenario.list.waning.hetero[[i]]<- param.list
}




for(i in  c(1:length(scenario.list.waning.hetero))){
  print(scenario.list.waning.hetero[[i]]$scenario);
  print(start<- Sys.time()); 
  inits.gamma <- with(scenario.list.waning.hetero[[i]],{
    #adjust p.protect with waning
    if(is.infinite(trans.mean)){p.protect.adjust = 0}else{
      p.protect.adjust <- p.hightitre*pgamma(intro.time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean, lower.tail = FALSE)};
    
    return(list(
      L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected
      I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
      R0 = c(0,0), # number of initially recovered
      S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
    )
  })
  sim.out<- sim.multitypeSEIR_tleap_distwaning(scenario.list.waning.hetero[[i]],
                                               inits.gamma,
                                               gamma.waning.distribution)
  #record output and parameters
  op <- list(out = sim.out, pars = scenario.list.waning.hetero[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y%m%d"),"output",scenario.list.waning.hetero[[i]]$scenario,".RDS"))
  print(Sys.time()-start)
}






