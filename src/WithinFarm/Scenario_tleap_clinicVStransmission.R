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

#baseline parameters
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

#scenarios clinic versus transmission ####
# Clinical protection ####
scenario.list.clinic.vs.trans <-list()

for(i in c(1:6)){
  for(j in c(1:6)){
    for(k in c(1:6)){
    
  param.list <- param.list.baseline.layer;
  param.list$runs <- 100;#number of runs per combination
  param.list$trans.mean <-  1000; #no waning
  param.list$p.hightitre <- c(0,0.1,0.2,0.3,0.4,0.5)[i] #this only affects the transmission
  param.list$pdie <- c(c(0.001, 0.005, 0.01, 0.05, 0.1,0.5)[j],c(0.001, 0.005, 0.01, 0.05, 0.1,0.5)[k]) #level of clinical protection is equal for both high and low titre birds
  param.list$scenario <- paste0("layerClinic1",param.list$pdie[1]*1000 ,"Clinic2",param.list$pdie[2]*1000,"Phigh",param.list$p.hightitre*100);
  scenario.list.clinic.vs.trans[[(i-1)*36+ (j-1)*6 +k]]<- param.list
    }
  }
}




for(i in  c(1:length(scenario.list.clinic.vs.trans))){
  print(scenario.list.clinic.vs.trans[[i]]$scenario);
  print(start<- Sys.time());
  inits.gamma <- with(scenario.list.clinic.vs.trans[[i]],{
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
  sim.out<- sim.multitypeSEIR_tleap_distwaning(scenario.list.clinic.vs.trans[[i]],
                                               inits.gamma,
                                               gamma.waning.distribution)
  #record output and parameters
  op <- list(out = sim.out, pars = scenario.list.clinic.vs.trans[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y%m%d"),"output",scenario.list.clinic.vs.trans[[i]]$scenario,".RDS"))
  print(Sys.time()-start)
}
