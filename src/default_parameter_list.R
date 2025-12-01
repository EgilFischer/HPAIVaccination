#default parameter list that will work with all function
default_parameter_list <- list(
  #Simulation settings####
  scenario="default", 
  runs = 100, #number of runs
  deltat = 0.005, #time step; nb: 0.01 is faster but not accurate, 0.001 is very slow and only mild improvement in accuracy
  #Flock settings####
  length.round = 3*30,#length of the round. In these simulations we will only look at the first three months after introduction
  age_at_d0 = 120, #unit: days; model starts at t = 0 when chickens are 120 days/17 wks old (important for egg production curve and vx buildup)
  N0 = 100, #population size
  #Vaccination settings ####
  itypes = 2, #types
  p.protect = 0., #proportion of the flock that is vaccinated
  #Introduction settings ####
  no.introduction= 1, #initially infected at introduction. Titre depends on the fraction in each of the titre groups.
  introduction.at.random = TRUE, #chance process based on number of susceptibles per itype
  intro.time = 0, #time during the round at which the infection is introduced
  #Infection settings ####
  beta = matrix(c(1.99, 1.99,0.23,0.23), ncol = 2), #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
  latency.period = c(1,1), #latency period of 1 day
  k.latency =2, #k parameter Erlang distribution
  infectious.period = c(4.28,3.59), #Duration infectious period 
  k.infectious = 2, #k parameter Erlang distribution, k=1 gives sim results comparable to pmaj predictions for exponential 
  trans.mean.wane = 450, #only considers transition from type 2 to 1; for waning off set to 20000
  trans.var.wane = 50, #only considers transition from type 2 to 1
  trans.mean.buildup = 12, #mean time to build-up immunity
  trans.var.buildup = 1, #variance of time to build-up immunity
  pdie = c(0.4, 0.01),#probability of dying at end of infectious period
  mortRate = 2e-4, #per capita death rate for layers (Hobbelen et al. 2020; DOI: https://doi.org/10.1038/s41598-020-68623-w)
  #Egg settings ####
  eh = 0.57, # daily egg-laying rate for healthy chickens (S, E, R)
  ei = 0.285,# daily egg-laying rate for infected chickens (I) (50% reduction in egg laying)
  disfigured = 0.1, # daily removal rate for infected eggs (10% do not pass inspection)
  pickup_time = 3 #how often eggs are picked up 
)