#########################################################
#                                                        
#                  Simulate scenarios                               
#                                                        
#                  Author:Egil Fischer                              
#                  Contact:           e.a.j.fischer@uu.nl                   
#                  Creation date: 2-5-2023                         
#########################################################



#define number of types
itypes = 2;

# SCENARIOS #

#Layers#
#Baseline parameters for layer flock ####
param.list.baseline.layer <- list(
  scenario = "baseline_Layer", #scenario
  runs =10, #number of runs
  max.time = 17*30,#length of the run
  itypes = itypes, #type
  N0 = 45000, #population size
  initial= 10 , #initially infected
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2, #Variance infectious period
  transRate = matrix(c(0,0.012,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.95,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005 #per capita death rate #Mortality events
)
# 90% high titre parameters for layer flock ####
param.list.90percHighTitre.layer <- param.list.baseline.layer;
param.list.90percHighTitre.layer$scenario <- "HighTitre90perc_layer";
param.list.90percHighTitre.layer$p.hightitre <- 0.9;


# 90% high titre parameters no waning for broiler flock ####
param.list.90percHighTitreNoWaning.layer <- param.list.baseline.layer;
param.list.90percHighTitreNoWaning.layer$scenario <- "HighTitre90percNoWaning_layer";
param.list.90percHighTitreNoWaning.layer$p.hightitre <- 0.9;
param.list.90percHighTitreNoWaning.layer$transRate <-transRate = matrix(c(0,0.0,0.0,0), nrow = itypes);

#Broilers#
#Baseline parameters for broiler flock ####
param.list.baseline.broiler <- list(
  scenario = "baseline_broiler", #scenario
  runs =10, #number of runs
  max.time = 42,#length of the run
  itypes = itypes, #type
  N0 = 75000, #population size
  initial= 10 , #initially infected
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2, #Variance infectious period
  transRate = matrix(c(0,0.012,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.95,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005 #per capita death rate #Mortality events
)

# 90% high titre parameters for broiler flock ####
param.list.90percHighTitre.broiler <- param.list.baseline.broiler;
param.list.90percHighTitre.broiler$scenario <- "HighTitre90perc_broiler";
param.list.90percHighTitre.broiler$p.hightitre <- 0.9;


# 90% high titre parameters no waning for broiler flock ####
param.list.90percHighTitreNoWaning.broiler <- param.list.baseline.broiler;
param.list.90percHighTitreNoWaning.broiler$scenario <- "HighTitre90percNoWaning_broiler";
param.list.90percHighTitreNoWaning.broiler$p.hightitre <- 0.9;
param.list.90percHighTitreNoWaning.broiler$transRate <- matrix(c(0,0.0,0.0,0), nrow = itypes);


#do simulations#
output.baseline.broiler <- simulate.multitypeSIR(param.list.baseline.broiler)
output.90percHighTitre.broiler <- simulate.multitypeSIR(param.list.90percHighTitre.broiler)
output.90percHighTitreNoWaning.broiler <- simulate.multitypeSIR(param.list.90percHighTitreNoWaning.broiler)

#visualize
plot.output(output.baseline.broiler,c("I.1","I.2","R.1","R.2"), "Broiler base line")
plot.output(output.baseline.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Broiler base line")

plot.output(output.90percHighTitre.broiler,c("I.1","I.2","R.1","R.2"), "Broiler 90% high titre")
plot.output(output.90percHighTitre.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Broiler 90% high titre")

plot.output(output.90percHighTitreNoWaning.broiler,c("I.1","I.2","R.1","R.2"), "Broiler 90% high titre / no waning")
plot.output(output.90percHighTitreNoWaning.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Broiler 90% high titre / no waning")

#detection times



