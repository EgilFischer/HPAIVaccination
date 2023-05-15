#########################################################
#                                                        
#                  Simulate scenarios                               
#                                                        
#                  Author:Egil Fischer                              
#                  Contact:           e.a.j.fischer@uu.nl                   
#                  Creation date: 2-5-2023                         
#########################################################

#load libraries
source("loadLibraries.R") 

#define number of types
itypes = 2;

# SCENARIOS #

#Layers#
#Baseline parameters for layer flock ####
param.list.baseline.layer <- list(
  scenario = "baseline_Layer", #scenario
  runs = 40, #number of runs
  max.time = 17*30,#length of the run
  itypes = itypes, #type
  N0 = 45000, #population size
  initial= 10 , #initially infected
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2, #Variance infectious period
  transRate = matrix(c(0,0.038,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.95,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005 #per capita death rate #Mortality events
)
# 60% high titre parameters for layer flock ####
param.list.60percHighTitre.layer <- param.list.baseline.layer;
param.list.60percHighTitre.layer$scenario <- "HighTitre60percLowWaning_layer";
param.list.60percHighTitre.layer$p.hightitre <- 0.6;


# 60% high titre parameters no waning for broiler flock ####
param.list.60percHighTitreNoWaning.layer <- param.list.baseline.layer;
param.list.60percHighTitreNoWaning.layer$scenario <- "HighTitre60percNoWaning_layer";
param.list.60percHighTitreNoWaning.layer$p.hightitre <- 0.6;
param.list.60percHighTitreNoWaning.layer$transRate <- matrix(c(0,0.0,0.0,0), nrow = itypes);


# 60% high titre parameters low waning for broiler flock ####
param.list.60percHighTitreLowWaning.layer <- param.list.baseline.layer;
param.list.60percHighTitreLowWaning.layer$scenario <- "HighTitre60percNoWaning_layer";
param.list.60percHighTitreLowWaning.layer$p.hightitre <- 0.6;
param.list.60percHighTitreLowWaning.layer$transRate <- matrix(c(0,0.12,0.0,0), nrow = itypes);



#probability of a major outbreak
q1q2.baseline.layer <- q1q2(with(param.list.baseline.layer,list(N0 = N0,
                                                 gamma = 1/infectious.period,
                                                 beta = beta,
                                                 mu = mortRate)))

n <- param.list.baseline.layer$initial
ggplot(data =q1q2.baseline.layer)+
  geom_path(aes(p,q1^(n), colour = "Low Titre"))+
  geom_path(aes(p,q2^(n), colour = "High Titre"))+
  geom_path(aes(p,(q1^(n*(1-p))* q2^(n*p)), colour = "Low & High titre ratio"))+
  geom_path(aes(p,Rv, colour = "R"))+
  xlab("Proportion with high titre")+
  ylab("Probability of a minor outbreak")+ 
  scale_colour_manual(name = paste("Introduction by",n,"birds \n"),
                      #labels = c("Low Titre","Type 1 & 2","Type 2 only"),
                      values = c("red","blue","grey","black"))#+
  #theme(legend.position="bottom")
ggsave("./output/figures/probMinorOutbreakLayer.png" )

#do simulations#
output.baseline.layer <- simulate.multitypeSIR(param.list.baseline.layer)
output.60percHighTitre.layer <- simulate.multitypeSIR(param.list.60percHighTitre.layer)
output.60percHighTitreNoWaning.layer <- simulate.multitypeSIR(param.list.60percHighTitreNoWaning.layer)

#load simulations
output.baseline.layer <- load.sims("./output/baseline_layer", interval = 0.1)$output
output.60percHighTitre.layer <- load.sims("./output/HighTitre60perc_layer", interval = 0.1)$output
output.60percHighTitreLowWaning.layer <- load.sims("./output/HighTitre60percLowWaning_layer", interval = 0.1)$output
output.60percHighTitreNoWaning.layer <- load.sims("./output/HighTitre60percNoWaning_layer", interval = 0.1)$output

#remove error in runs
output.baseline.layer[which(output.baseline.layer$time == 0 & output.baseline.layer$N == 45000, arr.ind = TRUE),"run"]<- c(1:41)
output.baseline.layer <- output.baseline.layer%>%filter(run!=41);
output.60percHighTitre.layer[which(output.60percHighTitre.layer$time == 0 & output.60percHighTitre.layer$N == 45000, arr.ind = TRUE),"run"]<- c(1:10)
output.60percHighTitreLowWaning.layer[which(output.60percHighTitreLowWaning.layer$time == 0 & output.60percHighTitreLowWaning.layer$N == 45000, arr.ind = TRUE),"run"]<- c(1:40)
output.60percHighTitreNoWaning.layer[which(output.60percHighTitreNoWaning.layer$time == 0 & output.60percHighTitreNoWaning.layer$N == 45000, arr.ind = TRUE),"run"]<- c(1:11)
output.60percHighTitreNoWaning.layer <- output.60percHighTitreNoWaning.layer%>%filter(run!=11);

# #visualize
# plot.output(output.baseline.layer,c("I.1","I.2","R.1","R.2"), "Layer base line")
# ggsave("./output/figures/baselinelayer.png")
# plot.output(output.baseline.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer base line")
# ggsave("./output/figures/baselinelayerdeaths.png")
# gc()
# plot.output(output.60percHighTitre.layer,c("I.1","I.2","R.1","R.2"), "Layer 60% high titre")
# ggsave("./output/figures/HighTitre60preclayer.png")
# plot.output(output.60percHighTitre.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 60% high titre")
# ggsave("./output/figures/HighTitre60preclayerdeaths.png")
# gc()
# plot.output(output.60percHighTitreLowWaning.layer,c("I.1","I.2","R.1","R.2"), "Layer 60% high titre / high waning")
# ggsave("./output/figures/HighTitre60precLowWaninglayer.png")
# plot.output(output.60percHighTitreLowWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 60% high titre / high waning")
# ggsave("./output/figures/HighTitre60precLowWaninglayerdeaths.png")
# gc()
# plot.output(output.60percHighTitreNoWaning.layer,c("I.1","I.2","R.1","R.2"), "Layer 60% high titre / no waning")
# ggsave("./output/figures/HighTitre60preclayerNoWaning.png")
# plot.output(output.60percHighTitreNoWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 60% high titre / no waning")
# ggsave("./output/figures/HighTitre60preclayerNoWaningdeaths.png")
# gc()

#detection times by passive surveillance
det.times <- detection.times(output.baseline.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2") ,detection.time.threshold.subsequent,1,0.005*45000, n= 2)
det.times <- cbind(det.times,detection.times(output.60percHighTitre.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2") ,detection.time.threshold.subsequent,1,0.005*45000, n= 2)$det.time)
det.times <- cbind(det.times,detection.times(output.60percHighTitreNoWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2") ,detection.time.threshold.subsequent,1,0.005*45000, n= 2)$det.time)
det.times <- cbind(det.times,detection.times(output.60percHighTitreLowWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2") ,detection.time.threshold.subsequent,1,0.005*45000, n= 2)$det.time)
names(det.times)<- c("run","baseline","60HighTitre","60HighTitreNoWaning","60HighTitreLowWaning")
det.times

#surveillance
surveillance.baseline.layer <-detection.times.surveillance(output.baseline.layer, c("DR.1","DR.2"),time.interval = 1,ints = 2,threshold = 0.005*45000,c("DR.1","DR.2"), pfarm = 1., panimal = .001,se = 0.9,Detectables.incidence = TRUE);
for(i in c(2:10)){
  surveillance.baseline.layer <- rbind(surveillance.baseline.layer,detection.times.surveillance(output.baseline.layer, c("DR.1","DR.2"),time.interval = 1,ints = 2,threshold = 0.005*45000,c("DR.1","DR.2"), pfarm = 0.1, panimal = .001,se = 0.95))
}
ggplot(surveillance.baseline.layer)+geom_histogram(aes(ac.det.time))+geom_histogram(aes(pas.det.time, fill = 'red'))
sum(surveillance.baseline.layer$ac.succes)/length(surveillance.baseline.layer$ac.succes)

surveillance.60percHighTitreNoWaning.layer <-detection.times.surveillance(output.60percHighTitreNoWaning.layer, c("DR.1","DR.2"),time.interval = 1,ints = 2,threshold = 0.005*45000,c("DR.1","DR.2"), pfarm = 0.1, panimal = .0001,se = 0.95);
for(i in c(2:10)){
  surveillance.60percHighTitreNoWaning.layer <- rbind(surveillance.60percHighTitreNoWaning.layer,detection.times.surveillance(output.60percHighTitreNoWaning.layer, c("DR.1","DR.2"),time.interval = 1,ints = 2,threshold = 0.005*45000,c("DR.1","DR.2"), pfarm = 0.1, panimal = .001,se = 0.9))
}
ggplot(surveillance.60percHighTitreNoWaning.layer)+geom_histogram(aes(ac.det.time))+geom_histogram(aes(pas.det.time, fill = 'red'))
sum(surveillance.60percHighTitreNoWaning.layer$ac.succes)/length(surveillance.60percHighTitreNoWaning.layer$ac.succes)

#human exposure
humanexposure <- human.exposure.total.multiple.runs(output.baseline.layer,pars.list.baseline.layer$beta,det.times$baseline)
humanexposure <- cbind(humanexposure, 
                       human.exposure.total.multiple.runs(output.60percHighTitre.layer,pars.list.60percHighTitre.layer$beta,det.times$'60HighTitre'))
humanexposure <- cbind(humanexposure, 
                       human.exposure.total.multiple.runs(output.60percHighTitreNoWaning.layer,pars.list.60percHighTitreNoWaning.layer$beta,det.times$'60HighTitreNoWaning'))
humanexposure <- cbind(humanexposure, 
                       human.exposure.total.multiple.runs(output.60percHighTitreLowWaning.layer,pars.list.60percHighTitreLowWaning.layer$beta,det.times$'60HighTitreLowWaning'))

names(humanexposure)<- c("run","baseline_total","baseline_detection",
                         "run","HighTitre60perc_total","HighTitre60perc_detection",
                         "run","HighTitre60percNoWaning_total","HighTitre60percNoWaning_detection",
                         "run","HighTitre60percLowWaning_total","HighTitre60percLowWaning_detection")

humanexposure.detection.ratio <- rbind(data.frame(scenario = "baseline",
  ratio = reshape2::melt(outer(humanexposure$baseline_detection,humanexposure$baseline_detection,"/"))$value),
  data.frame(scenario = "HighTitre60perc",
  ratio = reshape2::melt(outer(humanexposure$HighTitre60perc_detection,humanexposure$baseline_detection,"/"))$value),
  data.frame(scenario = "HighTitre60percNoWaning",
  ratio = reshape2::melt(outer(humanexposure$HighTitre60percNoWaning_detection,humanexposure$baseline_detection,"/"))$value),
  data.frame(scenario = "HighTitre60percLowWaning",
  ratio = reshape2::melt(outer(humanexposure$HighTitre60percLowWaning_detection,humanexposure$baseline_detection,"/"))$value))

ggplot(humanexposure.detection.ratio) +
  geom_histogram(aes(ratio,after_stat(density), fill = scenario))+
  xlab("Risk ratio of exposure (reference baseline)") + geom_vline(xintercept = 1)+facet_grid(scenario~.)

#prepare for plotting
plot.humanexposure <- reshape2::melt(humanexposure, id.vars = c("run"))
plot.humanexposure$scenario <- sapply(str_split(plot.humanexposure$variable, pattern = "_"),"[[",1)
plot.humanexposure$totdet <- sapply(str_split(plot.humanexposure$variable, pattern = "_"),"[[",2)

ggplot(data = plot.humanexposure)+geom_histogram(aes(value, fill = scenario ))+
  facet_grid(scenario~totdet)+ggtitle("Cumulative human exposure Layers")

#Broilers####
#Baseline parameters for broiler flock ####
param.list.baseline.broiler <- list(
  scenario = "baseline_broiler", #scenario
  runs = 40, #number of runs
  max.time = 42,#length of the run
  itypes = itypes, #type
  N0 = 75000, #population size
  initial= 10 , #initially infected
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2, #Variance infectious period
  transRate = matrix(c(0,0.038,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.95,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005 #per capita death rate #Mortality events
)

# 60% high titre parameters for broiler flock ####
param.list.60percHighTitre.broiler <- param.list.baseline.broiler;
param.list.60percHighTitre.broiler$scenario <- "HighTitre60percLowWaning_broiler";
param.list.60percHighTitre.broiler$p.hightitre <- 0.6;


# 60% high titre parameters no waning for broiler flock ####
param.list.60percHighTitreNoWaning.broiler <- param.list.baseline.broiler;
param.list.60percHighTitreNoWaning.broiler$scenario <- "HighTitre60percNoWaning_broiler";
param.list.60percHighTitreNoWaning.broiler$p.hightitre <- 0.6;
param.list.60percHighTitreNoWaning.broiler$transRate <- matrix(c(0,0.0,0.0,0), nrow = itypes);


#do simulations#
output.baseline.broiler <- simulate.multitypeSIR(param.list.baseline.broiler)
output.60percHighTitre.broiler <- simulate.multitypeSIR(param.list.60percHighTitre.broiler)
output.60percHighTitreNoWaning.broiler <- simulate.multitypeSIR(param.list.60percHighTitreNoWaning.broiler)

#visualize
plot.output.sparse(output.baseline.broiler,c("I.1","I.2","R.1","R.2"), "Broiler base line",.2)
ggsave("./output/figures/baselinebroiler.png")
plot.output.sparse(output.baseline.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Broiler base line",.2)
ggsave("./output/figures/baselinebroilerdeaths.png")

plot.output.sparse(output.60percHighTitre.broiler,c("I.1","I.2","R.1","R.2"), "Broiler 60% high titre",0.1)
ggsave("./output/figures/HighTitre60precbroiler.png")
plot.output.sparse(output.60percHighTitre.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Broiler 60% high titre",0.1)
ggsave("./output/figures/HighTitre60precbroilerdeaths.png")

plot.output.sparse(output.60percHighTitreNoWaning.broiler,c("I.1","I.2","R.1","R.2"), "Broiler 60% high titre / no waning",.1)
ggsave("./output/figures/HighTitre60precbroilerNoWaning.png")
plot.output(output.60percHighTitreNoWaning.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Broiler 60% high titre / no waning")
ggsave("./output/figures/HighTitre60precbroilerNoWaningdeaths.png")

#detection times
det.times <- detection.times(output.baseline.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2") ,detection.time.threshold.subsequent,1,0.005*75000, n= 2)
det.times <- cbind(det.times,detection.times(output.60percHighTitre.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2") ,detection.time.threshold.subsequent,1,0.005*75000, n= 2)$det.time)
det.times <- cbind(det.times,detection.times(output.60percHighTitreNoWaning.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2") ,detection.time.threshold.subsequent,1,0.005*75000, n= 2)$det.time)
names(det.times)<- c("run","baseline","60HighTitre","60HighTitreNoWaning")
det.times
