#########################################################
#                                                        
#                  Processing of simulated scenarios                               
#                                                        
#                  Author:        E.A.J. Fischer                       
#                  Contact:       e.a.j.fischer@uu.nl                       
#                  Creation date  12-5-2023                       
#########################################################



#source required 
source("loadLibraries.R") 
source("postprocesSimulations.R")

#load simulations
output.baseline.layer <- load.sims("./output/baseline_layer", interval = 0.1)$output
output.60percHighTitre.layer <- load.sims("./output/HighTitre60perc_layer", interval = 0.1)$output
output.60percHighTitreLowWaning.layer <- load.sims("./output/HighTitre60percLowWaning_layer", interval = 0.1)$output
output.60percHighTitreNoWaning.layer <- load.sims("./output/HighTitre60percNoWaning_layer", interval = 0.1)$output

#load parameters (assuming all parameter lists are equal within a folder)
pars.baseline.layer <- load.sims("./output/baseline_layer", params = TRUE)$pars[[1]]
pars.60percHighTitre.layer <- load.sims("./output/HighTitre60perc_layer", params = TRUE)$pars[[1]]
pars.60percHighTitreLowWaning.layer <- load.sims("./output/HighTitre60percLowWaning_layer", params = TRUE)$pars[[1]]
pars.60percHighTitreNoWaning.layer <- load.sims("./output/HighTitre60percNoWaning_layer", params = TRUE)$pars[[1]]


#visualize
plot.output(output.baseline.layer,c("I.1","I.2","R.1","R.2"), "Layer base line")
ggsave("./output/figures/baselinelayer.png")
plot.output(output.baseline.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer base line")
ggsave("./output/figures/baselinelayerdeaths.png")
gc()
plot.output(output.60percHighTitre.layer,c("I.1","I.2","R.1","R.2"), "Layer 60% high titre")
ggsave("./output/figures/HighTitre60preclayer.png")
plot.output(output.60percHighTitre.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 60% high titre")
ggsave("./output/figures/HighTitre60preclayerdeaths.png")
gc()
plot.output(output.60percHighTitreLowWaning.layer,c("I.1","I.2","R.1","R.2"), "Layer 60% high titre / low waning")
ggsave("./output/figures/HighTitre60precLowWaninglayer.png")
plot.output(output.60percHighTitreLowWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 60% high titre / low waning")
ggsave("./output/figures/HighTitre60precLowWaninglayerdeaths.png")
gc()
plot.output(output.60percHighTitreNoWaning.layer,c("I.1","I.2","R.1","R.2"), "Layer 60% high titre / no waning")
ggsave("./output/figures/HighTitre60preclayerNoWaning.png")
plot.output(output.60percHighTitreNoWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 60% high titre / no waning")
ggsave("./output/figures/HighTitre60preclayerNoWaningdeaths.png")
gc()

#surveillance
N0 <- pars.baseline.layer$N0
reps <- 10;
rm(surveillance.baseline.layer)
surveillance.baseline.layer <-detection.times.surveillance(output.baseline.layer, 
                                               c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                               time.interval = 1,
                                               ints = 2,
                                               threshold = 0.005*N0,
                                               c("DI.1","DI.2","DR.1","DR.2"), 
                                               pfarm = 1., 
                                               panimal = .001,
                                               se = 0.9,
                                               Detectables.incidence = TRUE,
                                               reps = reps);


det.times <- surveillance.baseline.layer%>%select("run","pas.det.time")
names(det.times)<-c("run","baseline")
ggplot(data = surveillance.baseline.layer )+
  geom_histogram(aes(pas.det.time),colour = "black", fill = "black", alpha = 0.5, binwidth = 1.0)+
  geom_histogram(aes(ac.det.time),colour = "black", fill = "orange",alpha = 0.5, binwidth = 1.0);
sum(surveillance.baseline.layer$ac.succes,na.rm = TRUE)/length(surveillance.baseline.layer$ac.succes)
sum(is.na(surveillance.baseline.layer$ac.succes))/length(surveillance.baseline.layer$ac.succes)

N0 <- pars.60percHighTitre.layer$N0
rm(surveillance.60percHighTitre.layer)
surveillance.60percHighTitre.layer <-detection.times.surveillance(output.60percHighTitre.layer, 
                                                           c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                           time.interval = 1,
                                                           ints = 2,
                                                           threshold = 0.005*N0,
                                                           c("DI.1","DI.2","DR.1","DR.2"), 
                                                           pfarm = 1., 
                                                           panimal = .001,
                                                           se = 0.9,
                                                           Detectables.incidence = TRUE,
                                                           reps = reps);



det.times <- cbind(det.times,surveillance.60percHighTitre.layer%>%select("pas.det.time"))
names(det.times)<-c("run","baseline","HighTitre60perc")
ggplot(data = surveillance.60percHighTitre.layer )+
  geom_histogram(aes(pas.det.time),colour = "black", fill = "black", alpha = 0.5, binwidth = 1.0)+
  geom_histogram(aes(ac.det.time),colour = "black", fill = "orange",alpha = 0.5, binwidth = 1.0)
sum(surveillance.60percHighTitre.layer$ac.succes,na.rm = TRUE)/length(surveillance.60percHighTitre.layer$ac.succes)
sum(is.na(surveillance.60percHighTitre.layer$ac.succes))/length(surveillance.60percHighTitre.layer$ac.succes)


N0 <- pars.60percHighTitreLowWaning.layer$N0
rm(surveillance.60percHighTitreLowWaning.layer)
surveillance.60percHighTitreLowWaning.layer <-detection.times.surveillance(output.60percHighTitreLowWaning.layer, 
                                                                  c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                                  time.interval = 1,
                                                                  ints = 2,
                                                                  threshold = 0.005*N0,
                                                                  c("DI.1","DI.2","DR.1","DR.2"), 
                                                                  pfarm = 1., 
                                                                  panimal = .001,
                                                                  se = 0.9,
                                                                  Detectables.incidence = TRUE,
                                                                  reps = reps);



det.times <- cbind(det.times,surveillance.60percHighTitreLowWaning.layer%>%select("pas.det.time"))
names(det.times)<-c("run","baseline","HighTitre60perc","HighTitre60percLowWaning")
ggplot(data = surveillance.60percHighTitreLowWaning.layer )+
  geom_histogram(aes(pas.det.time),colour = "black", fill = "black", alpha = 0.5, binwidth = 1.0)+
  geom_histogram(aes(ac.det.time),colour = "black", fill = "orange",alpha = 0.5, binwidth = 1.0);
sum(surveillance.60percHighTitreLowWaning.layer$ac.succes,na.rm = TRUE)/length(surveillance.60percHighTitreLowWaning.layer$ac.succes)
sum(is.na(surveillance.60percHighTitreLowWaning.layer$ac.succes))/length(surveillance.60percHighTitreLowWaning.layer$ac.succes)


N0 <- pars.60percHighTitreNoWaning.layer$N0
rm(surveillance.60percHighTitreNoWaning.layer)
surveillance.60percHighTitreNoWaning.layer <-detection.times.surveillance(output.60percHighTitreNoWaning.layer, 
                                                                           c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                                           time.interval = 1,
                                                                           ints = 2,
                                                                           threshold = 0.005*N0,
                                                                           c("DI.1","DI.2","DR.1","DR.2"), 
                                                                           pfarm = 1., 
                                                                           panimal = .001,
                                                                           se = 0.9,
                                                                           Detectables.incidence = TRUE,
                                                                           reps = reps);



det.times <- cbind(det.times,surveillance.60percHighTitreNoWaning.layer%>%select("pas.det.time"))
names(det.times)<-c("run","baseline","HighTitre60perc","HighTitre60percLowWaning","HighTitre60percNoWaning")
ggplot(data = surveillance.60percHighTitreNoWaning.layer )+
  geom_histogram(aes(pas.det.time),    colour = "black", 
                 fill = "black", 
                 alpha = 0.5, binwidth = 1.0)+
  geom_histogram(aes(ac.det.time),colour = "black", 
                 fill = "orange",
                 alpha = 0.5, binwidth = 1.0)
sum(surveillance.60percHighTitreNoWaning.layer$ac.succes,na.rm = TRUE)/length(surveillance.60percHighTitreNoWaning.layer$ac.succes)
sum(is.na(surveillance.60percHighTitreNoWaning.layer$ac.succes))/length(surveillance.60percHighTitreNoWaning.layer$ac.succes)


#human exposure
humanexposure <- human.exposure.total.multiple.runs(output.baseline.layer,pars.baseline.layer$beta,det.times$baseline)
humanexposure <- cbind(humanexposure, 
                       human.exposure.total.multiple.runs(output.60percHighTitre.layer,pars.60percHighTitre.layer$beta,det.times$HighTitre60perc))
humanexposure <- cbind(humanexposure, 
                       human.exposure.total.multiple.runs(output.60percHighTitreNoWaning.layer,pars.60percHighTitreNoWaning.layer$beta,det.times$HighTitre60percNoWaning))
humanexposure <- cbind(humanexposure, 
                       human.exposure.total.multiple.runs(output.60percHighTitreLowWaning.layer,pars.60percHighTitreLowWaning.layer$beta,det.times$HighTitre60percLowWaning))

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
