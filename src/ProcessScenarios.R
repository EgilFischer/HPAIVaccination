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
output.90percHighTitre.layer <- load.sims("./output/HighTitre60perc_layer", interval = 0.1)$output
output.60percHighTitreHighWaning.layer <- load.sims("./output/HighTitre60percHighWaning_layer", interval = 0.1)$output
output.60percHighTitreNoWaning.layer <- load.sims("./output/HighTitre60percNoWaning_layer", interval = 0.1)$output

#load parameters (assuming all parameter lists are equal within a folder)
pars.baseline.layer <- load.sims("./output/baseline_layer", params = TRUE)$pars[[1]]
pars.60percHighTitre.layer <- load.sims("./output/HighTitre60perc_layer", params = TRUE)$pars[[1]]
pars.60percHighTitreHighWaning.layer <- load.sims("./output/HighTitre60percHighWaning_layer", params = TRUE)$pars[[1]]
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
plot.output(output.60percHighTitreHighWaning.layer,c("I.1","I.2","R.1","R.2"), "Layer 60% high titre / High waning")
ggsave("./output/figures/HighTitre60precHighWaninglayer.png")
plot.output(output.60percHighTitreHighWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 60% high titre / High waning")
ggsave("./output/figures/HighTitre60precHighWaninglayerdeaths.png")
gc()
plot.output(output.60percHighTitreNoWaning.layer,c("I.1","I.2","R.1","R.2"), "Layer 60% high titre / no waning")
ggsave("./output/figures/HighTitre60preclayerNoWaning.png")
plot.output(output.60percHighTitreNoWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 60% high titre / no waning")
ggsave("./output/figures/HighTitre60preclayerNoWaningdeaths.png")
gc()

#surveillance
N0 <- pars.baseline.layer$N0
reps <- 10;
rm(surveillance.layer)
surveillance.layer <-cbind(detection.times.surveillance(output.baseline.layer, 
                                               c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                               time.interval = 1,
                                               ints = 2,
                                               threshold = 0.005*N0,
                                               c("DI.1","DI.2","DR.1","DR.2"), 
                                               pfarm = 1., 
                                               panimal = .001,
                                               se = 0.9,
                                               Detectables.incidence = TRUE,
                                               reps = reps),
                                    scenario = "baseline_layer");

N0 <- pars.60percHighTitre.layer$N0
surveillance.layer <-rbind(surveillance.layer, cbind(detection.times.surveillance(output.60percHighTitre.layer, 
                                                           c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                           time.interval = 1,
                                                           ints = 2,
                                                           threshold = 0.005*N0,
                                                           c("DI.1","DI.2","DR.1","DR.2"), 
                                                           pfarm = 1., 
                                                           panimal = .001,
                                                           se = 0.9,
                                                           Detectables.incidence = TRUE,
                                                           reps = reps),
                                                     scenario = "HighTitre60perc_layer"));


N0 <- pars.60percHighTitreHighWaning.layer$N0
surveillance.layer <-rbind(surveillance.layer, cbind(detection.times.surveillance(output.60percHighTitreHighWaning.layer, 
                                                                  c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                                  time.interval = 1,
                                                                  ints = 2,
                                                                  threshold = 0.005*N0,
                                                                  c("DI.1","DI.2","DR.1","DR.2"), 
                                                                  pfarm = 1., 
                                                                  panimal = .001,
                                                                  se = 0.9,
                                                                  Detectables.incidence = TRUE,
                                                                  reps = reps),
                                                     scenario = "HighTitre60percHighWaning"));




N0 <- pars.60percHighTitreNoWaning.layer$N0

surveillance.layer <-rbind(surveillance.layer, cbind(detection.times.surveillance(output.60percHighTitreNoWaning.layer, 
                                                                           c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                                           time.interval = 1,
                                                                           ints = 2,
                                                                           threshold = 0.005*N0,
                                                                           c("DI.1","DI.2","DR.1","DR.2"), 
                                                                           pfarm = 1., 
                                                                           panimal = .001,
                                                                           se = 0.9,
                                                                           Detectables.incidence = TRUE,
                                                                           reps = reps),  
                                                     scenario = "HighTitre60percNoWaning"));


scenario.label <- list("baseline_layer" = "Baseline",
                       "HighTitre60perc_layer" = "60% - Slow",
                       "HighTitre60percHighWaning" = "60% - Fast",
                       "HighTitre60percNoWaning" = "60% - No")


ggplot(data = surveillance.layer%>%filter(scenario != "baseline_layer") )+
  geom_histogram(aes(x = pas.det.time, y=..count../sum(..count..)),    colour = "black", 
                  fill = "black", 
                  alpha = 0.5, binwidth = 1.0)+
   geom_histogram(aes(ac.det.time, y=..count../sum(..count..)),colour = "black", 
                  fill = "orange",
                  alpha = 0.5, binwidth = 1.0)+
  geom_density(aes(min.det.time,after_stat(count/length(min.det.times$run))),
                colour = "red",
                linewidth = 1,
               bw = 1)+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  ggtitle("Active and passive surveillance")+
  facet_grid(scenario~., labeller =  function(variable,value){
    return(scenario.label[value])
  })
ggsave("./output/figures/scenarios_layer_surveillance.png")


#human exposure
rm(humanexposure.min)
sim.runs <- length(unique(output.baseline.layer$run))
#as we used multiple repeats to determine active surveillance repeat
humanexposure.min <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,pars.baseline.layer$beta,surveillance.layer%>%filter(scenario == "baseline_layer"), var = "min.det.time", rep = 10), scenario = "baseline")
humanexposure.min <- rbind(humanexposure.min,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitre.layer,
                                                                    pars.60percHighTitre.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60perc_layer"), var = "min.det.time", rep = 10), scenario = "HighTitre60perc"))
humanexposure.min <- rbind(humanexposure.min,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitreHighWaning.layer,
                                                                    pars.60percHighTitreHighWaning.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60percHighWaning"), var = "min.det.time", rep = 10), scenario = "HighTitre60percHighWaning"))
humanexposure.min <- rbind(humanexposure.min,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitreNoWaning.layer,
                                                                    pars.60percHighTitreNoWaning.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60percNoWaning"), var = "min.det.time", rep = 10), scenario = "HighTitre60percNoWaning"))


baseline.average.det <- mean(unlist(humanexposure.min%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.min%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.min$ratio.det <- humanexposure.min$detection.exposure/baseline.average.det
humanexposure.min$ratio.tot <- humanexposure.min$total.exposure/baseline.average.tot


ggplot(humanexposure.min)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))
scenario.label <- list("baseline" = "Baseline",
                       "HighTitre60perc" = "60% - Slow",
                       "HighTitre60percHighWaning" = "60% - Fast",
                       "HighTitre60percNoWaning" = "60% - No")
ggplot(humanexposure.min) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~., labeller =  function(variable, value){
    return(scenario.label[value])
  }) +ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer.png")

#average ratios
humanexposure.min%>%reframe(.by = scenario,
  mean = mean(ratio.det),
  min = min(ratio.det),
  max = max(ratio.det),
  perc25 = quantile(ratio.det,0.25),
  perc75 =quantile(ratio.det,0.75)
  )


#passive detection
humanexposure.pas <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,pars.baseline.layer$beta,surveillance.layer%>%filter(scenario == "baseline_layer"), var = "pas.det.time", rep = 1), scenario = "baseline")
humanexposure.pas <- rbind(humanexposure.pas,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitre.layer,
                                                                    pars.60percHighTitre.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60perc_layer"), var = "pas.det.time", rep = 1), scenario = "HighTitre60perc"))
humanexposure.pas <- rbind(humanexposure.pas,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitreHighWaning.layer,
                                                                    pars.60percHighTitreHighWaning.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60percHighWaning"), var = "pas.det.time", rep = 1), scenario = "HighTitre60percHighWaning"))
humanexposure.pas <- rbind(humanexposure.pas,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitreNoWaning.layer,
                                                                    pars.60percHighTitreNoWaning.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60percNoWaning"), var = "pas.det.time", rep = 1), scenario = "HighTitre60percNoWaning"))


baseline.average.det <- mean(unlist(humanexposure.pas%>%filter(scenario = "baseline_layer")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.pas%>%filter(scenario = "baseline_layer")%>%select("total.exposure")))
humanexposure.pas$ratio.det <- humanexposure.pas$detection.exposure/baseline.average.det
humanexposure.pas$ratio.tot <- humanexposure.pas$total.exposure/baseline.average.tot


ggplot(humanexposure.pas)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))
scenario.label <- list("baseline" = "Baseline",
                       "HighTitre60perc" = "60% - Slow",
                       "HighTitre60percHighWaning" = "60% - Fast",
                       "HighTitre60percNoWaning" = "60% - No")
ggplot(humanexposure.pas) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~., labeller =  function(variable, value){
    return(scenario.label[value])
  }) +ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer.pas.png")


#active detection
humanexposure.ac <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,pars.baseline.layer$beta,surveillance.layer%>%filter(scenario == "baseline_layer"), var = "ac.det.time", rep = 10), scenario = "baseline")
humanexposure.ac <- rbind(humanexposure.ac,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitre.layer,
                                                                    pars.60percHighTitre.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60perc_layer"), var = "ac.det.time", rep = 10), scenario = "HighTitre60perc"))
humanexposure.ac <- rbind(humanexposure.ac,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitreHighWaning.layer,
                                                                    pars.60percHighTitreHighWaning.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60percHighWaning"), var = "ac.det.time", rep = 10), scenario = "HighTitre60percHighWaning"))
humanexposure.ac <- rbind(humanexposure.ac,
                           cbind(human.exposure.total.multiple.runs(output.60percHighTitreNoWaning.layer,
                                                                    pars.60percHighTitreNoWaning.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre60percNoWaning"), var = "ac.det.time", rep = 10), scenario = "HighTitre60percNoWaning"))


baseline.average.det <- mean(unlist(humanexposure.ac%>%filter(scenario = "baseline_layer")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.ac%>%filter(scenario = "baseline_layer")%>%select("total.exposure")))
humanexposure.ac$ratio.det <- humanexposure.ac$detection.exposure/baseline.average.det
humanexposure.ac$ratio.tot <- humanexposure.ac$total.exposure/baseline.average.tot


ggplot(humanexposure.ac)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))
scenario.label <- list("baseline" = "Baseline",
                       "HighTitre60perc" = "60% - Slow",
                       "HighTitre60percHighWaning" = "60% - Fast",
                       "HighTitre60percNoWaning" = "60% - No")
ggplot(humanexposure.ac) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~., labeller =  function(variable, value){
    return(scenario.label[value])
  }) +ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer.ac.png")

# 90% high titre ####

#load simulations

output.90percHighTitre.layer <- load.sims("./output/HighTitre90perc_layer", interval = 0.1)$output
#output.60percHighTitreHighWaning.layer <- load.sims("./output/HighTitre60percHighWaning_layer", interval = 0.1)$output
output.90percHighTitreNoWaning.layer <- load.sims("./output/HighTitre90percNoWaning_layer", interval = 0.1)$output

#load parameters (assuming all parameter lists are equal within a folder)
pars.90percHighTitre.layer <- load.sims("./output/HighTitre90perc_layer", params = TRUE)$pars[[1]]
#pars.90percHighTitreHighWaning.layer <- load.sims("./output/HighTitre90percHighWaning_layer", params = TRUE)$pars[[1]]
pars.90percHighTitreNoWaning.layer <- load.sims("./output/HighTitre90percNoWaning_layer", params = TRUE)$pars[[1]]


#visualize
plot.output(output.90percHighTitre.layer,c("I.1","I.2","R.1","R.2"), "Layer 90% high titre")
ggsave("./output/figures/HighTitre90preclayer.png")
plot.output(output.90percHighTitre.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 90% high titre")
ggsave("./output/figures/HighTitre90preclayerdeaths.png")
gc()
# plot.output(output.90percHighTitreHighWaning.layer,c("I.1","I.2","R.1","R.2"), "Layer 90% high titre / High waning")
# ggsave("./output/figures/HighTitre90precHighWaninglayer.png")
# plot.output(output.90percHighTitreHighWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 90% high titre / High waning")
# ggsave("./output/figures/HighTitre90precHighWaninglayerdeaths.png")
# gc()
plot.output(output.90percHighTitreNoWaning.layer,c("I.1","I.2","R.1","R.2"), "Layer 90% high titre / no waning")
ggsave("./output/figures/HighTitre90preclayerNoWaning.png")
plot.output(output.90percHighTitreNoWaning.layer,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Layer 90% high titre / no waning")
ggsave("./output/figures/HighTitre90preclayerNoWaningdeaths.png")
gc()

#surveillance
N0 <- pars.baseline.layer$N0
reps <- 10;
rm(surveillance.layer)
surveillance.layer <-cbind(detection.times.surveillance(output.baseline.layer, 
                                                        c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                        time.interval = 1,
                                                        ints = 2,
                                                        threshold = 0.005*N0,
                                                        c("DI.1","DI.2","DR.1","DR.2"), 
                                                        pfarm = 1., 
                                                        panimal = .001,
                                                        se = 0.9,
                                                        Detectables.incidence = TRUE,
                                                        reps = reps),
                           scenario = "baseline_layer");

N0 <- pars.90percHighTitre.layer$N0
surveillance.layer <-rbind(surveillance.layer, cbind(detection.times.surveillance(output.90percHighTitre.layer, 
                                                                                  c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                                                  time.interval = 1,
                                                                                  ints = 2,
                                                                                  threshold = 0.005*N0,
                                                                                  c("DI.1","DI.2","DR.1","DR.2"), 
                                                                                  pfarm = 1., 
                                                                                  panimal = .001,
                                                                                  se = 0.9,
                                                                                  Detectables.incidence = TRUE,
                                                                                  reps = reps),
                                                     scenario = "HighTitre90perc_layer"));


# N0 <- pars.90percHighTitreHighWaning.layer$N0
# surveillance.layer <-rbind(surveillance.layer, cbind(detection.times.surveillance(output.90percHighTitreHighWaning.layer, 
#                                                                                   c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
#                                                                                   time.interval = 1,
#                                                                                   ints = 2,
#                                                                                   threshold = 0.005*N0,
#                                                                                   c("DI.1","DI.2","DR.1","DR.2"), 
#                                                                                   pfarm = 1., 
#                                                                                   panimal = .001,
#                                                                                   se = 0.9,
#                                                                                   Detectables.incidence = TRUE,
#                                                                                   reps = reps),
#                                                      scenario = "HighTitre90percHighWaning"));
# 



N0 <- pars.90percHighTitreNoWaning.layer$N0

surveillance.layer <-rbind(surveillance.layer, cbind(detection.times.surveillance(output.90percHighTitreNoWaning.layer, 
                                                                                  c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                                                  time.interval = 1,
                                                                                  ints = 2,
                                                                                  threshold = 0.005*N0,
                                                                                  c("DI.1","DI.2","DR.1","DR.2"), 
                                                                                  pfarm = 1., 
                                                                                  panimal = .001,
                                                                                  se = 0.9,
                                                                                  Detectables.incidence = TRUE,
                                                                                  reps = reps),  
                                                     scenario = "HighTitre90percNoWaning"));


scenario.label <- list("baseline_layer" = "Baseline",
                       "HighTitre90perc_layer" = "90% - Slow",
#                       "HighTitre90percHighWaning" = "90% - Fast",
                       "HighTitre90percNoWaning" = "90% - No")


ggplot(data = surveillance.layer%>%filter(scenario != "baseline_layer") )+
  geom_histogram(aes(x = pas.det.time, y=..count../sum(..count..)),    colour = "black", 
                 fill = "black", 
                 alpha = 0.5, binwidth = 1.0)+
  geom_histogram(aes(ac.det.time, y=..count../sum(..count..)),colour = "black", 
                 fill = "orange",
                 alpha = 0.5, binwidth = 1.0)+
  geom_density(aes(min.det.time,after_stat(count/length(min.det.times$run))),
               colour = "red",
               linewidth = 1,
               bw = 1)+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  ggtitle("Active and passive surveillance")+
  facet_grid(scenario~., labeller =  function(variable,value){
    return(scenario.label[value])
  })
ggsave("./output/figures/scenarios_layer_surveillance90perc.png")


#human exposure
rm(humanexposure.min)
sim.runs <- length(unique(output.baseline.layer$run))
#as we used multiple repeats to determine active surveillance repeat
humanexposure.min <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,pars.baseline.layer$beta,surveillance.layer%>%filter(scenario == "baseline_layer"), var = "min.det.time", rep = 10), scenario = "baseline")
humanexposure.min <- rbind(humanexposure.min,
                           cbind(human.exposure.total.multiple.runs(output.90percHighTitre.layer,
                                                                    pars.90percHighTitre.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre90perc_layer"), var = "min.det.time", rep = 10), scenario = "HighTitre90perc"))
# humanexposure.min <- rbind(humanexposure.min,
#                            cbind(human.exposure.total.multiple.runs(output.90percHighTitreHighWaning.layer,
#                                                                     pars.90percHighTitreHighWaning.layer$beta,
#                                                                     surveillance.layer%>%filter(scenario == "HighTitre90percHighWaning"), var = "min.det.time", rep = 10), scenario = "HighTitre90percHighWaning"))
humanexposure.min <- rbind(humanexposure.min,
                           cbind(human.exposure.total.multiple.runs(output.90percHighTitreNoWaning.layer,
                                                                    pars.90percHighTitreNoWaning.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre90percNoWaning"), var = "min.det.time", rep = 10), scenario = "HighTitre90percNoWaning"))


baseline.average.det <- mean(unlist(humanexposure.min%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.min%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.min$ratio.det <- humanexposure.min$detection.exposure/baseline.average.det
humanexposure.min$ratio.tot <- humanexposure.min$total.exposure/baseline.average.tot


ggplot(humanexposure.min)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))
scenario.label <- list("baseline" = "Baseline",
                       "HighTitre90perc" = "90% - Slow",
  #                     "HighTitre90percHighWaning" = "90% - Fast",
                       "HighTitre90percNoWaning" = "90% - No")
ggplot(humanexposure.min) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~., labeller =  function(variable, value){
    return(scenario.label[value])
  }) +ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer90perc.png")

#average ratios
humanexposure.min%>%reframe(.by = scenario,
                            mean = mean(ratio.det),
                            min = min(ratio.det),
                            max = max(ratio.det),
                            perc25 = quantile(ratio.det,0.25),
                            perc75 =quantile(ratio.det,0.75)
)


#passive detection
humanexposure.pas <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,pars.baseline.layer$beta,surveillance.layer%>%filter(scenario == "baseline_layer"), var = "pas.det.time", rep = 1), scenario = "baseline")
humanexposure.pas <- rbind(humanexposure.pas,
                           cbind(human.exposure.total.multiple.runs(output.90percHighTitre.layer,
                                                                    pars.90percHighTitre.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre90perc_layer"), var = "pas.det.time", rep = 1), scenario = "HighTitre90perc"))
# humanexposure.pas <- rbind(humanexposure.pas,
#                            cbind(human.exposure.total.multiple.runs(output.90percHighTitreHighWaning.layer,
#                                                                     pars.90percHighTitreHighWaning.layer$beta,
#                                                                     surveillance.layer%>%filter(scenario == "HighTitre90percHighWaning"), var = "pas.det.time", rep = 1), scenario = "HighTitre90percHighWaning"))
humanexposure.pas <- rbind(humanexposure.pas,
                           cbind(human.exposure.total.multiple.runs(output.90percHighTitreNoWaning.layer,
                                                                    pars.90percHighTitreNoWaning.layer$beta,
                                                                    surveillance.layer%>%filter(scenario == "HighTitre90percNoWaning"), var = "pas.det.time", rep = 1), scenario = "HighTitre90percNoWaning"))


baseline.average.det <- mean(unlist(humanexposure.pas%>%filter(scenario = "baseline_layer")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.pas%>%filter(scenario = "baseline_layer")%>%select("total.exposure")))
humanexposure.pas$ratio.det <- humanexposure.pas$detection.exposure/baseline.average.det
humanexposure.pas$ratio.tot <- humanexposure.pas$total.exposure/baseline.average.tot


ggplot(humanexposure.pas)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))
scenario.label <- list("baseline" = "Baseline",
                       "HighTitre90perc" = "90% - Slow",
#                       "HighTitre90percHighWaning" = "90% - Fast",
                       "HighTitre90percNoWaning" = "90% - No")
ggplot(humanexposure.pas) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~., labeller =  function(variable, value){
    return(scenario.label[value])
  }) +ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer90perc.pas.png")


#active detection
humanexposure.ac <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,pars.baseline.layer$beta,surveillance.layer%>%filter(scenario == "baseline_layer"), var = "ac.det.time", rep = 10), scenario = "baseline")
humanexposure.ac <- rbind(humanexposure.ac,
                          cbind(human.exposure.total.multiple.runs(output.90percHighTitre.layer,
                                                                   pars.90percHighTitre.layer$beta,
                                                                   surveillance.layer%>%filter(scenario == "HighTitre90perc_layer"), var = "ac.det.time", rep = 10), scenario = "HighTitre90perc"))
# humanexposure.ac <- rbind(humanexposure.ac,
#                           cbind(human.exposure.total.multiple.runs(output.90percHighTitreHighWaning.layer,
#                                                                    pars.90percHighTitreHighWaning.layer$beta,
#                                                                    surveillance.layer%>%filter(scenario == "HighTitre90percHighWaning"), var = "ac.det.time", rep = 10), scenario = "HighTitre90percHighWaning"))
humanexposure.ac <- rbind(humanexposure.ac,
                          cbind(human.exposure.total.multiple.runs(output.90percHighTitreNoWaning.layer,
                                                                   pars.90percHighTitreNoWaning.layer$beta,
                                                                   surveillance.layer%>%filter(scenario == "HighTitre90percNoWaning"), var = "ac.det.time", rep = 10), scenario = "HighTitre90percNoWaning"))


baseline.average.det <- mean(unlist(humanexposure.ac%>%filter(scenario = "baseline_layer")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.ac%>%filter(scenario = "baseline_layer")%>%select("total.exposure")))
humanexposure.ac$ratio.det <- humanexposure.ac$detection.exposure/baseline.average.det
humanexposure.ac$ratio.tot <- humanexposure.ac$total.exposure/baseline.average.tot


ggplot(humanexposure.ac)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))
scenario.label <- list("baseline" = "Baseline",
                       "HighTitre90perc" = "90% - Slow",
  #                     "HighTitre90percHighWaning" = "90% - Fast",
                       "HighTitre90percNoWaning" = "90% - No")
ggplot(humanexposure.ac) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~., labeller =  function(variable, value){
    return(scenario.label[value])
  }) +ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer90perc.ac.png")
