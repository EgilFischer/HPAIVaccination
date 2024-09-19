#source required 
source("./src/loadLibraries.R") 
source("./src/postprocesSimulations.R")
source("./src/detectionModule.R")
source("./src/probMajorOutbreak.R")
source("./src/LabelsForPlots.R")

#Baseline ####
baseline.layer<- readRDS("./output/coverage/20231210outputlayerSize32000Vac0.RDS")
#load baseline parameters (assuming all parameter lists are equal within a folder)
pars.baseline.layer <- baseline.layer$pars
#calculate mean and variance of infectious period
pars.baseline.layer$variance.infectious.period <- with(pars.baseline.layer, k.infectious/(k.infectious/infectious.period)^2)
pars.baseline.layer$variance.infectious.period <- with(pars.baseline.layer, k.infectious/(k.infectious/infectious.period)^2)


#trans vs clinic ###
#Layer scenarios with %-high titre####
layer.clintrans <- load.sims("./output/clinicJan24/")
output.layer.clintrans <- layer.clintrans$output
pars.layer.clintrans <- layer.clintrans$pars

#visualize ####
#plot the baseline size####
plot.data <- NULL;
for(j in c(1:length(output.layer.clintrans))){
  tmp <- cbind(data.frame(output.layer.clintrans[[j]]),
               scenario = pars.layer.clintrans[[j]]$scenario,
               clin1 = pars.layer.clintrans[[j]]$pdie[[1]],
               clin2 = pars.layer.clintrans[[j]]$pdie[[2]],
               phigh = pars.layer.clintrans[[j]]$p.hightitre)
    plot.data<- rbind(plot.data,tmp)
}
lf.plot.data<- plot.data%>%reshape2::melt(id.vars = c("time","run","dt","scenario", "clin1","clin2","phigh") )

#create a data set with the total of infected animals (e.g. sum(N) - sum(S)) and maximum dead animals
plot.data.end <- lf.plot.data%>%
  reframe(.by = c(run, time, scenario,clin1,clin2, phigh),
                                        n = sum(value[variable %in% c("N.1","N.2")]),
                                        s = sum(value[variable %in% c("S.1","S.2")]),
                                        dead = sum(value[variable %in% c("DS.1","DL.1","DI.1","DR.1","DS.2","DL.2","DI.2","DR.2")]))%>%
  reframe(.by = c(run,  scenario,clin1,clin2, phigh),
          infected = max(s)+10-min(s),
          dead = max(dead)
          )
saveRDS(plot.data.end, "plot.data.end.RDS")
plot.data.end<-readRDS("plot.data.end.RDS")
#summarize plot end data
plot.data.end.sum <- plot.data.end%>%reframe(.by=c(scenario,clin1,clin2, phigh),
                                             meanInfected = mean(infected),
                                             medianInfected = median(infected),
                                             meanDead = mean(dead),
                                             medianDead = median(dead))

head(plot.data.end.sum)
saveRDS(plot.data.end.sum, "plot.data.end.sum.RDS")

#visualize
ggplot(data = plot.data.end.sum, aes(x = clin1,  y = medianDead, colour = as.factor(clin2)))+
  geom_point()+
  facet_grid(phigh~.)


#improve this test data.
reframe.lf.plot.data <- lf.plot.data%>%reframe(.by = c(run, time, scenario,clin1,clin2, phigh, value),
                                        Deaths1 = variable %in% c("DI.1","DR.1")) %>%reframe(.by = c(run, time, scenario,clin1,clin2, phigh, Deaths1), sum(value))

ggplot(reframe.lf.plot.data%>%filter(phigh == 0.5))+
  geom_path(aes(x = time, y = `sum(value)`, colour = c("High titre","Low titre")[as.integer(Deaths1)+1],  group = paste0(run,Deaths1)))+
  facet_nested(clin1 ~ clin2) +theme(legend.position = "none")


ggplot(data = lf.plot.data%>%filter(variable %in% c("I.1","I.2","R.1","R.2"))%>%arrange(time,run) )+
  geom_path(aes(x = time, y = value, colour = variable,  group = paste0(run,variable)))+
 # scale_colour_manual(values = c("darkblue","darkred"),labels = c("Broiler","Layer"),name = "Type")+
  #facet_grid(type+as.factor(size) ~ variable, labeller = labeller(variable = SIR_labels)) 
  facet_nested(clin~ phigh)#, labeller = labeller(variable = SIR_labels)) 
ggsave("./output/figures/baselineoutbreaklayerclintransIR.png")

ggplot(data = lf.plot.data%>%filter(variable %in% c("DI.1","DI.2","DR.1","DR.2")&run <=1)%>%arrange(time,run)%>%reframe(.by = run,
                                                                                                                         time = time,
                                                                                                                         Death.1 = DI.1+DR.1,
                                                                                                                         Death.2 = DI.2 +DR.2) )+
  geom_path(aes(x = time, y = value, colour = variable,  group = paste0(run,variable)))+
  # scale_colour_manual(values = c("darkblue","darkred"),labels = c("Broiler","Layer"),name = "Type")+
  #facet_grid(type+as.factor(size) ~ variable, labeller = labeller(variable = SIR_labels)) 
  facet_nested(clin~ phigh)#, labeller = labeller(variable = SIR_labels)) 
ggsave("./output/figures/baselineoutbreaklayerclintransDeaths.png")


#detection of baseline scenarios####
reps <- 1;
surveillance.layer.clintrans<-NULL
for(i in c(1:length(output.layer.clintrans))){
  tmp <- cbind(repeat.detection.time.surveillance(as.data.frame(output.layer.clintrans[[i]]),
                                                  reps = reps,
                                                  deaths.vars  = c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                  time.interval.pas = 1,
                                                  threshold = 0.005*pars.layer.clintrans[[i]]$N0,
                                                  ints = 2,
                                                  detectables.vars = c("DI.1","DI.2","DR.1","DR.2"),
                                                  se = 0.99,
                                                  time.interval.ac =7,
                                                  init.ac ="rand",
                                                  detectables.incidence = TRUE,
                                                  pfarm = 1., 
                                                  panimal = 1.),
               scenario = pars.layer.clintrans[[i]]$scenario,
               max.time = pars.layer.clintrans[[i]]$max.time,
               size = pars.layer.clintrans[[i]]$N0,
               vaccination = pars.layer.clintrans[[i]]$p.hightitre,
               introTime =0)
  surveillance.layer.clintrans <-if(!exists("surveillance.layer.clintrans")|is.null(surveillance.layer.clintrans)){tmp}else {rbind(surveillance.layer.clintrans,tmp)}
  
}