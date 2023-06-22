#########################################################
#                                                        
#                  Processing of simulated scenarios                               
#                                                        
#                  Author:        E.A.J. Fischer                       
#                  Contact:       e.a.j.fischer@uu.nl                       
#                  Creation date  12-5-2023                       
#########################################################
#source required 
source("./src/loadLibraries.R") 
source("./src/postprocesSimulations.R")
source("./src/detectionModule.R")
source("./src/probMajorOutbreak.R")

# Layer baseline ####
#load baseline simulations 
rm(output.baseline.layer)
output.baseline.layer <- load.sims("./output/baseline_layer", interval = 0.001)$output


#load baseline parameters (assuming all parameter lists are equal within a folder)
pars.baseline.layer <- load.sims("./output/baseline_layer", params = TRUE)$pars[[1]]

#probability of a major outbreak ####
q1q2.baseline.layer <- q1q2(param.list.baseline.layer)
threshold <- q1q2.baseline.layer%>%select("p","Rv")%>%filter(Rv<=1)%>%first
  
n <- param.list.baseline.layer$initial

q1q2plot <- ggplot(data =q1q2.baseline.layer)+
  geom_path(aes(p,q1^(n), colour = "Low Titre"))+
  geom_path(aes(p,q2^(n), colour = "High Titre"))+
  geom_path(aes(p,(q1^(n*(1-p))* q2^(n*p)), colour = "Low & High titre ratio"))+
  geom_path(aes(p,Rv, colour = "R"))+
  xlab("Proportion with high titre")+
  ylab("Probability of a minor outbreak")+ 
  scale_colour_manual(name = paste("Introduction by",n,"birds \n"),
                      labels = c("High Titre","High and Low Titre","Low Titre","R"),
                      values = c("red","blue","grey","black"))+
  theme()
ggsave("./output/figures/probMinorOutbreakLayer.png" ,q1q2plot)

df = data.frame(scenario  = sapply(c(1:12), FUN = function(i){unlist(scenario.list[[i]][c("scenario")])}),
                p.hightitre = sapply(c(1:12), FUN = function(i){unlist(scenario.list[[i]][c("p.hightitre")])}))

phight0 = unique(df$p.hightitre)
#indicate the status at t = 0 of the scenarios
q1q2plot.anno <- ggplot(data =q1q2.baseline.layer)+
  geom_path(aes(p,q1^(n), colour = "Low Titre"))+
  geom_path(aes(p,q2^(n), colour = "High Titre"))+
  geom_path(aes(p,(q1^(n*(1-p))* q2^(n*p)), colour = "Low & High titre ratio"))+
  geom_path(aes(p,Rv, colour = "R"))+
  geom_vline(xintercept = phight0, linetype = "dotted")+
  xlab("Proportion with high titre")+
  ylab("Probability of a minor outbreak")+ 
  scale_colour_manual(name = paste("Introduction by",n,"birds"),
                      labels = c("High Titre","High and Low Titre","Low Titre","R"),
                      values = c("red","blue","grey","black"))+
  theme()
ggsave("./output/figures/probMinorOutbreakLayerAnno.png" )



#Layer scenarios ####
output.layer <- lapply(c(1:12),function(i){load.sims(paste0("./output/",gsub(scenario.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.1)$output})
pars.layer <- lapply(c(1:12),function(i){load.sims(paste0("./output/",gsub(scenario.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.1)}$pars[[1]])

#visualize ####
for(i in c(1:12)){
  show(plot.output(output.layer[[i]],c("I.1","I.2","R.1","R.2"), scenario.list[[i]]$scenario))
  ggsave(paste0("./output/figures/", gsub(scenario.list[[i]]$scenario,pattern = "[.]", replacement = ""),".png"))
  gc()  
}


#surveillance ###
reps <- 100;
surveillance.layer.baseline <- cbind(repeat.detection.time.surveillance(output.baseline.layer,
                                   reps = reps,
                                   deaths.vars  = c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                   time.interval.pas = 1,
                                   threshold = 0.005*pars.baseline.layer$N0,
                                   ints = 2,
                                   detectables.vars = c("DI.1","DI.2","DR.1","DR.2"),
                                   se = 0.99,
                                   time.interval.ac =7,
                                   init.ac =0,
                                   detectables.incidence = TRUE,
                                   pfarm = 1., 
                                   panimal = 1.),
scenario = pars.baseline.layer$scenario)
rm(surveillance.layer)
for(i in c(1:12)){
  tmp <- cbind(repeat.detection.time.surveillance(output.layer[[i]],
                                                  reps = reps,
                                                  deaths.vars  = c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                  time.interval.pas = 1,
                                                  threshold = 0.005*pars.layer[[i]]$N0,
                                                  ints = 2,
                                                  detectables.vars = c("DI.1","DI.2","DR.1","DR.2"),
                                                  se = 0.99,
                                                  time.interval.ac =7,
                                                  init.ac ="rand",
                                                  detectables.incidence = TRUE,
                                                  pfarm = 1., 
                                                  panimal = 1.),
               scenario = pars.layer[[i]]$scenario,
               WaningRate = max(pars.layer[[i]]$transRate),
               introTime = pars.baseline.layer$max.time-max(pars.layer[[i]]$max.time))
  surveillance.layer <-if(!exists("surveillance.layer")){tmp}else {rbind(surveillance.layer,tmp)}
  
}
surveillance.layer<- cbind(surveillance.layer,surveillance.layer$scenario%>%gsub(pattern ="layer_", replacement = "")%>%str_split_fixed(pattern = c("intro"),2)%>%as.data.frame())



ggplot(data = surveillance.layer.baseline)+
  geom_histogram(aes(x = pas.det.time, y=..count../sum(..count..)),    colour = "black", 
                 fill = "black", 
                 alpha = 0.5, binwidth = 1.0)+
  geom_histogram(aes(ac.det.time, y=..count../sum(..count..)),colour = "black", 
                 fill = "orange",
                 alpha = 0.5, binwidth = 1.0)+
  geom_density(aes(min.det.time,after_stat(..count../sum(..count..))),
               colour = "red",
               linewidth = .5,
               bw = 1)+
  xlim(0,15)+ylim(0,NA)+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  ggtitle("Active and passive surveillance")
ggsave("./output/figures/baseline_layer_surveillance.png")

ggplot(data = surveillance.layer%>%filter(scenario != "baseline_layer") )+
   geom_histogram(aes(x = pas.det.time, y=..count../sum(..count..)),    colour = "black", 
                   fill = "black", 
                   alpha = 0.5, binwidth = 1.0)+
    geom_histogram(aes(ac.det.time, y=..count../sum(..count..)),colour = "black", 
                   fill = "orange",
                   alpha = 0.5, binwidth = 1.0)+
    geom_density(aes(min.det.time,after_stat(..count../sum(..count..))),
                 colour = "red",
                 linewidth = .5,
                bw = 1)+
  xlim(0,NA)+ylim(0,NA)+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  ggtitle("Active and passive surveillance")+
   facet_grid(introTime~WaningRate
              #,labeller =  function(variable,value){return(scenario.label[value])}
     )
ggsave("./output/figures/scenarios_layer_surveillance.png")


#human exposure
rm(humanexposure.min)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.min <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                              pars.baseline.layer$beta,
                                                              surveillance.layer.baseline, var = "min.det.time"), scenario = "baseline")

for(k in c(1:12))
{
  humanexposure.min <- rbind(humanexposure.min,
                             cbind(human.exposure.total.multiple.runs(output.layer[[k]],
                                                                      pars.layer[[k]]$beta,
                                                                      surveillance.layer%>%filter(scenario == pars.layer[[k]]$scenario), var = "min.det.time"), scenario = pars.layer[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.min%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.min%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.min$ratio.det <- humanexposure.min$detection.exposure/baseline.average.det
humanexposure.min$ratio.tot <- humanexposure.min$total.exposure/baseline.average.tot


ggplot(humanexposure.min)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))

ggplot(humanexposure.min) +
  geom_histogram(aes(log10(ratio.det),after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = log10(1))+
  facet_grid(scenario~.)+#, labeller =  function(variable, value){
    #return(scenario.label[value])}) 
 ggtitle("Human exposure")+theme(legend.position = 'none')
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
rm(humanexposure.pas)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.pas <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                              pars.baseline.layer$beta,
                                                              surveillance.layer.baseline, var = "pas.det.time"), scenario = "baseline")

for(k in c(1:12))
{
  humanexposure.pas <- rbind(humanexposure.pas,
                             cbind(human.exposure.total.multiple.runs(output.layer[[k]],
                                                                      pars.layer[[k]]$beta,
                                                                      surveillance.layer%>%filter(scenario == pars.layer[[k]]$scenario), var = "pas.det.time"), scenario = pars.layer[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.pas%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.pas%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.pas$ratio.det <- humanexposure.pas$detection.exposure/baseline.average.det
humanexposure.pas$ratio.tot <- humanexposure.pas$total.exposure/baseline.average.tot


ggplot(humanexposure.pas)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))

ggplot(humanexposure.pas) +
  geom_histogram(aes(log10(ratio.det),after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = log10(1))+
  facet_grid(scenario~.)+#, labeller =  function(variable, value){
  #return(scenario.label[value])}) 
  ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayerPassive.png")

#average ratios
humanexposure.pas%>%reframe(.by = scenario,
                            mean = mean(ratio.det),
                            min = min(ratio.det),
                            max = max(ratio.det),
                            perc25 = quantile(ratio.det,0.25),
                            perc75 =quantile(ratio.det,0.75)
)

#active detection
rm(humanexposure.ac)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.ac <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                              pars.baseline.layer$beta,
                                                              surveillance.layer.baseline, var = "ac.det.time"), scenario = "baseline")

for(k in c(1:12))
{
  humanexposure.ac <- rbind(humanexposure.ac,
                             cbind(human.exposure.total.multiple.runs(output.layer[[k]],
                                                                      pars.layer[[k]]$beta,
                                                                      surveillance.layer%>%filter(scenario == pars.layer[[k]]$scenario), var = "ac.det.time"), scenario = pars.layer[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.ac%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.ac%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.ac$ratio.det <- humanexposure.ac$detection.exposure/baseline.average.det
humanexposure.ac$ratio.tot <- humanexposure.ac$total.exposure/baseline.average.tot


ggplot(humanexposure.ac)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))

ggplot(humanexposure.ac) +
  geom_histogram(aes(log10(ratio.det),after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = log10(1))+
  facet_grid(scenario~.)+#, labeller =  function(variable, value){
  #return(scenario.label[value])}) 
  ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayerActive.png")

#average ratios
humanexposure.ac%>%reframe(.by = scenario,
                            mean = mean(ratio.det),
                            min = min(ratio.det),
                            max = max(ratio.det),
                            perc25 = quantile(ratio.det,0.25),
                            perc75 =quantile(ratio.det,0.75)
)

