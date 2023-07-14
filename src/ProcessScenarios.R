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

#Baseline ####
#load baseline simulations 
rm(output.baseline.layer, output.baseline.broiler)
output.baseline.layer <- load.sims("./output/layerSize32000Vac0")$output
#output.baseline.broiler <- load.sims("./output/baseline_Broiler")$output

#load baseline parameters (assuming all parameter lists are equal within a folder)
pars.baseline.layer <- load.sims("./output/layerSize32000Vac0", params = TRUE)$pars[[1]]
#pars.baseline.broiler <- load.sims("./output/baseline_Broiler", params = TRUE)$pars[[1]]

#probability of a major outbreak ####
q1q2.baseline.layer <- q1q2(param.list.baseline.layer)
threshold <- q1q2.baseline.layer%>%select("p","Rv")%>%filter(Rv<=1)%>%first
  
n <- param.list.baseline.layer$initial

q1q2plot <- ggplot(data =q1q2.baseline.layer)+
  geom_path(aes(p,q1^(n), colour = "Low Titre"),linewidth = 1.0)+
  geom_path(aes(p,q2^(n), colour = "High Titre"),linewidth = 1.0)+
  geom_path(aes(p,(q1^(n*(1-p))* q2^(n*p)), colour = "Low & High titre ratio"),linewidth = 1.0)+
  geom_path(aes(p,Rv, colour = "R"),linewidth = 1.0)+
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
ggsave("./output/figures/probMinorOutbreakLayerAnno.png" 
       )


#visualize baseline ####
plot.output.grid(output = rbind(cbind(output.baseline.layer,scenario = "Layer")#,
                       #cbind(output.baseline.broiler,scenario = "Broiler")
                       ),
                 vars = c("I.1","I.2","DR.1","DR.2"), title = "Baseline", frac = 0.1)
ggsave("./output/figures/baselineLayer.png")




#surveillance ####
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
                                                                        init.ac ="rand",
                                                                        detectables.incidence = TRUE,
                                                                        pfarm = 1., 
                                                                        panimal = 1.),
                                     scenario = pars.baseline.layer$scenario)
# surveillance.broiler.baseline <- cbind(repeat.detection.time.surveillance(output.baseline.broiler,
#                                                                         reps = reps,
#                                                                         deaths.vars  = c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
#                                                                         time.interval.pas = 1,
#                                                                         threshold = 0.005*pars.baseline.broiler$N0,
#                                                                         ints = 2,
#                                                                         detectables.vars = c("DI.1","DI.2","DR.1","DR.2"),
#                                                                         se = 0.99,
#                                                                         time.interval.ac =7,
#                                                                         init.ac ="rand",
#                                                                         detectables.incidence = TRUE,
#                                                                         pfarm = 1., 
#                                                                         panimal = 1.),
#                                      scenario = pars.baseline.broiler$scenario)

ggplot(data = surveillance.layer.baseline)+
  geom_histogram(aes(x = pas.det.time, y=..count../sum(..count..), fill ="Passive"), colour = "black",
                 alpha = 0.5, binwidth = 1.0)+
  geom_histogram(aes(ac.det.time, y=..count../sum(..count..),fill = "Active"), colour = "black", 
                 alpha = 0.5, binwidth = 1.0)+
  geom_density(aes(min.det.time,after_stat(..density..), colour = "Minimum"),
               linewidth = 1.0,
               bw = 1)+
  xlim(0,15)+ylim(0,NA)+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  scale_fill_manual("Detection method",values = c("Passive" = "grey","Active"= "orange"))+
  scale_colour_manual("", values = c("Minimum" = "red"))+
  ggtitle("Active and passive surveillance Layers")
ggsave("./output/figures/baseline_layer_surveillance.png")

# ggplot(data = surveillance.broiler.baseline)+
#   geom_histogram(aes(x = pas.det.time, y=..count../sum(..count..), fill ="Passive", colour = "black" ),
#                  
#                  alpha = 0.5, binwidth = 1.0)+
#   geom_histogram(aes(ac.det.time, y=..count../sum(..count..),fill = "Active", colour = "black"), 
#                  
#                  alpha = 0.5, binwidth = 1.0)+
#    geom_density(aes(min.det.time,after_stat(..density..), colour = "Minimum"),
#                linewidth = 1.0,
#                bw = 1)+
#   xlim(0,15)+ylim(0,NA)+
#   xlab("Detection time")+
#   ylab("Proportion of runs")+
#   scale_fill_manual("Detection method",values = c("Passive" = "grey","Active"= "orange"))+
#   scale_colour_manual("", values = c("Minimum" = "red"))+
#   ggtitle("Active and passive surveillance Broilers")
# 
# ggsave("./output/figures/baseline_broiler_surveillance.png")



#Layer scenarios with %-high titre####
output.layer <- lapply(c(1:length(scenario.list.size.vaccination)),function(i){load.sims(paste0("./output/",gsub(scenario.list.size.vaccination[[i]]$scenario,pattern = "[.]", replacement = "")))$output})
pars.layer <- lapply(c(1:length(scenario.list.size.vaccination)),function(i){load.sims(paste0("./output/",gsub(scenario.list.size.vaccination[[i]]$scenario,pattern = "[.]", replacement = "")))}$pars[[1]])

#visualize ####
for(i in c(1:length(scenario.list.size.vaccination))){
  show(plot.output.sparse(output.layer[[i]],c("I.1","I.2","R.1","R.2"), scenario.list.size.vaccination[[i]]$scenario, frac= 0.1))
  ggsave(paste0("./output/figures/", gsub(scenario.list.size.vaccination[[i]]$scenario,pattern = "[.]", replacement = ""),".png"))
  gc()  
}


#surveillance ####
reps <- 100;
rm(surveillance.layer)
for(i in c(1:length(scenario.list.size.vaccination))){
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
surveillance.layer$vaccination <- unlist(str_split_fixed(surveillance.layer$scenario,pattern= c("Vac"),2))[,2]
surveillance.layer$size <-gsub(unlist(str_split_fixed(surveillance.layer$scenario,pattern= c("Vac"),2))[,1], pattern = "layerSize", replacement = "")
 



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
   facet_grid(vaccination~size
              #,labeller =  function(variable,value){return(scenario.label[value])}
     )
ggsave("./output/figures/scenarios_layer_surveillance.png")


dists <- data.frame(shape =c(), rate =c(), scenario =c())
for(its in unique(surveillance.layer$scenario))  {
  x = surveillance.layer%>%filter(scenario == its)%>%dplyr::select("min.det.time")%>%filter(min.det.time!=Inf)
  if(length(x$min.det.time)>1)
  {
  fit.gamma <- fitdistrplus::fitdist(data = x$min.det.time, distr = "gamma", method = "mle")
  sgamma <-summary(fit.gamma)
  dists<- rbind(dists,c(shape = as.numeric(sgamma$estimate[1]),rate = as.numeric(sgamma$estimate[2]), scenario = its))
  plot(fit.gamma)
  }else dists <- rbind(dists,c(shape =c(NA), rate =c(NA), scenario = its))
}
names(dists)<-c("shape","rate","scenario")
dists$shape <- as.numeric(dists$shape)
dists$rate <- as.numeric(dists$rate)
dists$mean <- dists$shape/dists$rate

#human exposure
rm(humanexposure.min.layer)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.min.layer <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                              pars.baseline.layer$beta,
                                                              surveillance.layer.baseline, var = "min.det.time"), scenario = "baseline")

for(k in c(1:12))
{
  humanexposure.min.layer <- rbind(humanexposure.min.layer,
                             cbind(human.exposure.total.multiple.runs(output.layer[[k]],
                                                                      pars.layer[[k]]$beta,
                                                                      surveillance.layer%>%filter(scenario == pars.layer[[k]]$scenario), var = "min.det.time"), scenario = pars.layer[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.min.layer%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.min.layer%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.min.layer$ratio.det <- humanexposure.min.layer$detection.exposure/baseline.average.det
humanexposure.min.layer$ratio.tot <- humanexposure.min.layer$total.exposure/baseline.average.tot


ggplot(humanexposure.min.layer)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))

ggplot(humanexposure.min.layer) +
  geom_histogram(aes(log10(ratio.det),after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = log10(1))+
  facet_grid(scenario~.)+#, labeller =  function(variable, value){
    #return(scenario.label[value])}) 
 ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer.png")

#average ratios
humanexposure.min.layer%>%reframe(.by = scenario,
  mean = mean(ratio.det),
  min = min(ratio.det),
  max = max(ratio.det),
  perc25 = quantile(ratio.det,0.25),
  perc75 =quantile(ratio.det,0.75)
  )


#passive detection
rm(humanexposure.pas.layer)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.pas.layer <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                              pars.baseline.layer$beta,
                                                              surveillance.layer.baseline, var = "pas.det.time"), scenario = "baseline")

for(k in c(1:12))
{
  humanexposure.pas.layer <- rbind(humanexposure.pas.layer,
                             cbind(human.exposure.total.multiple.runs(output.layer[[k]],
                                                                      pars.layer[[k]]$beta,
                                                                      surveillance.layer%>%filter(scenario == pars.layer[[k]]$scenario), var = "pas.det.time"), scenario = pars.layer[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.pas.layer%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.pas.layer%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.pas.layer$ratio.det <- humanexposure.pas.layer$detection.exposure/baseline.average.det
humanexposure.pas.layer$ratio.tot <- humanexposure.pas.layer$total.exposure/baseline.average.tot


ggplot(humanexposure.pas.layer)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))

ggplot(humanexposure.pas.layer) +
  geom_histogram(aes(log10(ratio.det),after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = log10(1))+
  facet_grid(scenario~.)+#, labeller =  function(variable, value){
  #return(scenario.label[value])}) 
  ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayerPassive.png")

#average ratios
humanexposure.pas.layer%>%reframe(.by = scenario,
                            mean = mean(ratio.det),
                            min = min(ratio.det),
                            max = max(ratio.det),
                            perc25 = quantile(ratio.det,0.25),
                            perc75 =quantile(ratio.det,0.75)
)

#active detection
rm(humanexposure.ac.layer)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.ac.layer <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                              pars.baseline.layer$beta,
                                                              surveillance.layer.baseline, var = "ac.det.time"), scenario = "baseline")

for(k in c(1:12))
{
  humanexposure.ac.layer <- rbind(humanexposure.ac.layer,
                             cbind(human.exposure.total.multiple.runs(output.layer[[k]],
                                                                      pars.layer[[k]]$beta,
                                                                      surveillance.layer%>%filter(scenario == pars.layer[[k]]$scenario), var = "ac.det.time"), scenario = pars.layer[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.ac.layer%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.ac.layer%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.ac.layer$ratio.det <- humanexposure.ac.layer$detection.exposure/baseline.average.det
humanexposure.ac.layer$ratio.tot <- humanexposure.ac.layer$total.exposure/baseline.average.tot


ggplot(humanexposure.ac.layer)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))

ggplot(humanexposure.ac.layer) +
  geom_histogram(aes(log10(ratio.det),after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = log10(1))+
  facet_grid(scenario~.)+#, labeller =  function(variable, value){
  #return(scenario.label[value])}) 
  ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayerActive.png")

#average ratios
humanexposure.ac.layer%>%reframe(.by = scenario,
                            mean = mean(ratio.det),
                            min = min(ratio.det),
                            max = max(ratio.det),
                            perc25 = quantile(ratio.det,0.25),
                            perc75 =quantile(ratio.det,0.75)
)



#Scenarios with clinical protection ####
scenario.clinprot.list <- list(list(scenario = "layerClinic50Phigh50"),list(scenario = "layerClinic50Phigh0"),list(scenario = "layerClinic05Phigh0"))
output.layer.clinprot <- lapply(c(1:3),function(i){load.sims(paste0("./output/",gsub(scenario.clinprot.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)$output})
pars.layer.clinprot <- lapply(c(1:3),function(i){load.sims(paste0("./output/",gsub(scenario.clinprot.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)}$pars[[1]])

#visualize ####
for(i in c(1:3)){
  show(plot.output(output.layer.clinprot[[i]],c("I.1","I.2","R.1","R.2"), scenario.clinprot.list[[i]]$scenario))
  ggsave(paste0("./output/figures/", gsub(scenario.clinprot.list[[i]]$scenario.clinprot,pattern = "[.]", replacement = ""),".png"))
  gc()  
}

plot.output.grid(rbind(cbind(output.layer.clinprot[[1]],scenario = pars.layer.clinprot[[1]]$scenario),
                       cbind(output.layer.clinprot[[2]],scenario = pars.layer.clinprot[[2]]$scenario),
                       cbind(output.layer.clinprot[[3]],scenario = pars.layer.clinprot[[3]]$scenario)), 
                 vars =c("I.1","I.2","R.1","R.2"),
                 title = "Clinical protection")
ggsave(file = "./output/figures/clinicalprotectionLayer.png")

#surveillance ####
reps <- 100;

rm(surveillance.layer.clinprot)
for(i in c(1:3)){
  tmp <- cbind(repeat.detection.time.surveillance(output.layer.clinprot[[i]],
                                                  reps = reps,
                                                  deaths.vars  = c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                  time.interval.pas = 1,
                                                  threshold = 0.005*pars.layer.clinprot[[i]]$N0,
                                                  ints = 2,
                                                  detectables.vars = c("DI.1","DI.2","DR.1","DR.2"),
                                                  se = 0.99,
                                                  time.interval.ac =7,
                                                  init.ac ="rand",
                                                  detectables.incidence = TRUE,
                                                  pfarm = 1., 
                                                  panimal = 1.),
               scenario = pars.layer.clinprot[[i]]$scenario,
               WaningRate = max(pars.layer.clinprot[[i]]$transRate),
               introTime = pars.baseline.layer$max.time-max(pars.layer.clinprot[[i]]$max.time))
  surveillance.layer.clinprot <-if(!exists("surveillance.layer.clinprot")){tmp}else {rbind(surveillance.layer.clinprot,tmp)}
  
}
surveillance.layer.clinprot<- cbind(surveillance.layer.clinprot,surveillance.layer.clinprot$scenario%>%gsub(pattern ="layer.clinprot_", replacement = "")%>%str_split_fixed(pattern = c("intro"),2)%>%as.data.frame())





ggplot(data = surveillance.layer.clinprot%>%filter(scenario != "baseline_layer.clinprot") )+
  geom_histogram(aes(x = pas.det.time, y=..count../sum(..count..)),    colour = "black", 
                 fill = "black", 
                 alpha = 0.5, binwidth = 1.0)+
  geom_histogram(aes(ac.det.time, y=..count../sum(..count..)),colour = "black", 
                 fill = "orange",
                 alpha = 0.5, binwidth = 1.0)+
  geom_density(aes(min.det.time,after_stat(..density..)),
               colour = "red",
               linewidth = 1.0,
               bw = 1)+
  xlim(0,NA)+ylim(0,NA)+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  ggtitle("Active and passive surveillance")+
  facet_grid(scenario~.
             #,labeller =  function(variable,value){return(scenario.label[value])}
  )
ggsave("./output/figures/scenarios_layer.clinprot_surveillance.png")


#human exposure ####
rm(humanexposure.min.layer)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.min.layer <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                              pars.baseline.layer$beta,
                                                              surveillance.layer.baseline, var = "min.det.time"), scenario = "baseline")

for(k in c(1:3)){
  humanexposure.min.layer <- rbind(humanexposure.min.layer,
                             cbind(human.exposure.total.multiple.runs(output.layer.clinprot[[k]],
                                                                      pars.layer.clinprot[[k]]$beta,
                                                                      surveillance.layer.clinprot%>%filter(scenario == pars.layer.clinprot[[k]]$scenario), var = "min.det.time"), scenario = pars.layer.clinprot[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.min.layer%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.min.layer%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.min.layer$ratio.det <- humanexposure.min.layer$detection.exposure/baseline.average.det
humanexposure.min.layer$ratio.tot <- humanexposure.min.layer$total.exposure/baseline.average.tot


ggplot(humanexposure.min.layer)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))

scenario.label <- list(baseline = "Baseline",
                       layerClinic50Phigh50 = "Clinic 50% \n High titre 50%",
                       layerClinic50Phigh0 = "Clinic 50% \n High titre 0%",
                       layerClinic0.5Phigh0 = "Clinic 0.5% \n High titre 0%")
ggplot(humanexposure.min.layer) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = scenario),binwidth = .5)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~., labeller =  function(variable, value){
    return(scenario.label[value])}) +
  ggtitle("Human exposure layers \n passive and active detection")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer.clinprot.min.png")

#average ratios
humanexposure.min.layer%>%reframe(.by = scenario,
                            mean = mean(ratio.det),
                            min = min(ratio.det),
                            max = max(ratio.det),
                            perc25 = quantile(ratio.det,0.25),
                            perc75 =quantile(ratio.det,0.75)
)


#passive detection####
rm(humanexposure.pas.layer)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.pas.layer <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                              pars.baseline.layer$beta,
                                                              surveillance.layer.baseline, var = "pas.det.time"), 
                                 scenario = "baseline")

for(k in c(1:3)){
  humanexposure.pas.layer <- rbind(humanexposure.pas.layer,
                                   cbind(human.exposure.total.multiple.runs(output.layer.clinprot[[k]],
                                                                            pars.layer.clinprot[[k]]$beta,
                                                                            surveillance.layer.clinprot%>%filter(scenario == pars.layer.clinprot[[k]]$scenario), var = "pas.det.time"), scenario = pars.layer.clinprot[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.pas.layer%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.pas.layer%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.pas.layer$ratio.det <- humanexposure.pas.layer$detection.exposure/baseline.average.det
humanexposure.pas.layer$ratio.tot <- humanexposure.pas.layer$total.exposure/baseline.average.tot


ggplot(humanexposure.pas.layer)+geom_point(aes(x = detection.time, y = ratio.det, colour = scenario))
scenario.label <- list(baseline = "Baseline",
                       layerClinic50Phigh50 = "Clinic 50% \n High titre 50%",
                       layerClinic50Phigh0 = "Clinic 50% \n High titre 0%",
                       layerClinic0.5Phigh0 = "Clinic 0.5% \n High titre 0%")
ggplot(humanexposure.pas.layer) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = scenario),binwidth = .5)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~., labeller =  function(variable, value){
  return(scenario.label[value])}) +
  ggtitle("Human exposure layers passive detection")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer.clinprot.pas.png")

#average ratios
humanexposure.pas.layer%>%reframe(.by = scenario,
                                  mean = mean(ratio.det),
                                  min = min(ratio.det),
                                  max = max(ratio.det),
                                  perc25 = quantile(ratio.det,0.25),
                                  perc75 =quantile(ratio.det,0.75)
)

#active detection
rm(humanexposure.ac.layer)
#as we used multiple repeats to determine active surveillance repeat
humanexposure.ac.layer <- cbind(human.exposure.total.multiple.runs(output.baseline.layer,
                                                             pars.baseline.layer$beta,
                                                             surveillance.layer.baseline, var = "ac.det.time"), scenario = "baseline")

for(k in c(1:3))
{
  humanexposure.ac.layer <- rbind(humanexposure.ac.layer,
                            cbind(human.exposure.total.multiple.runs(output.layer.clinprot[[k]],
                                                                     pars.layer.clinprot[[k]]$beta,
                                                                     surveillance.layer.clinprot%>%filter(scenario == pars.layer.clinprot[[k]]$scenario), var = "ac.det.time"), scenario = pars.layer.clinprot[[k]]$scenario))
  
}



baseline.average.det <- mean(unlist(humanexposure.ac.layer%>%filter(scenario == "baseline")%>%select("detection.exposure")))
baseline.average.tot <- mean(unlist(humanexposure.ac.layer%>%filter(scenario == "baseline")%>%select("total.exposure")))
humanexposure.ac.layer$ratio.det <- humanexposure.ac.layer$detection.exposure/baseline.average.det
humanexposure.ac.layer$ratio.tot <- humanexposure.ac.layer$total.exposure/baseline.average.tot


ggplot(humanexposure.ac.layer)+
  geom_point(aes(x = detection.time, y = cum.I.1+cum.I.2, colour = scenario))+facet_grid(scenario~.)

ggplot(humanexposure.ac.layer) +
  geom_histogram(aes(log10(ratio.det),after_stat(.15*density), fill = scenario),binwidth = .15)+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = log10(1))+
  facet_grid(scenario~.)+#, labeller =  function(variable, value){
  #return(scenario.label[value])}) 
  ggtitle("Human exposure")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer.clinprotActive.png")

#average ratios
humanexposure.ac.layer%>%reframe(.by = scenario,
                           mean = mean(ratio.det),
                           min = min(ratio.det),
                           max = max(ratio.det),
                           perc25 = quantile(ratio.det,0.25),
                           perc75 =quantile(ratio.det,0.75)
)




