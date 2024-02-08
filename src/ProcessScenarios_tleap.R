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

#labels for plotting ####
itype.labels <- c("I.1" ="I.1",
                  "I.2" = "I.2",
                  "R.1"="R.1",
                  "R.2"="R.2",
                  "DR.1"= "Dead 1",
                  "DR.2" = "Dead 2")
SIR_labels <- c("S.1" ="S",
                "I.1"="I",
                "R.1"="R",
                "DR.1"= "Dead",
                "DI.1" = "Dead")
labels.type <-c("46"  = "Broiler", "510" = "Layer")
labels.size <-c("15000"  = "Small", "32000" = "Medium", "64000" = "Large", "20000"  = "Small", "38000" = "Medium", "73000" = "Large")
scenario.label.clin <- c(layerSize32000Vac0 = "Baseline",
                         layerClinic50Phigh50 = "Clinic 50% \n High titre 50%",
                         layerClinic50Phigh0 = "Clinic 50% \n High titre 0%",
                         layerClinic0.1Phigh0 = "Clinic 0.1% \n High titre 0%")


scenario.label.wane <- c(layerSize32000Vac0 = "Baseline",
                    layerStartTime0 = "T0 = 0",
                    layerStartTime50 = "T0 =  50",
                    layerStartTime100 = "T0 = 100",
                    layerStartTime200 = "T0 = 200",
                    layerStartTime300 = "T0 = 300",
                    layerStartTime400 = "T0 = 400",
                    layerStartTime500 = "T0 = 500")

scenario.levels.label.wane <- names(scenario.label.wane)

scenario.label.waneDifStrain <- c(layerSize32000Vac0 = "Baseline",
                         layerStartTime0DifStrain = "T0 = 0",
                         layerStartTime50DifStrain = "T0 =  50",
                         layerStartTime100DifStrain = "T0 = 100",
                         layerStartTime200DifStrain = "T0 = 200",
                         layerStartTime300DifStrain = "T0 = 300",
                         layerStartTime400DifStrain = "T0 = 400",
                         layerStartTime500DifStrain = "T0 = 500")

scenario.levels.label.waneDifStrain <- names(scenario.label.waneDifStrain)



detection.method.label <- c(pas.det.time = "Passive",
                            ac.det.time = "Active",
                            min.det.time = "Minimum")



#Baseline ####
baseline.layer<- readRDS("./output/coverage/20231210outputlayerSize32000Vac0.RDS")
baseline.broiler<- readRDS("./output/coverage/20231210outputbroilerSize38000Vac0.RDS")
#load baseline simulations 
output.baseline.layer <- baseline.layer$out
output.baseline.broiler <- baseline.broiler$out

#load baseline parameters (assuming all parameter lists are equal within a folder)
pars.baseline.layer <- baseline.layer$pars
pars.baseline.broiler <- baseline.broiler$pars
#calculate mean and variance of infectious period
pars.baseline.layer$variance.infectious.period <- with(pars.baseline.layer, k.infectious/(k.infectious/infectious.period)^2)
pars.baseline.layer$variance.infectious.period <- with(pars.baseline.broiler, k.infectious/(k.infectious/infectious.period)^2)



#probability of a major outbreak ####
q1q2.baseline.layer <- q1q2(pars.baseline.layer)
threshold <- q1q2.baseline.layer%>%select("p","Rv")%>%filter(Rv<=1)%>%first
  
n <- param.list.baseline.layer$initial

q1q2plot <- ggplot(data =q1q2.baseline.layer)+
  geom_path(aes(p,q1^(n), colour = "Low Titre"),linewidth = 1.0)+
  geom_path(aes(p,q2^(n), colour = "High Titre"),linewidth = 1.0)+
  geom_path(aes(p, (q1^(n*(1-p))* q2^(n*p)), colour = "Low & High titre ratio"),linewidth = 1.0)+
  geom_path(aes(p,Rv, colour = "R"),linewidth = 1.0, linetype = "dashed")+
  xlab("Proportion with high titre")+
  ylab("Probability of a minor outbreak")+ 
  scale_colour_manual(name = paste("Introduction by",n,"birds"),
                      labels = c("High Titre","High and Low Titre","Low Titre","R"),
                      values = c("red","darkgrey","blue","black"))+
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "R"))+
  theme()
ggsave("./output/figures/probMinorOutbreakLayer.png" ,q1q2plot, scale = 1.23)



##########################################
#Layer and Broiler different sizes       #
##########################################
scenarios.tmp <-  data.frame(
  type = c(rep("layer",3),rep("broiler",3)),
  size = c(15000,32000,64000,
           20000,38000,73000))

scenario.list.size.type <- mapply(FUN = function(type, size){list(data.frame(scenario  = paste0(type,"Size",size,"Vac0")))}, scenarios.tmp$type,scenarios.tmp$size)

output.size.type <- lapply(c(1:length(scenario.list.size.type)),function(i){readRDS(paste0("./output/coverage/20231210output",gsub(scenario.list.size.type[[i]]$scenario,pattern = "[.]", replacement = ""),".RDS"))$out})
pars.size.type <- lapply(c(1:length(scenario.list.size.type)),function(i){readRDS(paste0("./output/coverage/20231210output",gsub(scenario.list.size.type[[i]]$scenario,pattern = "[.]", replacement = ""),".RDS"))}$pars)

#plot the baseline size####
plot.data <- NULL;
for(j in c(1:length(scenario.list.size.type))){
  tmp <- cbind(data.frame(output.size.type[[j]]),
               scenario = scenario.list.size.type[[j]]$scenario,
               size = scenarios.tmp$size[[j]],
               type = firstup( scenarios.tmp$type[[j]]))
  #tmp <- sample_n(tmp,ceiling(0.01*length(tmp$time)))
  plot.data<- rbind(plot.data,tmp)
}
lf.plot.data<- plot.data%>%reshape2::melt(id.vars = c("time","run","dt","scenario", "size","type") )
inf.lf.plot.data<- lf.plot.data%>%filter(variable %in% c("I.1","I.2","R.1","R.2"))

ggplot(data = lf.plot.data%>%filter(variable %in% c("S.1","I.1","R.1","DI.1"))%>%arrange(time,run) )+
  geom_step(aes(x = time, y = value, colour = type,group = as.factor(run)))+
  scale_colour_manual(values = c("darkblue","darkred"),labels = c("Broiler","Layer"),name = "Type")+
  #facet_grid(type+as.factor(size) ~ variable, labeller = labeller(variable = SIR_labels)) 
  facet_nested(type*as.factor(size) ~ variable, labeller = labeller(variable = SIR_labels)) 
ggsave("./output/figures/baselineoutbreaklayerbroiler.png")

#detection of baseline scenarios####
reps <- 100;
surveillance.size.type<-NULL
for(i in c(1:length(scenario.list.size.type))){
  tmp <- cbind(repeat.detection.time.surveillance(as.data.frame(output.size.type[[i]]),
                                                  reps = reps,
                                                  deaths.vars  = c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                  time.interval.pas = 1,
                                                  threshold = 0.005*pars.size.type[[i]]$N0,
                                                  ints = 2,
                                                  detectables.vars = c("DI.1","DI.2","DR.1","DR.2"),
                                                  se = 0.99,
                                                  time.interval.ac =7,
                                                  init.ac ="rand",
                                                  detectables.incidence = TRUE,
                                                  pfarm = 1., 
                                                  panimal = 1.),
               scenario = pars.size.type[[i]]$scenario,
               max.time = pars.size.type[[i]]$max.time,
               size = pars.size.type[[i]]$N0,
               vaccination = pars.size.type[[i]]$p.hightitre,
               introTime =0)
  surveillance.size.type <-if(!exists("surveillance.size.type")|is.null(surveillance.size.type)){tmp}else {rbind(surveillance.size.type,tmp)}
  
}


#save surveillance times 
saveRDS(surveillance.size.type, file = "./output/surveillanceSizeType.RDS")
#read surveillance times
surveillance.size.type <- readRDS(file = "./output/surveillanceSizeType.RDS")


#write code to show that size does not matter for min detetion times
ggplot(data = surveillance.size.type%>%
         select(c("run","rep","size","max.time","scenario","pas.det.time","ac.det.time","min.det.time"))%>%
         reshape2::melt(id.vars = c("run","rep","size","max.time","scenario")) )+
  geom_histogram(aes(x = as.numeric(value), y=after_stat(count)/sum(after_stat(count)), fill = variable),    colour = "black", 
                 alpha = 0.5, binwidth = 1.0,position = "identity")+
  xlim(0,NA)+ylim(0,NA)+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"), name = "Detection")+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  #ggtitle("Surveillance")+
  facet_grid(size + max.time~variable, labeller = labeller(max.time = labels.type, variable = detection.method.label))
ggsave("./output/figures/scenarios_sizetype_surveillance.png", scale = 1.23)

surveillance.size.type%>%reframe(.by = scenario,
                                 mean = mean(pas.det.time), 
                                 var = var(pas.det.time),
                                 q1 = quantile(pas.det.time,c(0.025)), 
                                 q2 = quantile(pas.det.time,c(0.975)))

#fit distributions ####
add.vars <- list(size = list(layerSize15000Vac0 = 15000,layerSize32000Vac0 = 32000,layerSize64000Vac0 = 64000,
                             broilerSize20000Vac0 = 20000,broilerSize38000Vac0 = 38000,broilerSize73000Vac0 = 73000),
                 type = list(layerSize15000Vac0 = "LAYER",layerSize32000Vac0 = "LAYER",layerSize64000Vac0 = "LAYER",
                             broilerSize20000Vac0 ="BROILER",broilerSize38000Vac0 = "BROILER",broilerSize73000Vac0 = "BROILER"))
dists<- fit.distribution(surveillance.size.type, added.vars = add.vars)

saveRDS(dists, "./output/detDist_sizetype.RDS")
dists <- readRDS("./output/detDist_sizetype.RDS")
ggplot(dists)+geom_point(aes(size, mean, colour = detection.method))+facet_grid(type~.)
fit.pas.mean.sizetype <- lm(mean ~ size + type , data = dists%>%filter(detection.method =="pas.det.time"))
drop1(fit.pas.mean.sizetype)
summary(fit.pas.mean.sizetype)
saveRDS(fit.pas.mean.sizetype, "./output/detMean_sizetype_passive.RDS")
fit.pas.mean.sizetype<- readRDS("./output/detMean_sizetype_passive.RDS")

fit.min.mean.sizetype <- lm(mean ~ size + type , data = dists%>%filter(detection.method =="min.det.time"))
drop1(fit.min.mean.sizetype)
summary(fit.min.mean.sizetype)
saveRDS(fit.min.mean.sizetype, "./output/detMean_sizetype_minimum.RDS")
fit.min.mean.sizetype <- readRDS("./output/detMean_sizetype_minimum.RDS")

fit.pas.var.sizetype <- lm(var ~ size + type  , data = dists%>%filter(detection.method =="pas.det.time"))
drop1(fit.pas.var.sizetype)
summary(fit.pas.var.sizetype)
saveRDS(fit.pas.var.sizetype, "./output/detVar_sizetype_passive.RDS")
fit.pas.var.sizetype <- readRDS("./output/detVar_sizetype_passive.RDS")

fit.min.var.sizetype <- lm(var ~ size + type , data = dists%>%filter(detection.method =="min.det.time"))
drop1(fit.min.var.sizetype)
summary(fit.min.var.sizetype)
fit.min.var.sizetype <- lm(var ~ size , data = dists%>%filter(detection.method =="min.det.time"))
drop1(fit.min.var.sizetype)
summary(fit.min.var.sizetype)

saveRDS(fit.min.var.sizetype, "./output/detVar_sizetype_minimum.RDS")
fit.min.var.sizetype <- readRDS("./output/detVar_sizetype_minimum.RDS")

#human exposure ####
humanexposure.size.type  <- NULL

for(k in c(1:length(output.size.type))){
  humanexposure.size.type <- rbind(humanexposure.size.type,
                               cbind(human.exposure.detection.multiple.runs(as.data.frame(output.size.type[[k]]),
                                                                            pars.size.type[[k]]$beta,
                                                                            surveillance.size.type%>%filter(scenario == pars.size.type[[k]]$scenario)),
                                     scenario = pars.size.type[[k]]$scenario,
                                     size = pars.size.type[[k]]$N0,
                                     vaccination = pars.size.type[[k]]$p.hightitre))
  
}
saveRDS(humanexposure.size.type, "./output/exposure_sizetype.RDS")
humanexposure.size.type<- readRDS( "./output/exposure_sizetype.RDS")

ggplot(humanexposure.size.type)+
  geom_histogram(aes(x = detection.exposure, fill = scenario))+
  facet_nested(size ~.)+
  ggtitle("Unvaccinated")

dists.exposure <- fit.distribution(humanexposure.size.type, long.form = TRUE, scale = 1000,  added.vars = add.vars)
ggplot(dists.exposure)+geom_point(aes(size, mean, colour = detection.method))+facet_grid(type~.)
saveRDS(dists.exposure, "./output/exposureDist_sizetype.RDS")
dists.exposure <-readRDS("./output/exposureDist_sizetype.RDS")

fit.pas.mean.sizetype <- lm(mean ~ size + type , data = dists.exposure%>%filter(detection.method =="pas.det.time"))
drop1(fit.pas.mean.sizetype)
summary(fit.pas.mean.sizetype)
ggplot()+
  geom_point(aes(size, mean, colour =as.factor(type)),dists.exposure%>%filter(detection.method =="pas.det.time"))+
  geom_point(aes(size, val), data.frame(size =dists.exposure%>%filter(detection.method =="pas.det.time")%>%select("size"), val =predict(fit.pas.mean.sizetype) ))


saveRDS(fit.pas.mean.sizetype, "./output/exposureMean_sizetype_passive.RDS")

fit.pas.mean.sizetype<-readRDS("./output/exposureMean_sizetype_passive.RDS")

fit.min.mean.sizetype <- lm(mean ~ size + type , data = dists.exposure%>%filter(detection.method =="min.det.time"))
drop1(fit.min.mean.sizetype)
summary(fit.min.mean.sizetype)


saveRDS(fit.min.mean.sizetype, "./output/exposureMean_sizetype_minimum.RDS")
fit.min.mean.sizetype<-readRDS("./output/exposureMean_sizetype_minimum.RDS")

fit.pas.var.sizetype <- lm(sqrt(var) ~ size + type , data = dists.exposure%>%filter(detection.method =="pas.det.time"))
drop1(fit.pas.var.sizetype)
summary(fit.pas.var.sizetype)


saveRDS(fit.pas.var.sizetype, "./output/exposureSD_sizetype_passive.RDS")




##########################################
#Layer different sizes and % high-titre  #
##########################################
#Layer scenarios with %-high titre####
scenarios.tmp <-  expand.grid(size = c(15000,
                                       32000
                                       ,64000), vac = c(0,50,60,70,80,90))

scenario.list.size.vaccination <- mapply(FUN = function(size,vac){list(data.frame(scenario  = paste0("layerSize",size,"Vac",vac)))}, scenarios.tmp$size,scenarios.tmp$vac)

output.layer <- lapply(c(1:length(scenario.list.size.vaccination)),function(i){readRDS(paste0("./output/coverage/20231210output",gsub(scenario.list.size.vaccination[[i]]$scenario,pattern = "[.]", replacement = ""),".RDS"))$out})
pars.layer <- lapply(c(1:length(scenario.list.size.vaccination)),function(i){readRDS(paste0("./output/coverage/20231210output",gsub(scenario.list.size.vaccination[[i]]$scenario,pattern = "[.]", replacement = ""),".RDS"))}$pars)

#plot the baseline size####
plot.data <- NULL;
for(j in c(1:length(scenario.list.size.vaccination))){
  plot.data<- rbind(plot.data,cbind(data.frame(output.layer[[j]]),
                                    scenario = scenario.list.size.vaccination[[j]]$scenario))
}
plot.data$size <- (plot.data$scenario%>%gsub(pattern = "layerSize", replacement = "")%>%str_split(pattern= c("Vac"))%>%unlist%>%matrix(nrow =  2))[1,]
plot.data$vac <- (plot.data$scenario%>%gsub(pattern = "layerSize", replacement = "")%>%str_split(pattern= c("Vac"))%>%unlist%>%matrix(nrow =  2))[2,]
sel.plot.data <- plot.data%>%select(c("time","run","scenario","size","vac","I.1","I.2","R.1","R.2","DR.1","DR.2"))%>%reshape2::melt(id.vars = c("time","run","scenario","size","vac"))
sel.plot.data <- (data.frame(sel.plot.data)%>%group_by(scenario)%>%sample_n(ceiling(0.05*length(time))))%>%ungroup

ggplot(sel.plot.data%>%filter(size == 32000))+
  geom_step(aes(x = time, y = value, colour = variable, group = run))+
  labs(x = "Time", y = "Number of birds")+
  facet_grid(vac~variable,labeller = labeller(variable = itype.labels) , scales= c("free_y"))
ggsave("./output/figures/baselineoutbreaklayer.png")


#surveillance ####
reps <- 10;
rm(surveillance.layer)
for(i in c(1:length(scenario.list.size.vaccination))){
  tmp <- cbind(repeat.detection.time.surveillance(as.data.frame(output.layer[[i]]),
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
               size = pars.layer[[i]]$N0,
               vaccination = pars.layer[[i]]$p.hightitre,
               introTime =0)
  surveillance.layer <-if(!exists("surveillance.layer")){tmp}else {rbind(surveillance.layer,tmp)}
  
}
#add for plotting and processing scenario name split into size and vaccination and add parameters for size and vaccination. 
surveillance.layer<- cbind(surveillance.layer,surveillance.layer$scenario%>%gsub(pattern ="layerSize", replacement = "")%>%str_split_fixed(pattern = c("Vac"),2)%>%as.data.frame())
surveillance.layer$vaccination <- unlist(str_split_fixed(surveillance.layer$scenario,pattern= c("Vac"),2))[,2]
surveillance.layer$size <-gsub(unlist(str_split_fixed(surveillance.layer$scenario,pattern= c("Vac"),2))[,1], pattern = "layerSize", replacement = "")
#save surveillance times 
saveRDS(surveillance.layer, file = "./output/surveillanceLayer.RDS")

surveillance.layer <- readRDS(file = "./output/surveillanceLayer.RDS")

#write code to show that size does not matter for min detection times

ggplot(data = surveillance.layer%>%
                  select(c("run","rep","size","vaccination","scenario","pas.det.time","ac.det.time","min.det.time"))%>%
         reshape2::melt(id.vars = c("run","rep","size","vaccination","scenario")) )+
   geom_histogram(aes(x = as.numeric(value), y=after_stat(count)/sum(after_stat(count)), fill = variable),    colour = "black", 
                           alpha = 0.5, binwidth = 1.0,position = "identity")+
  xlim(0,NA)+ylim(0,NA)+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "red",min.det.time= "orange"),
                    labels = detection.method.label)+
  xlab("Detection time")+
  ylab("Proportion of runs")+
#  ggtitle("Surveillance")+
   facet_nested(vaccination~size * variable, labeller = labeller(variable = detection.method.label) )+theme(legend.position = "none")
ggsave("./output/figures/scenarios_layer_surveillance.png", scale = 1.23)

#fit surveillance detection  
scenario.names<- surveillance.layer%>%select(scenario, size, vaccination)%>%unique
add.vars <- list(size = as.list(setNames( scenario.names$size,scenario.names$scenario)),
                 vac= as.list(setNames( scenario.names$vaccination,scenario.names$scenario)))


dists <- fit.distribution(surveillance.layer, added.vars = add.vars)
dists$size <- as.numeric(dists$size)
dists$vac <- as.numeric(dists$vac)
saveRDS(dists, "./output/detDist_sizevac.RDS")
dists <- readRDS("./output/detDist_sizevac.RDS")

fit.pas.mean.sizevac<-stats::lm(mean ~vac + size, data = dists%>%filter(detection.method=="pas.det.time"& mean >0))
drop1(fit.pas.mean.sizevac)
summary(fit.pas.mean.sizevac)

fit.pas.mean.sizevac<-stats::lm(mean ~vac , data = dists%>%filter(detection.method=="pas.det.time"& mean >0))
drop1(fit.pas.mean.sizevac)
summary(fit.pas.mean.sizevac)
saveRDS(fit.pas.mean.sizevac, "./output/detMean_sizevac_passive.RDS")
fit.pas.mean.sizevac<- readRDS( "./output/detMean_sizevac_passive.RDS")

fit.pas.sd.sizevac<-stats::lm(sqrt(var) ~vac + size, data = dists%>%filter(detection.method=="pas.det.time"& var >0))
drop1(fit.pas.sd.sizevac)
summary(fit.pas.sd.sizevac)

fit.sd.sizevac<-stats::lm(sqrt(var) ~vac , data = dists%>%filter(detection.method=="pas.det.time"& var >0))
drop1(fit.sd.sizevac)
summary(fit.sd.sizevac)
saveRDS(fit.pas.sizevac, "./output/detSD_sizevac_passive.RDS")


fit.min.mean.sizevac<-stats::lm(mean ~vac + size, data = dists%>%filter(detection.method=="min.det.time"& mean >0))
drop1(fit.min.mean.sizevac)
summary(fit.min.mean.sizevac)

fit.min.mean.sizevac<-stats::lm(mean ~vac , data = dists%>%filter(detection.method=="min.det.time"& mean >0))
drop1(fit.min.mean.sizevac)
summary(fit.min.mean.sizevac)
saveRDS(fit.min.mean.sizevac, "./output/detMean_sizevac_minimum.RDS")

fit.min.sd.sizevac<-stats::lm(sqrt(var) ~vac + size, data = dists%>%filter(detection.method=="min.det.time"& var >0))
drop1(fit.min.sd.sizevac)
summary(fit.min.sd.sizevac)

fit.sd.sizevac<-stats::lm(sqrt(var) ~vac , data = dists%>%filter(detection.method=="min.det.time"& var >0))
drop1(fit.sd.sizevac)
summary(fit.sd.sizevac)
saveRDS(fit.sd.sizevac, "./output/detSD_sizevac_minimum.RDS")



#human exposure ####
rm(humanexposure.layer)
#as we used multiple repeats to deterpase active surveillance repeat
humanexposure.layer  <- NULL

for(k in c(1:length(output.layer))){
  humanexposure.layer <- rbind(humanexposure.layer,
                                   cbind(human.exposure.detection.multiple.runs(output.layer[[k]],
                                                                            pars.layer[[k]]$beta,
                                                                            surveillance.layer%>%filter(scenario == pars.layer[[k]]$scenario)),
                                         scenario = pars.layer[[k]]$scenario,
                                         size = pars.layer[[k]]$N0,
                                         vaccination = pars.layer[[k]]$p.hightitre))
  
}

#calculate baseline exposure for each farm size
baseline.average <- humanexposure.layer%>%
  filter(vaccination == 0 & detection.method == "pas.det.time")%>%
  select(c("size","detection.exposure","total.exposure"))%>%reframe(.by = size,
                                                                                                                                           mean.det.exposure = mean(detection.exposure),
                                                                                                                                           mean.tot.exposure = mean(total.exposure))
#join with exposure 
humanexposure.layer <- left_join(x = humanexposure.layer, y = baseline.average, by = "size")

#calculate ratio
humanexposure.layer$ratio.det <- humanexposure.layer$detection.exposure/humanexposure.layer$mean.det.exposure
humanexposure.layer$ratio.tot <- humanexposure.layer$total.exposure/humanexposure.layer$mean.tot.exposure

saveRDS(humanexposure.layer, file = "./output/humanexposureLayer.RDS")
humanexposure.layer<- readRDS( file = "./output/humanexposureLayer.RDS")



ggplot(humanexposure.layer) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = detection.method),binwidth = .15,alpha = 0.5,colour = "black",position = "identity")+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  scale_x_log10()+
  #lims(x = c(0,3))+
  geom_vline(xintercept = 1)+
  facet_grid((as.numeric(vaccination)*100)~size)
  #ggtitle("Human exposure layers")#+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer.png")

# 
# saveRDS(humanexposure.layer%>%reframe(.by = c(scenario,detection.method),
#                                            mean = mean(detection.exposure)), file = "./output/MEANexposurefitLayerVaccin.RDS")
humanexposure.layer.dists <- fit.distribution(humanexposure.layer, 
                                              long.form = TRUE, long.var = "detection.exposure",
                                              scale = 1000, 
                                              added.vars = add.vars)

fit.exp.mean.pas <- lm(log10(mean) ~  as.numeric(size) +as.numeric(vac), 
          data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time"))

summary(fit.exp.mean.pas)
drop1(fit.exp.mean.pas)
predict(fit.exp.mean.pas, newdata = data.frame(size = 10, vac = 0.0))
predict(fit.exp.mean.pas, newdata = data.frame(size = 15000, vac = 0.0))
predict(fit.exp.mean.pas)
ggplot()+
  geom_point(aes(vac, log10(mean), colour =as.factor(size)),humanexposure.layer.dists%>%filter(detection.method == "pas.det.time"))+
  geom_path(aes(vac, val), data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "pas.det.time")%>%select("vac"), val =predict(fit.exp.mean.pas) ))

#requires another approach
#1. linear decrease from 0 to 50 % vaccination
#2. fixed value from 50 to 70% vaccination
#3 exponential decrease from 70 to 100% vaccination

fit.exp.mean.pas.1 <- lm(log10(mean) ~  as.numeric(size) +as.numeric(vac), 
                       data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) <= 50))

fit.exp.mean.pas.2 <- lm(log10(mean) ~  as.numeric(size) +as.numeric(vac), 
                         data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 50 & as.numeric(vac) <= 70))
fit.exp.mean.pas.2 <- lm(log10(mean) ~  as.numeric(size) , 
                         data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 50 & as.numeric(vac) <= 70))

fit.exp.mean.pas.3 <- lm(log10(mean) ~  as.numeric(size)+as.numeric(vac),  
                         data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 70 & as.numeric(vac) <= 90))

ggplot()+
  geom_point(aes(vac, log10(mean), colour =as.factor(size)),humanexposure.layer.dists%>%filter(detection.method == "pas.det.time"))+
  geom_path(aes(vac, val), data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) <= 50)%>%select("vac"), val =predict(fit.exp.mean.pas.1) ))+
  geom_path(aes(vac, val), data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 50 & as.numeric(vac) <= 70)%>%select("vac"), val =predict(fit.exp.mean.pas.2) ))+
  geom_path(aes(vac, val), data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 70 & as.numeric(vac) <= 90)%>%select("vac"), val =predict(fit.exp.mean.pas.3) ))



saveRDS(list(fit.exp.mean.pas.1,fit.exp.mean.pas.2,fit.exp.mean.pas.3), file = "./output/exposureMeanfitLayerPas.RDS")

#variance
fit.exp.var.pas <- lm(log10(var) ~  as.numeric(size) +as.numeric(vac), 
                       data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time"))

summary(fit.exp.var.pas)
drop1(fit.exp.var.pas)
ggplot()+
  geom_point(aes(vac, log10(var), colour =as.factor(size)),humanexposure.layer.dists%>%filter(detection.method == "pas.det.time"))+
  geom_path(aes(vac, val), data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "pas.det.time")%>%select("vac"), val =predict(fit.exp.var.pas) ))

#requires another approach
#1. linear decrease from 0 to 50 % vaccination
#2. fixed value from 50 to 70% vaccination
#3 exponential decrease from 70 to 100% vaccination

fit.exp.var.pas.1 <- lm(log10(var) ~  as.numeric(size) +as.numeric(vac), 
                         data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) <= 50))

fit.exp.var.pas.2 <- lm(log10(var) ~  as.numeric(size) +as.numeric(vac), 
                         data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 50 & as.numeric(vac) <= 70))

fit.exp.var.pas.3 <- lm(log10(var) ~  as.numeric(size)+as.numeric(vac),  
                         data = humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 70 & as.numeric(vac) <= 90))

ggplot()+
  geom_point(aes(vac, log10(var), colour =as.factor(size)),humanexposure.layer.dists%>%filter(detection.method == "pas.det.time"))+
  geom_path(aes(vac, val), rbind(data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) <= 50)%>%select("vac"), val =predict(fit.exp.var.pas.1)) ,
   data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 50 & as.numeric(vac) <= 70)%>%select("vac"), val =predict(fit.exp.var.pas.2) ),
  data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "pas.det.time", as.numeric(vac) >= 70 & as.numeric(vac) <= 90)%>%select("vac"), val =predict(fit.exp.var.pas.3) )))


saveRDS(list(fit.exp.var.pas.1,fit.exp.var.pas.2,fit.exp.var.pas.3), file = "./output/exposureVarfitLayerPas.RDS")

#minimum detection ####
fit.exp.mean.min <- lm(log10(mean) ~  as.numeric(size) +as.numeric(vac), 
                       data = humanexposure.layer.dists%>%filter(detection.method == "min.det.time"))

summary(fit.exp.mean.min)
drop1(fit.exp.mean.min)
ggplot()+
  geom_point(aes(vac, log10(mean), colour =as.factor(size)),humanexposure.layer.dists%>%filter(detection.method == "min.det.time"))+
  geom_path(aes(vac, val), data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "min.det.time")%>%select("vac"), val =predict(fit.exp.mean.min) ))


saveRDS(list(fit.exp.mean.min), file = "./output/exposureMeanfitLayerMin.RDS")

#variance
fit.exp.var.min <- lm(log10(var) ~  as.numeric(size) +as.numeric(vac), 
                      data = humanexposure.layer.dists%>%filter(detection.method == "min.det.time"))

summary(fit.exp.var.min)
drop1(fit.exp.var.min)
ggplot()+
  geom_point(aes(vac, log10(var), colour =as.factor(size)),humanexposure.layer.dists%>%filter(detection.method == "min.det.time"))+
  geom_path(aes(vac, val), data.frame(vac =humanexposure.layer.dists%>%filter(detection.method == "min.det.time")%>%select("vac"), val =predict(fit.exp.var.min) ))


saveRDS(fit.exp.var.min, file = "./output/exposureVarfitLayerMin.RDS")

#average ratios####
humanexposure.layer%>%reframe(.by = c(scenario,detection.method),
                                  mean = mean(ratio.det),
                                  min = min(ratio.det),
                                  max = max(ratio.det),
                                  perc25 = quantile(ratio.det,0.25),
                                  perc75 =quantile(ratio.det,0.75)
)





########################################
#Scenarios with clinical protection ####
########################################
#add the 32 000 animal farm without vaccination as baseline
scenario.clinprot.list <- list(list(scenario = "coverage/20231210outputlayerSize32000Vac0"),
                               list(scenario = "clinic/20231210outputlayerClinic50Phigh50"),
                               list(scenario = "clinic/20231210outputlayerClinic50Phigh0"),
                               list(scenario = "clinic/20231210outputlayerClinic0.1Phigh0"))
output.layer.clinprot <- lapply(c(1:4),function(i){as.data.frame(readRDS(paste0("./output/",scenario.clinprot.list[[i]]$scenario,".RDS"))$out)})
pars.layer.clinprot <- lapply(c(1:4),function(i){readRDS(paste0("./output/",scenario.clinprot.list[[i]]$scenario,".RDS"))$pars})

#visualize ####


plot.output.grid(rbind(cbind(output.layer.clinprot[[1]],scenario = pars.layer.clinprot[[1]]$scenario),
                       cbind(output.layer.clinprot[[2]],scenario = pars.layer.clinprot[[2]]$scenario),
                       cbind(output.layer.clinprot[[3]],scenario = pars.layer.clinprot[[3]]$scenario),
                       cbind(output.layer.clinprot[[4]],scenario = pars.layer.clinprot[[4]]$scenario)), 
                 vars =c("I.1","I.2","R.1","R.2","DR.1","DR.2"),itype.label = itype.labels,scenario.label = scenario.label.clin,
                 title = "Clinical protection")
ggsave(file = "./output/figures/clinicalprotectionLayer.png",scale = 1.23)

#surveillance ####
reps <- 10;

rm(surveillance.layer.clinprot)
for(i in c(1:length(scenario.clinprot.list))){
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
unique(surveillance.layer.clinprot$scenario)

saveRDS(surveillance.layer.clinprot,"./output/surveillancelayerclinprot.RDS")
surveillance.layer.clinprot<-readRDS("./output/surveillancelayerclinprot.RDS")

ggplot(data = surveillance.layer.clinprot%>%
         select(c("run","rep","scenario","pas.det.time","min.det.time","ac.det.time"))%>%
         reshape2::melt(id.vars = c("run","rep","scenario")) )+
  geom_histogram(aes(x = as.numeric(value), y=after_stat(count/sum(count)), fill = variable),    colour = "black", 
                 alpha = 0.5, binwidth = 1.0,position = "identity")+
  xlim(0,NA)+ylim(0,NA)+
  scale_fill_manual(values=c(pas.det.time = "black",min.det.time= "red",ac.det.time = "orange"),
                    labels = c(pas.det.time = "Passive",min.det.time= "Minimum",ac.det.time = "Active"))+
  xlab("Detection time")+
  ylab("Proportion of runs")+
#  ggtitle("Surveillance")+
  facet_grid(scenario~variable,labeller = labeller(scenario = scenario.label.clin, variable = detection.method.label))

ggsave("./output/figures/scenarios_layer.clinprot_surveillance.png")

dists.clinprot.surveillance <- fit.distribution(surveillance.layer.clinprot)

surveillance.layer.clinprot%>%
  select(c("run","rep","scenario","pas.det.time","ac.det.time","min.det.time"))%>%reframe(.by = scenario,
                                                                                          mean.pas = mean(pas.det.time),
                                                                                          mean.min = mean(min.det.time))

saveRDS(dists.clinprot.surveillance, "./output/distsDetClinProt.rds")
dists.clinprot.surveillance<-readRDS( "./output/distsDetClinProt.rds")

#human exposure Clinical protection ####
#as we used multiple repeats to determine active surveillance repeat###
humanexposure.layer.clin <- NULL

for(k in c(1:length(scenario.clinprot.list))){
  humanexposure.layer.clin <- rbind(humanexposure.layer.clin,
                             cbind(human.exposure.detection.multiple.runs(output.layer.clinprot[[k]],
                                                                      pars.layer.clinprot[[k]]$beta,
                                                                      surveillance.layer.clinprot%>%filter(scenario == pars.layer.clinprot[[k]]$scenario)), 
                                   scenario = pars.layer.clinprot[[k]]$scenario))
  
}

#calculate baseline exposure for each farm size
baseline.average.clin <- humanexposure.layer.clin%>%
  filter(scenario == "layerSize32000Vac0" & detection.method == "pas.det.time")%>%
  select(c("scenario","detection.exposure","total.exposure"))%>%reframe(
    .by = scenario,
                                                                    mean.det.exposure = mean(detection.exposure),
                                                                    mean.tot.exposure = mean(total.exposure))
#join with exposure 
humanexposure.layer.clin <- cbind(humanexposure.layer.clin, baseline.average.clin[,c(2,3)])

#calculate ratio
humanexposure.layer.clin$ratio.det <- humanexposure.layer.clin$detection.exposure/humanexposure.layer.clin$mean.det.exposure
humanexposure.layer.clin$ratio.tot <- humanexposure.layer.clin$total.exposure/humanexposure.layer.clin$mean.tot.exposure

saveRDS(humanexposure.layer.clin, file = "./output/humanexposureLayerClin.RDS")
humanexposure.layer.clin <- readRDS(file = "./output/humanexposureLayerClin.RDS")
dists.clinprot.exposure <- fit.distribution(humanexposure.layer.clin, 
                                            long.form= TRUE, 
                                            long.var = "detection.exposure",
                                            scale = 1000)

ggplot(humanexposure.layer.clin) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = detection.method),binwidth = .15,alpha = 0.5,colour = "black",position = "identity")+
  scale_fill_manual(values=c(pas.det.time = "black",min.det.time= "red",ac.det.time = "orange"),
                    labels = c(pas.det.time = "Passive",min.det.time= "Minimum",ac.det.time = "Active"))+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  lims(x = c(0,3))+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~detection.method,labeller = labeller(scenario = scenario.label.clin,
                                            detection.method = detection.method.label))
#+  ggtitle("Human exposure layers")#+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer_clinprot.png")



#average ratios

humanexposure.layer.clin%>%reframe(.by = c(scenario,detection.method),
                            mean = mean(ratio.det),
                            min = min(ratio.det),
                            max = max(ratio.det),
                            perc25 = quantile(ratio.det,0.25),
                            perc75 =quantile(ratio.det,0.75)
)
saveRDS(humanexposure.layer.clin%>%reframe(.by = c(scenario,detection.method),
                                           mean = mean(detection.exposure)), file = "./output/MEANexposurefitLayerVaccinClin.RDS")
dists.clinprot.exposure



###################################################
#Waning immunity - same strain           ####
###################################################
#add the 32 000 animal farm without vaccination as baseline
scenario.wane.list <- list(list(scenario = "layerSize32000Vac0"),
                           list(scenario = "layerStartTime0"),
                           list(scenario = "layerStartTime50"),
                           list(scenario = "layerStartTime100"),
                           list(scenario = "layerStartTime200"),
                           list(scenario = "layerStartTime300"),
                           list(scenario = "layerStartTime400"),
                           list(scenario = "layerStartTime500"))
scenario.names<- surveillance.layer.wane%>%select(scenario)%>%unique
phigh <- c(0,(1-pgamma(c(0,50,100,200,300,400,500),36,0.07)))
add.vars <- list(phigh= as.list(setNames( phigh,scenario.names$scenario)),
                 intro.time = as.list(setNames(c(0,0,50,100,200,300,400,500),scenario.names$scenario)))


output.layer.wane <- lapply(c(1:length(scenario.wane.list)),function(i){load.sims(paste0("./output/",gsub(scenario.wane.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)$output})
pars.layer.wane <- lapply(c(1:length(scenario.wane.list)),function(i){load.sims(paste0("./output/",gsub(scenario.wane.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)}$pars[[1]])



#visualize ##

plot.output.grid(rbind(cbind(output.layer.wane[[1]],scenario = pars.layer.wane[[1]]$scenario),
                       cbind(output.layer.wane[[2]],scenario = pars.layer.wane[[2]]$scenario),
                       cbind(output.layer.wane[[3]],scenario = pars.layer.wane[[3]]$scenario),
                       cbind(output.layer.wane[[4]],scenario = pars.layer.wane[[4]]$scenario),
                       cbind(output.layer.wane[[5]],scenario = pars.layer.wane[[5]]$scenario),
                       cbind(output.layer.wane[[6]],scenario = pars.layer.wane[[6]]$scenario),
                       cbind(output.layer.wane[[7]],scenario = pars.layer.wane[[7]]$scenario),
                       cbind(output.layer.wane[[8]],scenario = pars.layer.wane[[8]]$scenario)), 
                 vars =c("I.1","I.2","R.1","R.2","DR.1","DR.2"),
                 #title = "Waning Immunity", 
                 scales = "free_y" , itype.label = itype.labels, scenario.label = scenario.label.wane, scenario.levels = scenario.levels.label.wane )
ggsave(file = "./output/figures/waningLayer.png", scale =1.23)

#surveillance ###
reps <- 100;

rm(surveillance.layer.wane)
for(i in c(1:length(scenario.wane.list))){
  tmp <- cbind(repeat.detection.time.surveillance(output.layer.wane[[i]],
                                                  reps = reps,
                                                  deaths.vars  = c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                  time.interval.pas = 1,
                                                  threshold = 0.005*pars.layer.wane[[i]]$N0,
                                                  ints = 2,
                                                  detectables.vars = c("DI.1","DI.2","DR.1","DR.2"),
                                                  se = 0.99,
                                                  time.interval.ac =7,
                                                  init.ac ="rand",
                                                  detectables.incidence = TRUE,
                                                  pfarm = 1., 
                                                  panimal = 1.),
               scenario = pars.layer.wane[[i]]$scenario,
               introTime = pars.layer.wane[[i]]$intro.time)
  surveillance.layer.wane <-if(!exists("surveillance.layer.wane")){tmp}else {rbind(surveillance.layer.wane,tmp)}
  
}
surveillance.layer.wane<- cbind(surveillance.layer.wane,surveillance.layer.wane$scenario%>%gsub(pattern ="layer.wane_", replacement = "")%>%str_split_fixed(pattern = c("intro"),2)%>%as.data.frame())


saveRDS(surveillance.layer.wane,"./output/surveillancelayerwane.RDS")
surveillance.layer.wane<-readRDS("./output/surveillancelayerwane.RDS")
surveillance.layer.wane$scenario <- factor(surveillance.layer.wane$scenario, levels = scenario.levels.label.wane)

ggplot(data = surveillance.layer.wane%>%
         select(c("run","rep","scenario","pas.det.time","ac.det.time","min.det.time"))%>%
         reshape2::melt(id.vars = c("run","rep","scenario")) )+
  geom_histogram(aes(x = as.numeric(value), y=..count../sum(..count..), fill = variable),    colour = "black", 
                 alpha = 0.5, binwidth = 1.0,position = "identity")+
  xlim(0,NA)+ylim(0,NA)+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Detection time")+
  ylab("Proportion of runs")+
#  ggtitle("Surveillance")+
  facet_grid(scenario~variable,labeller = labeller(scenario = scenario.label.wane, 
                                                   variable = detection.method.label))


ggsave("./output/figures/scenarios_layer.wane_surveillance.png")

dists.wane <- fit.distribution(surveillance.layer.wane, added.vars = add.vars)

saveRDS(dists.wane, "./output/distsDetWane.rds")
dists.wane

surveillance.layer.wane%>%reframe(.by = scenario, mean.pas = mean(pas.det.time),mean.ac = mean(ac.det.time))

#human exposure Waning immunity####
#as we used multiple repeats to determine active surveillance repeat
humanexposure.layer.wane <- NULL

for(k in c(1:length(scenario.wane.list))){
  humanexposure.layer.wane <- rbind(humanexposure.layer.wane,
                                    cbind(human.exposure.detection.multiple.runs(output.layer.wane[[k]],
                                                                                 pars.layer.wane[[k]]$beta,
                                                                                 surveillance.layer.wane%>%filter(scenario == pars.layer.wane[[k]]$scenario)), 
                                          scenario = pars.layer.wane[[k]]$scenario,
                                          phigh = pars.layer.wane[[k]]$p.hightitre))
  
}

#calculate baseline exposure for each farm size
baseline.average.wane <- humanexposure.layer.wane%>%
  filter(scenario == "layerSize32000Vac0" & detection.method == "pas.det.time")%>%
  select(c("scenario","detection.exposure","total.exposure"))%>%reframe(
    .by = scenario,
    mean.det.exposure = mean(detection.exposure),
    mean.tot.exposure = mean(total.exposure))
#join with exposure 
humanexposure.layer.wane <- cbind(humanexposure.layer.wane, baseline.average.wane[,c("mean.det.exposure","mean.tot.exposure")])

#calculate ratio
humanexposure.layer.wane$ratio.det <- humanexposure.layer.wane$detection.exposure/humanexposure.layer.wane$mean.det.exposure
humanexposure.layer.wane$ratio.tot <- humanexposure.layer.wane$total.exposure/humanexposure.layer.wane$mean.tot.exposure

saveRDS(humanexposure.layer.wane, file = "./output/humanexposureLayerwane.RDS")
humanexposure.layer.wane<-readRDS(file = "./output/humanexposureLayerwane.RDS")

humanexposure.layer.wane$scenario <- factor(humanexposure.layer.wane$scenario, levels= scenario.levels.label.wane)
ggplot(humanexposure.layer.wane) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = detection.method),binwidth = .15,alpha = 0.5,colour = "black",position = "identity")+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,0,1))+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~detection.method,labeller = labeller(scenario = scenario.label.wane, detection.method = detection.method.label))+
  #ggtitle("Human exposure layers")+
  theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer_wane.png", dpi = 300, scale = 1.23)


#distributions human exposure####
dists.exposure.wane <- fit.distribution(humanexposure.layer.wane,
                                   long.form = TRUE, 
                                   long.var = "detection.exposure",
                                   scale = 1000, 
                                   added.vars = add.vars )
saveRDS(dists.exposure.wane,"./output/distsExposureWaneDifStrain.RDS")
dists.exposure.wane<- readRDS("./output/distsExposureWaneDifStrain.RDS")


fit.mean.exp.wane <- lm(mean ~  intro.time, data = dists.exposure.wane%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0"&intro.time > 300))
          
summary(fit.mean.exp.wane)
drop1(fit.mean.exp.wane)

ggplot()+
  geom_point(aes(intro.time, mean),dists.exposure.wane%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0"&intro.time > 0))+
  geom_line(aes(intro.time, val),color = "red", data.frame(intro.time=dists.exposure.wane%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0"&intro.time > 300)%>%select(intro.time), 
                                                           val =predict(fit.mean.exp.wane ) ))

saveRDS(fit.mean.exp.wane, file = "./output/exposurefitLayerVaccinWanePAS.RDS")
fit.mean.exp.wane<-readRDS(file = "./output/exposurefitLayerVaccinWanePAS.RDS")

fit.mean.exp.wane.min <- lm(mean ~  intro.time, data = dists.exposure.wane%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0"&intro.time <= 300))

summary(fit.mean.exp.wane.min)
drop1(fit.mean.exp.wane.min)

ggplot()+
  geom_point(aes(intro.time, mean),dists.exposure.wane%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0"&intro.time > 0))+
  geom_line(aes(intro.time, val),color = "red", data.frame(intro.time=dists.exposure.wane%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0"&intro.time > 0)%>%select(intro.time), 
                                                           val =predict(fit.mean.exp.wane.min ) ))


saveRDS(humanexposure.layer.wane%>%reframe(.by = c(scenario,detection.method),
                                           mean = mean(detection.exposure)), file = "./output/MEANexposurefitLayerVaccinWane.RDS")
saveRDS(fit, file = "./output/exposurefitLayerVaccinWaneMIN.RDS")
readRDS(file = "./output/exposurefitLayerVaccinWaneMIN.RDS")



#average ratios
humanexposure.layer.wane%>%reframe(.by = c(scenario,detection.method),
                                   mean = mean(ratio.det),
                                   min = min(ratio.det),
                                   max = max(ratio.det),
                                   perc25 = quantile(ratio.det,0.25),
                                   perc75 =quantile(ratio.det,0.75)
)


###################################################
#Waning immunity - different strain        ####
###################################################
#Waning immunity for a different strain####
#add the 32 000 animal farm without vaccination as baseline
scenario.waneDifStrain.list <- list(list(scenario = "layerSize32000Vac0"),
                           list(scenario = "layerStartTime0DifStrain"),
                           list(scenario = "layerStartTime50DifStrain"),
                           list(scenario = "layerStartTime100DifStrain"),
                           list(scenario = "layerStartTime200DifStrain"),
                           list(scenario = "layerStartTime300DifStrain"),
                           list(scenario = "layerStartTime400DifStrain"),
                           list(scenario = "layerStartTime500DifStrain"))
output.layer.waneDifStrain <- lapply(c(1:length(scenario.waneDifStrain.list)),function(i){load.sims(paste0("./output/",gsub(scenario.waneDifStrain.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)$output})
pars.layer.waneDifStrain <- lapply(c(1:length(scenario.waneDifStrain.list)),function(i){load.sims(paste0("./output/",gsub(scenario.waneDifStrain.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)}$pars[[1]])
pars.layer.waneDifStrain[[1]]$intro.time<-0



plot.output.grid(rbind(cbind(output.layer.waneDifStrain[[1]],scenario = pars.layer.waneDifStrain[[1]]$scenario),
                       cbind(output.layer.waneDifStrain[[2]],scenario = pars.layer.waneDifStrain[[2]]$scenario),
                       cbind(output.layer.waneDifStrain[[3]],scenario = pars.layer.waneDifStrain[[3]]$scenario),
                       cbind(output.layer.waneDifStrain[[4]],scenario = pars.layer.waneDifStrain[[4]]$scenario),
                       cbind(output.layer.waneDifStrain[[5]],scenario = pars.layer.waneDifStrain[[5]]$scenario),
                       cbind(output.layer.waneDifStrain[[6]],scenario = pars.layer.waneDifStrain[[6]]$scenario),
                       cbind(output.layer.waneDifStrain[[7]],scenario = pars.layer.waneDifStrain[[7]]$scenario),
                       cbind(output.layer.waneDifStrain[[8]],scenario = pars.layer.waneDifStrain[[8]]$scenario)), 
                 vars =c("I.1","I.2","R.1","R.2","DR.1","DR.2"),
                 #title = "Waning Immunity - Heterologous", 
                 scales = "free_y" , itype.label = itype.labels, scenario.label = scenario.label.waneDifStrain, scenario.levels = scenario.levels.label.waneDifStrain )
ggsave(file = "./output/figures/waningLayerDifStrain.png", scale =1.23)


#surveillance ####
reps <- 100;

rm(surveillance.layer.waneDifStrain)
for(i in c(1:length(scenario.waneDifStrain.list))){
  tmp <- cbind(repeat.detection.time.surveillance(output.layer.waneDifStrain[[i]],
                                                  reps = reps,
                                                  deaths.vars  = c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"),
                                                  time.interval.pas = 1,
                                                  threshold = 0.005*pars.layer.waneDifStrain[[i]]$N0,
                                                  ints = 2,
                                                  detectables.vars = c("DI.1","DI.2","DR.1","DR.2"),
                                                  se = 0.99,
                                                  time.interval.ac =7,
                                                  init.ac ="rand",
                                                  detectables.incidence = TRUE,
                                                  pfarm = 1., 
                                                  panimal = 1.),
               scenario = pars.layer.waneDifStrain[[i]]$scenario,
               introTime = pars.layer.waneDifStrain[[i]]$intro.time)
  surveillance.layer.waneDifStrain <-if(!exists("surveillance.layer.waneDifStrain")){tmp}else {rbind(surveillance.layer.waneDifStrain,tmp)}
  
}
surveillance.layer.waneDifStrain<- cbind(surveillance.layer.waneDifStrain,surveillance.layer.waneDifStrain$scenario%>%gsub(pattern ="layer.waneDifStrain_", replacement = "")%>%str_split_fixed(pattern = c("intro"),2)%>%as.data.frame())
surveillance.layer.waneDifStrain$scenario <- factor(surveillance.layer.waneDifStrain$scenario, levels = scenario.levels.label.waneDifStrain)
saveRDS(surveillance.layer.waneDifStrain,"./output/surveillancelayerwaneDifStrain.RDS")
#read saved outputs
surveillance.layer.waneDifStrain<-readRDS("./output/surveillancelayerwaneDifStrain.RDS")


ggplot(data = surveillance.layer.waneDifStrain%>%
         select(c("run","rep","scenario","pas.det.time","ac.det.time","min.det.time"))%>%
         reshape2::melt(id.vars = c("run","rep","scenario")) )+
  geom_histogram(aes(x = as.numeric(value), y=..count../sum(..count..), fill = variable),    colour = "black", 
                 alpha = 0.5, binwidth = 1.0,position = "identity")+
  xlim(0,NA)+ylim(0,NA)+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Detection time")+
  ylab("Proportion of runs")+
#  ggtitle("Surveillance")+
  facet_grid(scenario~variable,labeller = labeller(scenario = scenario.label.waneDifStrain, 
                                                   variable = detection.method.label))


ggsave("./output/figures/scenarios_layer.waneDifStrain_surveillance.png")

#fit distributions for detection ####
scenario.names<- surveillance.layer.waneDifStrain%>%select(scenario, introTime)%>%unique
phigh <- (1-pgamma(scenario.names$introTime,36,0.07))
phigh[1] <- 0 #baseline
add.vars <- list(phigh= as.list(setNames( phigh,scenario.names$scenario)),
                 intro.time = as.list(setNames(c(0,0,50,100,200,300,400,500),scenario.names$scenario)))

dists.wane.ds <- fit.distribution(surveillance.layer.waneDifStrain, added.vars = add.vars)
  
saveRDS(dists.wane.ds, "./output/detDist_waneds.RDS")
dists.wane.ds <- readRDS("./output/detDist_waneds.RDS")

fit.pas.mean.wane.ds<-stats::lm(mean ~ 0+phigh^150,offset = rep(9, nrow(dists.wane.ds%>%filter(detection.method=="pas.det.time"& mean >0))),  data = dists.wane.ds%>%filter(detection.method=="pas.det.time"& mean >0))
drop1(fit.pas.mean.wane.ds)
summary(fit.pas.mean.wane.ds)
ggplot()+
  geom_point(aes(phigh, mean),dists.wane.ds%>%filter(detection.method == "pas.det.time"& mean >0 & mean <20))+
  geom_point(aes(phigh, val),color = "red", data.frame(vac =dists.wane.ds%>%filter(detection.method == "pas.det.time"& mean >0)%>%select("phigh"), val =predict(fit.pas.mean.wane.ds) ))+
  geom_line(aes(phigh, phigh.1),color = "blue", data= data.frame(vac = c(dists.wane.ds%>%filter(detection.method == "pas.det.time"& mean >0)%>%select("phigh")), 
                                                        val = c(9+(42-9)*(dists.wane.ds%>%filter(detection.method == "pas.det.time"& mean >0)%>%select("phigh"))^500)) )


ggplot()+
  geom_point(aes(intro.time, mean),dists.wane.ds%>%filter(detection.method == "pas.det.time" & scenario != "layerSize32000Vac0" ))#+
  #geom_point(aes(intro.time, val),color = "red", data.frame(vac =dists.wane.ds%>%filter(detection.method == "pas.det.time"& mean >0)%>%select("phigh"), val =predict(fit.pas.mean.wane.ds) ))#+
#  geom_line(aes(intro.time, phigh.1),color = "blue", data= data.frame(vac = c(dists.wane.ds%>%filter(detection.method == "pas.det.time"& mean >0)%>%select("phigh")), 
 #                                                                     val = c(9+(42-9)*(dists.wane.ds%>%filter(detection.method == "pas.det.time"& mean >0)%>%select("phigh"))^500)) )         

fit.pas.mean.wane.ds<-stats::lm(log10(mean) ~ intro.time,  
                                data = dists.wane.ds%>%filter(detection.method=="pas.det.time"& intro.time >50))
  
drop1(fit.pas.mean.wane.ds)
summary(fit.pas.mean.wane.ds)

ggplot()+
  geom_point(aes(intro.time, mean),dists.wane.ds%>%filter(detection.method == "pas.det.time" &  intro.time >50 ))+
  geom_line(aes(intro.time, 10^val),color = "red", data.frame(vac =dists.wane.ds%>%filter(detection.method == "pas.det.time"& intro.time >50)%>%select("intro.time"), 
                                                                  val =predict(fit.pas.mean.wane.ds) ))
saveRDS(fit.pas.mean.wane.ds, "./output/detMean_wane.ds_passive.RDS")


fit.pas.sd.wane.ds<-stats::lm(sqrt(var) ~intro.time, data = dists.wane.ds%>%filter(detection.method=="pas.det.time"& var >0))
drop1(fit.pas.sd.wane.ds)
summary(fit.pas.sd.wane.ds)

ggplot()+
  geom_point(aes(intro.time, sqrt(var)),dists.wane.ds%>%filter(detection.method == "pas.det.time" & intro.time >50 ))+
  geom_line(aes(intro.time, val),color = "red", data.frame(vac =dists.wane.ds%>%filter(detection.method == "pas.det.time"& intro.time >50)%>%select("intro.time"), 
                                                              val =predict(fit.pas.sd.wane.ds) ))

drop1(fit.pas.sd.wane.ds)
summary(fit.pas.sd.wane.ds)
saveRDS(fit.pas.sd.wane.ds, "./output/detSD_wane.ds_passive.RDS")

fit.min.mean.wane.ds<-stats::lm(mean ~ intro.time,  
                                data = dists.wane.ds%>%filter(detection.method=="min.det.time"& intro.time >50))

drop1(fit.min.mean.wane.ds)
summary(fit.min.mean.wane.ds)

ggplot()+
  geom_point(aes(intro.time,mean),dists.wane.ds%>%filter(detection.method == "min.det.time" &  intro.time >50 ))+
  geom_line(aes(intro.time, val),color = "red", data.frame(vac =dists.wane.ds%>%filter(detection.method == "min.det.time"& intro.time >50)%>%select("intro.time"), 
                                                              val =predict(fit.min.mean.wane.ds) ))
saveRDS(fit.min.mean.wane.ds, "./output/detMean_wane.ds_minimum.RDS")


fit.min.sd.wane.ds<-stats::lm(sqrt(var) ~intro.time, data = dists.wane.ds%>%filter(detection.method=="min.det.time"& intro.time >50))
drop1(fit.min.sd.wane.ds)
summary(fit.min.sd.wane.ds)

ggplot()+
  geom_point(aes(intro.time, sqrt(var)),dists.wane.ds%>%filter(detection.method == "min.det.time" & intro.time >50 ))+
  geom_line(aes(intro.time, val),color = "red", data.frame(vac =dists.wane.ds%>%filter(detection.method == "min.det.time"& intro.time >50)%>%select("intro.time"), 
                                                           val =predict(fit.min.sd.wane.ds) ))

drop1(fit.min.sd.wane.ds)
summary(fit.min.sd.wane.ds)
saveRDS(fit.min.sd.wane.ds, "./output/detSD_wane.ds_minimum.RDS")



#human exposure Waning immunity####
#as we used multiple repeats to determine active surveillance repeat
humanexposure.layer.waneDifStrain <- NULL

for(k in c(1:length(scenario.waneDifStrain.list))){
  humanexposure.layer.waneDifStrain <- rbind(humanexposure.layer.waneDifStrain,
                                    cbind(human.exposure.detection.multiple.runs(output.layer.waneDifStrain[[k]],
                                                                                 pars.layer.waneDifStrain[[k]]$beta,
                                                                                 surveillance.layer.waneDifStrain%>%filter(scenario == pars.layer.waneDifStrain[[k]]$scenario)), 
                                          scenario = pars.layer.waneDifStrain[[k]]$scenario,
                                          phigh = pars.layer.waneDifStrain[[k]]$p.hightitre))
  
}

#calculate baseline exposure for each farm size
baseline.average.waneDifStrain <- humanexposure.layer.waneDifStrain%>%
  filter(scenario == "layerSize32000Vac0" & detection.method == "pas.det.time")%>%
  select(c("scenario","detection.exposure","total.exposure"))%>%reframe(
    .by = scenario,
    mean.det.exposure = mean(detection.exposure),
    mean.tot.exposure = mean(total.exposure))
#join with exposure 
humanexposure.layer.waneDifStrain <- cbind(humanexposure.layer.waneDifStrain, baseline.average.waneDifStrain[,c("mean.det.exposure","mean.tot.exposure")])

#calculate ratio
humanexposure.layer.waneDifStrain$ratio.det <- humanexposure.layer.waneDifStrain$detection.exposure/humanexposure.layer.waneDifStrain$mean.det.exposure
humanexposure.layer.waneDifStrain$ratio.tot <- humanexposure.layer.waneDifStrain$total.exposure/humanexposure.layer.waneDifStrain$mean.tot.exposure
humanexposure.layer.waneDifStrain$scenario <- factor(humanexposure.layer.waneDifStrain$scenario, levels= scenario.levels.label.waneDifStrain)
saveRDS(humanexposure.layer.waneDifStrain, file = "./output/humanexposureLayerwaneDifStrain.RDS")
humanexposure.layer.waneDifStrain<-readRDS(file = "./output/humanexposureLayerwaneDifStrain.RDS")


ggplot(humanexposure.layer.waneDifStrain) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = detection.method),
                 #binwidth = .15,
                 alpha = 0.5,colour = "black",position = "identity")+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  scale_x_log10(breaks=c(0.0001,0.001,0.01,0.1,0,1))+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~detection.method,labeller = labeller(scenario = scenario.label.waneDifStrain, detection.method = detection.method.label))+
#  ggtitle("Human exposure layers")+
  theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer_waneDifStrain.png", dpi = 300, scale = 1.23)

humanexposure.layer.waneDifStrain$fade.out <-humanexposure.layer.waneDifStrain$phigh >0.5

ggplot(humanexposure.layer.waneDifStrain%>%filter(detection.method =="min.det.time")%>% select(-"rep") %>%unique) + 
  geom_point(aes(x = detection.time, y =log10(detection.exposure), colour = scenario ))+facet_grid(scenario~., scales = "free")
  
  
#distributions human exposure####
dists.exposure <- fit.distribution(humanexposure.layer.waneDifStrain,
                                    long.form = TRUE, 
                                    long.var = "detection.exposure",
                                    scale = 1000, 
                                    added.vars = add.vars )
saveRDS(dists.exposure,"./output/distsExposureWaneDifStrain.RDS")
dists.exposure<- readRDS("./output/distsExposureWaneDifStrain.RDS")

fit.pas.exp.mean.wane.ds <- lm(log10(mean) ~  intro.time,
          data = dists.exposure%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0"))
summary(fit.pas.exp.mean.wane.ds )
drop1(fit.pas.exp.mean.wane.ds )
ggplot()+
  geom_point(aes(intro.time, mean),dists.exposure%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0"))+
  geom_line(aes(intro.time, 10^(val)),color = "red", data.frame(intro.time=dists.exposure%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0")%>%select(intro.time), 
                                                           val =predict(fit.pas.exp.mean.wane.ds ) ))
saveRDS(fit.pas.exp.mean.wane.ds , file = "./output/exposurefitLayerVaccinwaneDifStrainPAS.RDS")
fit.pas.exp.mean.wane.ds <- readRDS(file = "./output/exposurefitLayerVaccinwaneDifStrainPAS.RDS")

fit.pas.exp.var.wane.ds <- lm(sqrt(var) ~  intro.time,
                               data = dists.exposure%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0"))
summary(fit.pas.exp.var.wane.ds )
drop1(fit.pas.exp.var.wane.ds )
ggplot()+
  geom_point(aes(intro.time, sqrt(var)),dists.exposure%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0"))+
  geom_line(aes(intro.time, val),color = "red", data.frame(intro.time=dists.exposure%>%filter(detection.method=="pas.det.time"&scenario != "layerSize32000Vac0")%>%select(intro.time), 
                                                                val =predict(fit.pas.exp.var.wane.ds ) ))
saveRDS(fit.pas.exp.var.wane.ds, file = "./output/exposurefitLayerVaccinwaneDifStrainSDPAS.RDS")
fit.pas.exp.var.wane.ds<- readRDS(file = "./output/exposurefitLayerVaccinwaneDifStrainSDPAS.RDS")

fit.min.exp.mean.wane.ds <- lm(mean ~  intro.time,
                               data = dists.exposure%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0" &intro.time > 50))
summary(fit.min.exp.mean.wane.ds )
drop1(fit.min.exp.mean.wane.ds )
ggplot()+
  geom_point(aes(intro.time, mean),dists.exposure%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0"&intro.time > 50))+
  geom_line(aes(intro.time, val),color = "red", data.frame(intro.time=dists.exposure%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0"&intro.time > 50)%>%select(intro.time), 
                                                                val =predict(fit.min.exp.mean.wane.ds ) ))
saveRDS(fit.min.exp.mean.wane.ds, file = "./output/exposurefitLayerVaccinwaneDifStrainmin.RDS")


fit.min.exp.var.wane.ds <- lm(sqrt(var) ~  intro.time,
                              data = dists.exposure%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0"&intro.time > 50))
summary(fit.min.exp.var.wane.ds )
drop1(fit.min.exp.var.wane.ds )
ggplot()+
  geom_point(aes(intro.time, sqrt(var)),dists.exposure%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0"&intro.time > 50))+
  geom_line(aes(intro.time, val),color = "red", data.frame(intro.time=dists.exposure%>%filter(detection.method=="min.det.time"&scenario != "layerSize32000Vac0"&intro.time > 50)%>%select(intro.time), 
                                                                val =predict(fit.min.exp.var.wane.ds ) ))
saveRDS(fit.min.exp.var.wane.ds, file = "./output/exposurefitLayerVaccinwaneDifStrainSDmin.RDS")









#average ratios
humanexposure.layer.waneDifStrain%>%reframe(.by = c(scenario,detection.method),
                                   mean = mean(ratio.det),
                                   min = min(ratio.det),
                                   max = max(ratio.det),
                                   perc25 = quantile(ratio.det,0.25),
                                   perc75 =quantile(ratio.det,0.75)
)
