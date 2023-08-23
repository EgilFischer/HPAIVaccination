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
output.baseline.layer <- cbind(load.sims("./output/layerSize32000Vac0")$output)
output.baseline.broiler <- load.sims("./output/broilerSize38000Vac0")$output

#load baseline parameters (assuming all parameter lists are equal within a folder)
pars.baseline.layer <- load.sims("./output/layerSize32000Vac0", params = TRUE)$pars[[1]]
pars.baseline.broiler <- load.sims("./output/broilerSize38000Vac0", params = TRUE)$pars[[1]]

#probability of a major outbreak ####
q1q2.baseline.layer <- q1q2(pars.baseline.layer)
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
ggsave("./output/figures/probMinorOutbreakLayerAnno.png")


#visualize baseline #
plot.output.grid(output = rbind(cbind(output.baseline.layer,scenario = "Layer")#,
                       #cbind(output.baseline.broiler,scenario = "Broiler")
                       ),
                 vars = c("I.1","I.2","DR.1","DR.2"), title = "Baseline", frac = 0.1)
ggsave("./output/figures/baselineLayer.png")




#surveillance #
reps <- 100;
system.time(surveillance.layer.baseline <- cbind(repeat.detection.time.surveillance(output.baseline.layer,
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
                                     scenario = pars.baseline.layer$scenario))



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







#Layer and Broiler different sizes ####
scenarios.tmp <-  data.frame(
  type = c(rep("layer",3),rep("broiler",3)),
  size = c(15000,32000,64000,
           20000,38000,73000))

scenario.list.size.type <- mapply(FUN = function(type, size){list(data.frame(scenario  = paste0(type,"Size",size,"Vac0")))}, scenarios.tmp$type,scenarios.tmp$size)

output.size.type <- lapply(c(1:length(scenario.list.size.type)),function(i){load.sims(paste0("./output/",gsub(scenario.list.size.type[[i]]$scenario,pattern = "[.]", replacement = "")))$output})
pars.size.type <- lapply(c(1:length(scenario.list.size.type)),function(i){load.sims(paste0("./output/",gsub(scenario.list.size.type[[i]]$scenario,pattern = "[.]", replacement = "")))}$pars[[1]])

#plot the baseline size
plot.data <- NULL;
for(j in c(1:length(scenario.list.size.type))){
  tmp <- cbind(data.frame(output.size.type[[j]]),
               scenario = scenario.list.size.type[[j]]$scenario,
               size = scenarios.tmp$size[[j]],
               type = scenarios.tmp$type[[j]])
  tmp <- sample_n(tmp,ceiling(0.01*length(tmp$time)))
  plot.data<- rbind(plot.data,tmp)
}
lf.plot.data<- plot.data%>%reshape2::melt(id.vars = c("time","run","dt","scenario", "size","type") )
inf.lf.plot.data<- lf.plot.data%>%filter(variable %in% c("I.1","I.2","R.1","R.2"))

ggplot(data = lf.plot.data%>%filter(variable %in% c("S.1","I.1","R.1")))+
  geom_step(aes(x = time, y = value, colour = type,group = as.factor(run)))+
  scale_colour_manual(values = c("darkblue","darkred"),labels = c("Broiler","Layer"),name = "Type")+
  facet_grid(as.factor(size) ~ variable) 

ggsave("./output/figures/baselineoutbreaklayerbroiler.png")

#detection of baseline scenarios
reps <- 100;
rm(surveillance.size.type)
for(i in c(1:length(scenario.list.size.type))){
  tmp <- cbind(repeat.detection.time.surveillance(output.size.type[[i]],
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
  surveillance.size.type <-if(!exists("surveillance.size.type")){tmp}else {rbind(surveillance.size.type,tmp)}
  
}


#save surveillance times 
saveRDS(surveillance.layer, file = "./output/surveillanceSizeType.RDS")
#surveillance.size.type <- readRDS(file = "./output/surveillanceSizeType.RDS")


#write code to show that size does not matter for min detetion times
labels.type <-c("46"  = "Broiler", "510" = "Layer")
labels.size <-c("15000"  = "Small", "32000" = "Medium", "64000" = "Large", "20000"  = "Small", "38000" = "Medium", "73000" = "Large")
ggplot(data = surveillance.size.type%>%
         select(c("run","rep","size","max.time","scenario","pas.det.time","ac.det.time","min.det.time"))%>%
         reshape2::melt(id.vars = c("run","rep","size","max.time","scenario")) )+
  geom_histogram(aes(x = as.numeric(value), y=after_stat(count)/sum(after_stat(count)), fill = variable),    colour = "black", 
                 alpha = 0.5, binwidth = 1.0,position = "identity")+
  xlim(0,NA)+ylim(0,NA)+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  ggtitle("Surveillance")+
  facet_grid(size~max.time, labeller = labeller(max.time = labels.type))
ggsave("./output/figures/scenarios_sizetype_surveillance.png")

surveillance.size.type%>%reframe(.by = scenario,
                                 mean = mean(pas.det.time), 
                                 var = var(pas.det.time),
                                 q1 = quantile(pas.det.time,c(0.025)), 
                                 q2 = quantile(pas.det.time,c(0.975)))

#fit distributions ####
dists <- data.frame(shape =c(), rate =c(), scenario =c())
for(its in unique(surveillance.size.type$scenario))  {
  #for passive detection only use first repetition otherwise the data is inflated
  x = surveillance.size.type%>%
    filter(rep == 1)%>%
    filter(scenario == its)%>%dplyr::select("pas.det.time")%>%filter(pas.det.time!=Inf)
  
  
  if(length(x$pas.det.time)>1)
  {
    fit.gamma <- fitdistrplus::fitdist(data = x$pas.det.time, distr = "gamma", method = "mle")
    sgamma <-summary(fit.gamma)
    dists<- rbind(dists,c(shape = as.numeric(sgamma$estimate[1]),rate = as.numeric(sgamma$estimate[2]), scenario = its, pdetect = length(x$pas.det.time)/10))
    #plot(fit.gamma)
  }else dists <- rbind(dists,c(shape =c(0), rate =c(10^-5), scenario = its, pdetect = length(x$pas.det.time)/10))
}
names(dists)<-c("shape","rate","scenario","pdetect")
dists$shape <- as.numeric(dists$shape)
dists$rate <- as.numeric(dists$rate)
dists$mean <- dists$shape/dists$rate
dists$var <- dists$shape/dists$rate^2

dists


saveRDS(dists[,c("scenario","rate")], "./output/pasdetRate_sizetype.RDS")
saveRDS(dists[,c("scenario","shape")], "./output/pasdetShape_sizetype.RDS")



#Layer scenarios with %-high titre####
scenarios.tmp <-  expand.grid(size = c(15000,
                                       32000
                                       ,64000), vac = c(0,50,60,70,80,90))
#scenarios.tmp <- scenarios.tmp[c(1:15,18),]
scenario.list.size.vaccination <- mapply(FUN = function(size,vac){list(data.frame(scenario  = paste0("layerSize",size,"Vac",vac)))}, scenarios.tmp$size,scenarios.tmp$vac)

output.layer <- lapply(c(1:length(scenario.list.size.vaccination)),function(i){load.sims(paste0("./output/",gsub(scenario.list.size.vaccination[[i]]$scenario,pattern = "[.]", replacement = "")))$output})
pars.layer <- lapply(c(1:length(scenario.list.size.vaccination)),function(i){load.sims(paste0("./output/",gsub(scenario.list.size.vaccination[[i]]$scenario,pattern = "[.]", replacement = "")))}$pars[[1]])

#visualize ####
for(i in c(1:length(scenario.list.size.vaccination))){
  show(plot.output.sparse(output.layer[[i]],c("I.1","I.2","R.1","R.2"), scenario.list.size.vaccination[[i]]$scenario, frac= 0.05))
  ggsave(paste0("./output/figures/", gsub(scenario.list.size.vaccination[[i]]$scenario,pattern = "[.]", replacement = ""),".png"))
  gc()  
}
#plot the baseline size
plot.data <- NULL;
for(j in c(1:length(scenario.list.size.vaccination))){#} c(2,5,8,11,14,17)){
  plot.data<- rbind(plot.data,cbind(data.frame(output.layer[[j]]),scenario = scenario.list.size.vaccination[[j]]$scenario))
}
plot.data$size <- (plot.data$scenario%>%gsub(pattern = "layerSize", replacement = "")%>%str_split(pattern= c("Vac"))%>%unlist%>%matrix(nrow =  2))[1,]
plot.data$vac <- (plot.data$scenario%>%gsub(pattern = "layerSize", replacement = "")%>%str_split(pattern= c("Vac"))%>%unlist%>%matrix(nrow =  2))[2,]
sel.plot.data <- plot.data%>%select(c("time","run","scenario","size","vac","I.1","I.2","R.1","R.2"))%>%reshape2::melt(id.vars = c("time","run","scenario","size","vac"))
sel.plot.data <- data.frame(sel.plot.data)%>%sample_n(ceiling(0.01*length(sel.plot.data$time)))
ggplot(sel.plot.data)+
  geom_step(aes(x = time, y = value, colour = variable))+
  labs(x = "Time", y = "Number of birds")+
  facet_grid(vac~size)

pght <- plot.output.grid(plot.data,vars =c("I.1","I.2","R.1","R.2") ,title = "Layer",frac = 0.05)+
  theme(legend.position = "none")
pght
ggsave("./output/figures/baselineoutbreaklayer.png", pght)

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
#surveillance.layer <- readRDS(file = "./output/surveillanceLayer32000.RDS")
#surveillance.layer <- readRDS(file = "./output/surveillanceLayer.RDS")

#write code to show that size does not matter for min detetion times


ggplot(data = surveillance.layer%>%
                  select(c("run","rep","size","vaccination","scenario","pas.det.time","ac.det.time","min.det.time"))%>%
         reshape2::melt(id.vars = c("run","rep","size","vaccination","scenario")) )+
   geom_histogram(aes(x = as.numeric(value), y=after_stat(count)/sum(after_stat(count)), fill = variable),    colour = "black", 
                           alpha = 0.5, binwidth = 1.0,position = "identity")+
  xlim(0,NA)+ylim(0,NA)+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "red",min.det.time= "orange"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  ggtitle("Surveillance")+
   facet_grid(vaccination~size)
ggsave("./output/figures/scenarios_layer_surveillance.png")


dists <- data.frame(shape =c(), rate =c(), scenario =c())
for(its in unique(surveillance.layer$scenario))  {
  #for passive detection only use first repetition otherwise the data is inflated
  x = surveillance.layer%>%
    filter(rep == 1)%>%
    filter(scenario == its)%>%dplyr::select("pas.det.time")%>%filter(pas.det.time!=Inf)
  
  
  if(length(x$pas.det.time)>1)
  {
  fit.gamma <- fitdistrplus::fitdist(data = x$pas.det.time, distr = "gamma", method = "mle")
  sgamma <-summary(fit.gamma)
  dists<- rbind(dists,c(shape = as.numeric(sgamma$estimate[1]),rate = as.numeric(sgamma$estimate[2]), scenario = its, pdetect = length(x$pas.det.time)/10))
  #plot(fit.gamma)
  }else dists <- rbind(dists,c(shape =c(0), rate =c(10^-5), scenario = its, pdetect = length(x$pas.det.time)/10))
}
names(dists)<-c("shape","rate","scenario","pdetect")
dists$shape <- as.numeric(dists$shape)
dists$rate <- as.numeric(dists$rate)
dists$mean <- dists$shape/dists$rate
dists$var <- dists$shape/dists$rate^2
dists$vac <- gsub(x = dists$scenario,pattern = "layerSize",replacement = "")
dists$size <- as.numeric(sapply(dists$vac, function(x){str_split_1(x, pattern=c("Vac"))[1]}))
dists$vac <- as.numeric(sapply(dists$vac, function(x){str_split_1(x, pattern=c("Vac"))[2]}))
dists<- rbind(dists, last(dists))
dists[length(dists$shape),"vac"] <- 100
dists
dists

saveRDS(dists[,c("scenario","rate")], "./output/pasdetRate.RDS")
saveRDS(dists[,c("scenario","shape")], "./output/pasdetShape.RDS")


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
humanexposure.layer%>%filter(rep == 1 & detection.method=="pas.det.time"& is.finite(detection.time))%>%ggplot()+geom_point(aes(detection.time,detection.exposure,shape = as.factor(vaccination), colour = as.factor(size)))

fit <- lm(detection.exposure ~  size , 
          data = humanexposure.layer%>%filter(rep == 1 & detection.method=="pas.det.time"& is.finite(detection.time)& vaccination ==0))
summary(fit)
drop1(fit)
predict(fit, newdata = data.frame(size = 2000, vaccination = 0.5))
saveRDS(fit, file = "./output/exposurefitLayer.RDS")

fit <- lm(detection.exposure ~  size + vaccination, 
          data = humanexposure.layer%>%filter(rep == 1 & detection.method=="pas.det.time"& is.finite(detection.time)& vaccination >0))
summary(fit)
drop1(fit)
predict(fit, newdata = data.frame(size = 1000, vaccination = 0.5))
saveRDS(fit, file = "./output/exposurefitLayerVaccin.RDS")

#only passive
ggplot(humanexposure.layer%>%filter(detection.method =="pas.det.time")) +
  geom_histogram(aes(log10(ratio.det),after_stat(.15*density), fill = detection.method),binwidth = .15,alpha = 0.5,colour = "black",position = "identity")+
  labs(x = bquote(log[10] ~ 'Risk ratio of exposure  (reference = 0)')
       , y = "Proportion", fill = "Detection method")+
  geom_vline(xintercept = 0)+
  facet_grid(vaccination~size)+
             #,labeller = labeller(scenario = scenario.label))+
  ggtitle("Human exposure layers")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer_pas.png")

ggplot(data =humanexposure.layer) +
  geom_histogram( aes(log10(ratio.det),after_stat(.1*density), fill = detection.method),binwidth = .1,colour = "black",alpha = 0.5,position = "identity")+
  labs(x = bquote(log[10] ~ 'Risk ratio of exposure  (reference = 0)')
       , y = "Proportion", fill = "Detection method")+
  scale_fill_manual(values = c(min.det.time = "red",ac.det.time = "orange",pas.det.time = "black"),
                    labels = c(min.det.time ="Minimum",ac.det.time ="Active",pas.det.time ="Passive"))+
  geom_vline(xintercept = log10(1))+
  facet_grid(vaccination~size)+
  ggtitle("Human exposure")#+theme(legend.title =  element_text('Prop. high titre'))
ggsave("./output/figures/humanexposurelayer.png")

#average ratios
humanexposure.layer%>%reframe(.by = c(scenario,detection.method),
                                  mean = mean(ratio.det),
                                  min = min(ratio.det),
                                  max = max(ratio.det),
                                  perc25 = quantile(ratio.det,0.25),
                                  perc75 =quantile(ratio.det,0.75)
)



#Scenarios with clinical protection ####
#add the 32 000 animal farm without vaccination as baseline
scenario.clinprot.list <- list(list(scenario = "layerSize32000Vac0"),list(scenario = "layerClinic50Phigh50"),list(scenario = "layerClinic50Phigh0"),list(scenario = "layerClinic01Phigh0"))
output.layer.clinprot <- lapply(c(1:4),function(i){load.sims(paste0("./output/",gsub(scenario.clinprot.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)$output})
pars.layer.clinprot <- lapply(c(1:4),function(i){load.sims(paste0("./output/",gsub(scenario.clinprot.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)}$pars[[1]])



#visualize ##
scenario.label.clin <- c(layerSize32000Vac0 = "Baseline",
                    layerClinic50Phigh50 = "Clinic 50% \n High titre 50%",
                    layerClinic50Phigh0 = "Clinic 50% \n High titre 0%",
                    layerClinic0.1Phigh0 = "Clinic 0.1% \n High titre 0%")

detection.method.label <- c(pas.det.time = "Passive",
                            ac.det.time = "Active",
                            min.det.time = "Minimum")

for(i in c(1:length(scenario.clinprot.list))){
  show(plot.output(output.layer.clinprot[[i]],c("I.1","I.2","R.1","R.2"), scenario.clinprot.list[[i]]$scenario))
  ggsave(paste0("./output/figures/", gsub(scenario.clinprot.list[[i]]$scenario.clinprot,pattern = "[.]", replacement = ""),".png"))
  gc()  
}

plot.output.grid(rbind(cbind(output.layer.clinprot[[1]],scenario = pars.layer.clinprot[[1]]$scenario),
                       cbind(output.layer.clinprot[[2]],scenario = pars.layer.clinprot[[2]]$scenario),
                       cbind(output.layer.clinprot[[3]],scenario = pars.layer.clinprot[[3]]$scenario),
                       cbind(output.layer.clinprot[[4]],scenario = pars.layer.clinprot[[4]]$scenario)), 
                 vars =c("I.1","I.2","R.1","R.2"),
                 title = "Clinical protection")
ggsave(file = "./output/figures/clinicalprotectionLayer.png")

#surveillance ###
reps <- 100;

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


saveRDS(surveillance.layer.clinprot,"./output/surveillancelayerclinprot.RDS")
surveillance.layer.clinprot<-readRDS("./output/surveillancelayerclinprot.RDS")

ggplot(data = surveillance.layer.clinprot%>%
         select(c("run","rep","scenario","pas.det.time","ac.det.time","min.det.time"))%>%
         reshape2::melt(id.vars = c("run","rep","scenario")) )+
  geom_histogram(aes(x = as.numeric(value), y=after_stat(count/sum(count)), fill = variable),    colour = "black", 
                 alpha = 0.5, binwidth = 1.0,position = "identity")+
  xlim(0,NA)+ylim(0,NA)+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "red",min.det.time= "orange"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Detection time")+
  ylab("Proportion of runs")+
  ggtitle("Surveillance")+
  facet_grid(scenario~.,labeller = labeller(scenario = scenario.label.clin))

ggsave("./output/figures/scenarios_layer.clinprot_surveillance.png")

dists <- data.frame(shape =c(), rate =c(), scenario =c())
for(its in unique(surveillance.layer.clinprot$scenario))  {
  #for passive detection only use first repetition otherwise the data is inflated
  x = surveillance.layer.clinprot%>%
    filter(rep == 1)%>%
    filter(scenario == its)%>%dplyr::select("pas.det.time")%>%filter(pas.det.time!=Inf)
  
  
  if(length(x$pas.det.time)>1)
  {
    fit.gamma <- fitdistrplus::fitdist(data = x$pas.det.time, distr = "gamma", method = "mle")
    sgamma <-summary(fit.gamma)
    dists<- rbind(dists,c(shape = as.numeric(sgamma$estimate[1]),rate = as.numeric(sgamma$estimate[2]), scenario = its, pdetect = length(x$pas.det.time)/10))
    #plot(fit.gamma)
  }else dists <- rbind(dists,c(shape =c(0), rate =c(10^-5), scenario = its, pdetect = length(x$pas.det.time)/10))
}
names(dists)<-c("shape","rate","scenario","pdetect")
dists$shape <- as.numeric(dists$shape)
dists$rate <- as.numeric(dists$rate)
dists$mean <- dists$shape/dists$rate
dists$var <- dists$shape/dists$rate^2
dists


#plot those that are needed
ggplot(data.frame(x = c(1:20000)*.001), aes(x))+
  lapply(c(1:8), function(k){stat_function(fun=dgamma, args=list(shape=dists$shape[k], rate = dists$rate[k]))})+xlim(c(0,40))

#plot relation distribution parameters and vaccination
ggplot(data = dists)+
  geom_point(aes(x = size, y = mean, colour = as.factor(vac)))+
  geom_point(aes(x = size, y = var, colour = as.factor(vac)))

fit.size <-  lm(mean~size + I(size^2),dists)
drop1(fit.size)
fit.size <-  lm(mean~ + I(size),dists)
drop1(fit.size)

fit.shape <- lm(shape~vac + I(vac^2)+ I(vac^3),dists)
summary(fit.shape)

ggplot(data= dists,aes(x = vac, y = shape))+
  geom_smooth()+
  geom_point()+
  geom_path(aes(y = predict(fit.shape)))

fit.rate <- lm(rate~vac + I(vac^2)+ I(vac^3),dists)
summary(fit.rate)


ggplot(data= dists,aes(x = vac, y = rate))+
  geom_smooth()+
  geom_point()+
  geom_path(aes(y = predict(fit.rate)))

saveRDS(fit.rate, "./output/pasdetRate_clinprot_pasdet.RDS")
saveRDS(fit.shape, "./output/pasdetShape_clinprot_pasdet.RDS")

#human exposure Clinical protection ####
#as we used multiple repeats to determine active surveillance repeat
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
humanexposure.layer.clin <- left_join(humanexposure.layer.clin, baseline.average.clin)

#calculate ratio
humanexposure.layer.clin$ratio.det <- humanexposure.layer.clin$detection.exposure/humanexposure.layer.clin$mean.det.exposure
humanexposure.layer.clin$ratio.tot <- humanexposure.layer.clin$total.exposure/humanexposure.layer.clin$mean.tot.exposure

saveRDS(humanexposure.layer.clin, file = "./output/humanexposureLayerClin.RDS")
humanexposure.layer.clin <- readRDS(file = "./output/humanexposureLayerClin.RDS")


ggplot(humanexposure.layer.clin) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = detection.method),binwidth = .15,alpha = 0.5,colour = "black",position = "identity")+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  lims(x = c(0,3))+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~.,labeller = labeller(scenario = scenario.label.clin,
                                            detection.method = detection.method.label))+
  ggtitle("Human exposure layers")#+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer_clinprot.png")

#only passive
ggplot(humanexposure.layer.clin%>%filter(detection.method =="pas.det.time")) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = detection.method),binwidth = .15,alpha = 0.5,colour = "black",position = "identity")+
  #  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "red",min.det.time= "orange"),
  #                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  lims(x = c(0,3))+
  facet_grid(scenario~.,labeller = labeller(scenario = scenario.label.clin))+
  ggtitle("Human exposure layers")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer_clinprot_pas.png")

#average ratios

humanexposure.layer.clin%>%reframe(.by = c(scenario,detection.method),
                            mean = mean(ratio.det),
                            min = min(ratio.det),
                            max = max(ratio.det),
                            perc25 = quantile(ratio.det,0.25),
                            perc75 =quantile(ratio.det,0.75)
)



fit <- lm(detection.exposure ~  detection.time, 
          data = humanexposure.layer.clin%>%filter(rep == 1 & detection.method=="pas.det.time"& is.finite(detection.time)))
summary(fit)
drop1(fit)
predict(fit, newdata = data.frame(detection.time = 10))

saveRDS(fit, file = "./output/exposurefitLayerVaccinClinprot.RDS")

#Waning immunity ####
#add the 32 000 animal farm without vaccination as baseline
scenario.wane.list <- list(list(scenario = "layerSize32000Vac0"),list(scenario = "layerStartTime0"),list(scenario = "layerStartTime10"),list(scenario = "layerStartTime100"))
output.layer.wane <- lapply(c(1:4),function(i){load.sims(paste0("./output/",gsub(scenario.wane.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)$output})
pars.layer.wane <- lapply(c(1:4),function(i){load.sims(paste0("./output/",gsub(scenario.wane.list[[i]]$scenario,pattern = "[.]", replacement = "")), interval = 0.01)}$pars[[1]])



#visualize ##

scenario.label <- c(layerSize32000Vac0 = "Baseline",
                    layerStartTime0 = "T0 = 0",
                    layerStartTime10 = "T0 = 10",
                    layerStartTime100 = "T0 = 100")

detection.method.label <- c(pas.det.time = "Passive",
                            ac.det.time = "Active",
                            min.det.time = "Minimum")

for(i in c(1:length(scenario.wane.list))){
  show(plot.output(output.layer.wane[[i]],c("I.1","I.2","R.1","R.2"), scenario.wane.list[[i]]$scenario))
  ggsave(paste0("./output/figures/", gsub(scenario.wane.list[[i]]$scenario.wane,pattern = "[.]", replacement = ""),".png"))
  gc()  
}

plot.output.grid(rbind(cbind(output.layer.wane[[1]],scenario = pars.layer.wane[[1]]$scenario),
                       cbind(output.layer.wane[[2]],scenario = pars.layer.wane[[2]]$scenario),
                       cbind(output.layer.wane[[3]],scenario = pars.layer.wane[[3]]$scenario),
                       cbind(output.layer.wane[[4]],scenario = pars.layer.wane[[4]]$scenario)), 
                 vars =c("I.1","I.2","R.1","R.2"),
                 title = "Waning Immunity", scales = "free_y" , scenario.label = scenario.label)
ggsave(file = "./output/figures/waningLayer.png")

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
               WaningRate = max(pars.layer.wane[[i]]$transRate),
               introTime = pars.baseline.layer$max.time-max(pars.layer.wane[[i]]$max.time))
  surveillance.layer.wane <-if(!exists("surveillance.layer.wane")){tmp}else {rbind(surveillance.layer.wane,tmp)}
  
}
surveillance.layer.wane<- cbind(surveillance.layer.wane,surveillance.layer.wane$scenario%>%gsub(pattern ="layer.wane_", replacement = "")%>%str_split_fixed(pattern = c("intro"),2)%>%as.data.frame())


saveRDS(surveillance.layer.wane,"./output/surveillancelayerwane.RDS")
surveillance.layer.wane<-readRDS(surveillance.layer.wane,"./output/surveillancelayerwane.RDS")

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
  ggtitle("Surveillance")+
  facet_grid(scenario~.,labeller = labeller(scenario = scenario.label))


ggsave("./output/figures/scenarios_layer.wane_surveillance.png")

surveillance.layer.wane%>%reframe(.by = scenario, mean.pas = mean(pas.det.time),mean.ac = mean(ac.det.time))

#human exposure Waning immunity####
#as we used multiple repeats to determine active surveillance repeat
humanexposure.layer.wane <- NULL

for(k in c(1:length(scenario.wane.list))){
  humanexposure.layer.wane <- rbind(humanexposure.layer.wane,
                                    cbind(human.exposure.detection.multiple.runs(output.layer.wane[[k]],
                                                                                 pars.layer.wane[[k]]$beta,
                                                                                 surveillance.layer.wane%>%filter(scenario == pars.layer.wane[[k]]$scenario)), 
                                          scenario = pars.layer.wane[[k]]$scenario))
  
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

ggplot(humanexposure.layer.wane) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = detection.method),binwidth = .15,alpha = 0.5,colour = "black",position = "identity")+
  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "orange",min.det.time= "red"),
                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~.,labeller = labeller(scenario = scenario.label))+
  ggtitle("Human exposure layers")#+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer_wane.png")

#only passive
ggplot(humanexposure.layer.wane%>%filter(detection.method =="pas.det.time")) +
  geom_histogram(aes(ratio.det,after_stat(.15*density), fill = detection.method),binwidth = .15,alpha = 0.5,colour = "black",position = "identity")+
#  scale_fill_manual(values=c(pas.det.time = "black",ac.det.time = "red",min.det.time= "orange"),
#                    labels = c(pas.det.time = "Passive",ac.det.time = "Active",min.det.time= "Minimum"))+
  xlab("Risk ratio of exposure \n (reference baseline)") + 
  ylab("Proportion")+
  geom_vline(xintercept = 1)+
  facet_grid(scenario~.,labeller = labeller(scenario = scenario.label))+
  ggtitle("Human exposure layers")+theme(legend.position = 'none')
ggsave("./output/figures/humanexposurelayer_wane_pas.png")


#average ratios
humanexposure.layer.wane%>%reframe(.by = c(scenario,detection.method),
                                   mean = mean(ratio.det),
                                   min = min(ratio.det),
                                   max = max(ratio.det),
                                   perc25 = quantile(ratio.det,0.25),
                                   perc75 =quantile(ratio.det,0.75)
)
