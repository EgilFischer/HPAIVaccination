#########################################################
#                                                        
#                 Post-proces simulations                                
#                                                        
#                  Author: Egil Fischer                               
#                  Contact: e.a.j.fischer@uu.nl                             
#                  Creation date: 28-3-2023                         
#########################################################
# 
# 
# #include libraries ####
# packages <- c("ggplot2","dplyr","tidyverse")
# 
# 
# ## Now load or install&load all
# package.check <- lapply(
#   packages,
#   FUN = function(x) {
#     if (!require(x, character.only = TRUE)) {
#       install.packages(x, dependencies = TRUE)
#       library(x, character.only = TRUE)
#     }
#   }
# )
# 
# #import simulations
# for(i in c(1:10)){
#   
#   if(i == 1)output <- as.data.frame(readRDS(paste0("20230407outputlayerT1001.RDS"))$out)
#   output <- rbind(output,as.data.frame(readRDS(paste0("20230407outputlayerT100",i,".RDS"))$out))
# }
# 
# ##plot the output
plot.output <- function(output,vars,title = NULL ){
  ggplot(data =
           data.frame(output)%>%select(time, run,vars)%>%reshape2::melt(id.vars = c("time","run"),value.name = "prevalence",variable.name=c("itype")))+
        geom_step(aes(x = time, y = prevalence,colour = itype))+
      ylab("#number of birds")+facet_grid(.~itype, scales = "free")+ggtitle(title)
  
}
plot.output(output.baseline.broiler,c("S.1","S.2","I.1","I.2","R.1","R.2"), "Broiler base line")
plot.output(output.baseline.broiler,c("DS.1","DS.2","DI.1","DI.2","DR.1","DR.2"), "Broiler base line")

# 
# #human exposure ####
# #Assumption of exposure is that the rate at which humans get exposed to infection is proportionate to the infectivity towards chickens
# a <-0.001; #factor scaling towards human exposure
# beta.human <- a*beta
# 
# 
# output.df<- data.frame(output)
# exposure.human <- output.df%>%group_by(run)%>%reframe(time = time,
#                                                       I.1 = I.1,
#                                                       I.2= I.2,
#                                                       psurvdt = exp(-beta.human[1,1]*I.1*dt-beta.human[1,2]* I.2*dt))
# exposure.human <- exposure.human%>%group_by(run)%>%mutate(
#   pinf = 1 - cumprod(psurvdt))
# ggplot(data = exposure.human) + 
#   geom_path(aes(x = time, y = psurvdt, colour = factor(run)))
# ggplot(data = exposure.human) + 
#   geom_path(aes(x = time, y = pinf, colour = factor(run)))

#detection time ####
#detect if a certain number of R dead are found within a certain time interval

detection.time.threshold.interval <- function(times,Deaths, time.interval,threshold){
  #determine the number of deaths at the end of the interval
  Dinterval <- data.frame(t(sapply(FUN = function(t){c(time = t,deaths = max(Deaths[times < t]))}, X = seq(1, ceiling(max(times)),time.interval))));
  #define the increment during the interval
  Dinterval$deaths <- c(0, tail(Dinterval$deaths,-1) -  head(Dinterval$deaths,-1))
  #return the first time the threshold is passed /  return Inf if never
  if(max(Dinterval$deaths)<threshold){return(Inf)}
  return(min(Dinterval$time[Dinterval$deaths>=threshold], max(Dinterval$time)));
  
  
}

detection.time.threshold.subsequent <- function(times,Deaths, time.interval, threshold,n){
  #determine the number of deaths at the end of the interval
  Dinterval <- data.frame(t(sapply(FUN = function(t){c(time = t,deaths = max(Deaths[times < t]))}, X = seq(1, ceiling(max(times)),time.interval))));
  #define the increment during the interval
  Dinterval$deaths <- c(0, tail(Dinterval$deaths,-1) -  head(Dinterval$deaths,-1))
  if(max(Dinterval$deaths)<threshold){return(Inf)}
  #find index of first two consecutive times in which deaths exceed the threshold
  #determine whether it is above the threshold
  Dinterval$exceed <- as.numeric(Dinterval$deaths > threshold);
  index <- which(diff(Dinterval$exceed == 1, 2) == -2)
  return(min(Dinterval$time[index], max(Dinterval$time)));
  
}


detection.times <- function(output,vars, detection.function,time.interval, threshold, ...)
{
  detection.times <- data.frame(run =c(), det.time =c())
  for(i in c(1:max(output$run)))
  {
    tmp <- output%>%filter(run == i)%>%select("time", vars)%>%rowwise()%>%reframe(time = time,
                                                                      prevalence = rowSums(.[vars]))
    if(exists("detection.times")){
      detection.times <- rbind(detection.times, data.frame(run = i, det.time = detection.function(tmp$times, tmp$prevalence,time.interval,threshold,...)))}
      else{
      detection.times <- data.frame(run = i, det.time = detection.function(tmp$times, tmp$prevalence,time.interval,threshold,...))
    }
    
  }
  return(detection.times)
}

detection.times(output.baseline.broiler,c("DR.1","DR.2") ,detection.time.threshold.subsequent,1,0.005*75, n= 2)

det.times.interval <- data.frame(scenario =c(), run =c(), det.time =c())
det.times.subseq <- data.frame(scenario =c(), run =c(), det.time =c())
# add all death together
output$D = output$DS.1 + output$DL.1+ output$DI.1 + output$DR.1+output$DS.2 + output$DL.2+ output$DI.2 + output$DR.2
#
for(i in c(1:4)){
  for(j in c(1:10))  {
    det.time<-detection.time.threshold.interval(times = output$time[output$scenario == i & output$run ==j],D = output$D[output$scenario == i & output$run ==j], time.interval = 2,0.005*c(75000,45000,45000,45000)[i])
    det.times.interval <- rbind(det.times,data.frame(scenario = i, run = j, det.time = det.time))
    
    det.time<-detection.time.threshold.subsequent(times = output$time[output$scenario == i & output$run ==j],D = output$D[output$scenario == i & output$run ==j], time.interval = 1,n =2 , 0.005*c(75000,45000,45000,45000)[i])
    det.times.subseq <- rbind(det.times.subseq,data.frame(scenario = i, run = j, det.time = det.time))
  }
}


# 
# det.time <- detection.time.threshold.interval(times = output$time[output$run ==4],D = output$DR.1[output$run ==4]+output$DR.2[output$run ==4], time.interval = 1,0.05*75000)
# det.time
# 
# # #detect if a certain number of R animals are found 
# # detection.time.threshold <- function(times,cumD,threshold){
# #   return(head(times[cumD >= threshold],1))
# # }
# # 
# # detection.time.threshold(output$time[output$run ==1],output$DR.1[output$run ==1]+output$DR.2[output$run ==1], 0.05*20000)
# 
# output.df<- data.frame(output)%>%filter(time <= det.time)
# exposure.human <- output.df%>%group_by(run)%>%summarise(time = time,
#                                                         I.1 = I.1,
#                                                         I.2= I.2,psurvdt = exp(-beta.human[1,1]*I.1*dt-beta.human[1,2]* I.2*dt))
# exposure.human <-exposure.human%>%group_by(run) %>% mutate(pinf = 1 - cumprod(psurvdt))
# ggplot(data = exposure.human) + geom_path(aes(x = time, y = psurvdt, group = run))
# ggplot(data = exposure.human) + geom_path(aes(x = time, y = pinf, group = run))
