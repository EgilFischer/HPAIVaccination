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

#load simulations
load.sims <- function(path, interval = NULL){
  sims <- list.files(path, pattern = ".RDS");if(length(sims) == 0) stop("No simulations in this folder");
  
  if(is.null(interval)){
  output <- as.data.frame(readRDS(paste0(path,"/",sims[1]))$out)
   for(i in sims){
     if(!exists("output"))
        {output <-as.data.frame(readRDS(paste0(path,"/",i))$out)}else
        {output <- rbind(output,as.data.frame(readRDS(paste0(path,"/",i))$out))}
        }}else #merge into intervals
        {
          tmp <- as.data.frame(readRDS(paste0(path,"/",sims[1]))$out);
          tmp$interval.index <- tmp$time%/%interval;
          tmp$tround <- interval*(tmp$interval.index+sign(tmp$time))
          output <- tmp%>%group_by(tround)%>%slice(n())%>%ungroup
          for(i in sims){
            tmp <- as.data.frame(readRDS(paste0(path,"/",i))$out);
            tmp$interval.index <- tmp$time%/%interval;
            tmp$tround <- interval*(tmp$interval.index+sign(tmp$time))
            tmp <- tmp%>%group_by(tround)%>%slice(n())%>%ungroup
            if(!exists("output")) {output <-tmp}else{ 
            output <- rbind(output,tmp)}
          }
        }
  return(output)
}


# ##plot the output
plot.output <- function(output,vars,title = NULL ){
  ggplot(data =
           data.frame(output)%>%select(time, run,vars)%>%reshape2::melt(id.vars = c("time","run"),value.name = "prevalence",variable.name=c("itype")))+
        geom_step(aes(x = time, y = prevalence,colour = itype, group =run))+
      ylab("#number of birds")+facet_grid(.~itype, scales = "free")+ggtitle(title)
  
}

plot.output.sparse <- function(output,vars,title = NULL, frac = 0.5){
  out <- data.frame(output)%>%sample_n(round(frac*length(output$time)));
  return(plot.output(out, vars, title))
}

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

detection.time.threshold.subsequent <- function(times,Deaths, time.interval, threshold,ints ){
  #determine the number of deaths at the end of the interval
  Dinterval <- data.frame(t(sapply(FUN = function(t){c(time = t,deaths = max(Deaths[times < t]))}, X = seq(1, ceiling(max(times)),time.interval))));
  #define the increment during the interval
  Dinterval$deaths <- c(0, tail(Dinterval$deaths,-1) -  head(Dinterval$deaths,-1))
  if(max(Dinterval$deaths)<threshold){return(Inf)}
  #find index of first n consecutive times in which deaths exceed the threshold
  #determine whether it is above the threshold
  Dinterval$exceed <- as.numeric(Dinterval$deaths > threshold);
  index <- which(diff(cumsum(Dinterval$exceed),ints)==ints)[1]+ ints;
  return(min(Dinterval$time[index], max(Dinterval$time)));
  
}

#
detection.times <- function(output,vars, detection.function,time.interval, threshold, ...){

  
  
  output <-data.frame(output%>%select("time","run", vars)%>%reframe(times = time,
                                                            run = run,
                                                            prevalence = rowSums(.[vars])));
  for(i in c(1:max(output$run)))
  {
    if(exists("detection.times.output")){
      detection.times.output <- rbind(detection.times.output, data.frame(run = i, det.time = detection.function(output[output$run ==i,]$times, output[output$run ==i,]$prevalence,time.interval,threshold,...)))}
      else{
      detection.times.output <- data.frame(run = i, det.time = detection.function(output[output$run ==i,]$times, output[output$run ==i,]$prevalence,time.interval,threshold,...))
    }
    
  }
  
  return(detection.times.output)
}

#calculate the detection time in a surveillance system where detection also occurs when a threshold on dead animals is reached
detection.time.surveillance <- function(times,  #vector with times
                                        Deaths,  #vector with dead animals
                                        time.interval,  #time interval for passive detection
                                        ints,  #number of days with subsequent increased mortality
                                        threshold,  #threshold for passive detection
                                        Detectables, #matrix with size length(times) and groups of detectable animals
                                        pfarm, #probability that this farm is being monitored
                                        panimal, #probability per detectable animal to being monitored,
                                        se, #sensitivity of test per type of detectable animal,
                                        seed = NULL,
                                        roundTime = TRUE,
                                        ...
){
  set.seed(seed);
  #if this farm is not monitored return passive detection 
  pas.det.time <- detection.time.threshold.subsequent(times,Deaths, time.interval, threshold,ints);
  if(pfarm < runif(1)){return(data.frame(pas.det.time,ac.det.time = Inf, min.det.time = min(pas.det.time,Inf), ac.succes = FALSE))};
  
  #determine detection time by active surveillance
  ac.det.time <- first(times[rowSums(apply(Dets, c(1,2), function(x){rbinom(x,n =1, p = se*panimal)}))>1]);
  if(roundTime)ac.det.time <- round(ac.det.time/time.interval)*time.interval;
  #return
  tmp.det <- data.frame(pas.det.time,ac.det.time, min.det.time = min(pas.det.time,ac.det.time), ac.succes = pas.det.time>ac.det.time);
  return(tmp.det)
}

#detection.time.surveillance 
detection.times.surveillance <- function(output,  #output
                                      Deaths.vars,  #vector with dead animals
                                      time.interval,  #time interval for passive detection
                                      ints,  #number of days with subsequent increased mortality
                                      threshold,  #threshold for passive detection
                                      Detectables.vars, #matrix with size length(times) and groups of detectable animals
                                      pfarm, #probability that this farm is being monitored
                                      panimal, #probability per detectable animal to being monitored,
                                      se, #sensitivity of test per type of detectable animal,
                                      seed = NULL,
                                      roundTime = TRUE,
                                      Detectables.incidence = FALSE, #either a boolean or vector of length Detectables.vars with booleans indicating whether the values need to be transformed to incidence data. 
                                      ...
){
  
  tmp.output <-output%>%select("time","run", Deaths.vars)%>%reframe(run = run,
                                                                    times =time,
                                                                    deaths = rowSums(.[Deaths.vars]));
  tmp.output <-cbind(tmp.output, output%>%select(Detectables.vars))
  
  
  for(i in c(1:max(output$run)))
  {
    #temporary values of this run
    tmp <- tmp.output%>%filter(run ==i)
    #if transform to incidence is required 
   if(Detectables.incidence ==TRUE || length(Detectables.incidence)>1)
   {
    #allowing just to say true to all values
    if(length(Detectables.incidence)==1) {Detectables.incidence = rep(TRUE,length(Detectables.vars))}
    #transform those that require transform
    for(i in c(1:length(Detectables.incidence)))
    {
      #transform this particular value
      if(Detectables.incidence[i]){
        #transform output
        
        tmp[,Detectables.vars[i]] <-unlist(c(0, unlist(tail(tmp[,Detectables.vars[i]],-1) -  head(tmp[,Detectables.vars[i]],-1))))
        
      }
    }
  }
  
    
    
    if(exists("detection.times.output")){
      detection.times.output <- rbind(detection.times.output,data.frame(run = i, detection.time.surveillance(tmp$times, tmp$deaths,time.interval,ints,threshold,as.matrix(tmp[,Detectables.vars]),pfarm,panimal, se, seed, roundTime, ...)))}
    else{
      detection.times.output <- data.frame(run = i, detection.time.surveillance(tmp$times, tmp$deaths,time.interval,ints,threshold,as.matrix(tmp[,Detectables.vars]),pfarm,panimal, se, seed, roundTime, ...))
    }
    
  }
  return(detection.times.output)
}

# det.times.interval <- data.frame(scenario =c(), run =c(), det.time =c())
# det.times.subseq <- data.frame(scenario =c(), run =c(), det.time =c())
# # add all death together
# output$D = output$DS.1 + output$DL.1+ output$DI.1 + output$DR.1+output$DS.2 + output$DL.2+ output$DI.2 + output$DR.2
# #
# for(i in c(1:4)){
#   for(j in c(1:10))  {
#     det.time<-detection.time.threshold.interval(times = output$time[output$scenario == i & output$run ==j],D = output$D[output$scenario == i & output$run ==j], time.interval = 2,0.005*c(75000,45000,45000,45000)[i])
#     det.times.interval <- rbind(det.times,data.frame(scenario = i, run = j, det.time = det.time))
#     
#     det.time<-detection.time.threshold.subsequent(times = output$time[output$scenario == i & output$run ==j],D = output$D[output$scenario == i & output$run ==j], time.interval = 1,n =2 , 0.005*c(75000,45000,45000,45000)[i])
#     det.times.subseq <- rbind(det.times.subseq,data.frame(scenario = i, run = j, det.time = det.time))
#   }
# }


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

human.exposure.timeseries.multiple.runs <- function(output, beta.human)
  {
  exposure<- data.frame(output)%>% group_by(run)%>%reframe(time = time,
                                            I.1 = I.1,
                                            I.2= I.2,
                                            cum.exposure = cumsum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  return(exposure)
}


human.exposure.timeseries <- function(output, beta.human){
  exposure<- data.frame(output)%>% reframe(time = time,
                                               I.1 = I.1,
                                               I.2= I.2,
                                               cum.exposure = cumsum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  
}

human.exposure.total.multiple.runs<- function(output, beta.human,detection.time)
{
  exposure<- data.frame(output)%>% group_by(run)%>%reframe(total.exposure = sum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  exp.det<-c();
  for(i in exposure$run){
    exp.det<-rbind(exp.det, data.frame(output)%>%filter(run == i & time<= detection.time[i])%>%reframe(detection.exposure = sum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt)))
  }
  exposure<- cbind(exposure, data.frame(detection.exposure = exp.det))
  return(exposure)
}
# ggplot(data = exposure.human) + geom_path(aes(x = time, y = psurvdt, group = run))
# ggplot(data = exposure.human) + geom_path(aes(x = time, y = pinf, group = run))


