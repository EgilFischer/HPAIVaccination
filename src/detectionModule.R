#########################################################
#                                                        
#                 Detection module                            
#                                                        
#                  Author: Egil Fischer                               
#                  Contact: e.a.j.fischer@uu.nl                             
#                  Creation date: 28-3-2023                         
#########################################################

source("./src/loadLibraries.R") 

#detection time ####

#function to produce passive and active detection times for every run. 
#Determining detection times is repeated 'reps' times for each run
repeat.detection.time.surveillance <- function(output, #output consisting of one or more runs of the model
                                               reps, #number of times to repeat detection
                                               #passive detection input
                                               deaths.vars,  #names of columns with dead animals
                                               time.interval.pas,  #time interval for passive detection
                                               threshold,  #threshold for passive detection
                                               ints,  #number of days with subsequent increased mortality 
                                               #active detection input
                                               detectables.vars, #names of columns with detectable animals
                                               se, #sensitivity of test per type of detectable animal,
                                               time.interval.ac = 1,  #time interval for binning of detectable animals
                                               init.ac =0, #time when active monitoring starts
                                               detectables.incidence = FALSE, #either a boolean or vector of length Detectables.vars with booleans indicating whether the values need to be transformed to incidence data. 
                                               pfarm, #probability that this farm is being monitored
                                               panimal, #probability per detectable animal to being monitored,
                                               seed = NULL, #seed
                                               roundTime = 1,# times will be rounded to multiple of roundTime. If roundTime == Null no rounding
                                               ...
                                               
    )
{
  #iterate over runs
  for(i in c(1:max(output$run))){
    #arrange data of this run for detection
    det.data <- arrange.data.detection(output%>%filter(run==i), 
                                       #passive detection
                                       deaths.vars,
                                       #active detection
                                       detectable.vars,
                                       time.interval.ac,
                                       init.ac,
                                       detectables.incidence,
                                       #overall
                                       roundTime)
    #iterate over the number of reps
    for(j in c(1:reps))
    {
    #determine passive detection time
    pas.det.time <- passive.detection.time.threshold.subsequent(det.data$times.passive,
                                                                det.data$deaths,
                                                                time.interval.pas,
                                                                threshold,
                                                                ints)
    #determine active detection time
    ac.det.time <- active.detection.time(det.data$times.active,
                                         det.data$detectables,
                                         se,
                                         pfarm,
                                         panimal)
    #record
    if(!exists("surveillance.time")){surveillance.time<- data.frame(run =i,
                                                                   rep = j,
                                                                   pas.det.time,
                                                                   ac.det.time,
                                                                   min.det.time = min(pas.det.time,ac.det.time),
                                                                   ac.succes = pas.det.time>ac.det.time)} else{
                                    surveillance.time<-rbind(surveillance.time, data.frame(run =i,
                                                                                           rep = j,
                                                                                           pas.det.time,
                                                                                           ac.det.time,
                                                                                           min.det.time = min(pas.det.time,ac.det.time),
                                                                                           ac.succes = pas.det.time>ac.det.time))                          
                                                                   }
    }
    
  }
}





#Passive detection: detect if a certain number of R dead are found within a certain time interval
passive.detection.time.threshold.interval <- function(times,Deaths, time.interval.pas,threshold){
  #determine the number of deaths at the end of the interval
  Dinterval <- data.frame(t(sapply(FUN = function(t){c(time = t,deaths = max(Deaths[times < t]))}, X = seq(1, ceiling(max(times)),time.interval.pas))));
  #define the increment during the interval
  Dinterval$deaths <- c(0, tail(Dinterval$deaths,-1) -  head(Dinterval$deaths,-1))
  #return the first time the threshold is passed /  return Inf if never
  if(max(Dinterval$deaths)<threshold){return(Inf)}
  return(min(Dinterval$time[Dinterval$deaths>=threshold], max(Dinterval$time)));
  
  
}
#Passive detection: detect if a certain number of R dead are found within subsequent time intervals
passive.detection.time.threshold.subsequent <- function(times,Deaths, time.interval.pas, threshold,ints ){
  #determine the number of deaths at the end of the interval
  Dinterval <- data.frame(t(sapply(FUN = function(t){c(time = t,deaths = max(Deaths[times < t]))}, X = seq(1, ceiling(max(times)),time.interval.pas))));
  #define the increment during the interval
  Dinterval$deaths <- c(0, tail(Dinterval$deaths,-1) -  head(Dinterval$deaths,-1))
  if(max(Dinterval$deaths)<threshold){return(Inf)}
  #find index of first n consecutive times in which deaths exceed the threshold
  #determine whether it is above the threshold
  Dinterval$exceed <- as.numeric(Dinterval$deaths > threshold);
  index <- which(diff(cumsum(Dinterval$exceed),ints)==ints)[1]+ ints;
  return(min(Dinterval$time[index], max(Dinterval$time)));
  
}

#active detection time
active.detection.time <- function(times,detectables, se,pfarm = 1, panimal =1){
  if(length(se)!= length(detectables[1,])&length(se)>1)stop("Not all sensitivities are defined");
  #farm escapes detection by active sampling
  if(runif(1)>pfarm){return(Inf)}
  #based on sensitivity the probability of those that are detectable at a certain moment are used to draw from a binomial distribution
  #return the first time for which the number of successes is larger than one
  ac.det.time <- first(times[apply(X = detectables,1, FUN = function(x){sum(mapply(d= x,s = panimal*se, FUN = function(d,s){rbinom(size =d, n =1, p = s)}))})>1]);
  return(ac.det.time);
}

#arrange data in such a way that it can be used for detection functions ####
arrange.data.detection <- function(output, 
                                   death.vars,
                                   detectable.vars,
                                   time.interval.ac,
                                   init.ac,
                                   detectables.incidence = TRUE,
                                   roundTime =1){
  
  #arrange data for passive detection 
  tmp.passive  <- output%>%select(all_of(c("time", death.vars)))%>%reframe(times =time,
                                                                deaths = rowSums(.[death.vars]));
  if(!is.null(roundTime)){
    tmp.passive$times <- ceiling(tmp.passive$times/roundTime)*roundTime;
    tmp.passive <- tmp.passive%>%reframe(.by = times,
                                         deaths = max(deaths))
  }
  
  #arrange data for active detection
  tmp.active <- output%>%select(all_of(c("time",detectable.vars)))
  tmp.active$time <- ceiling(tmp.active$time/time.interval.ac)*time.interval.ac
  #bin the detectables to time.interval.ac times
  tmp.active <- unique(tmp.active%>%group_by(time) %>%mutate_at(vars(!matches(c("time")) ),max))
  #if incidence is required 
  if(detectables.incidence){
    tmp.active <-cbind(time = tail(tmp.active$time,-1), 
                       tail(tmp.active[detectable.vars],-1) - head(tmp.active[detectable.vars],-1))
    if(min(tmp.active[,detectable.vars])<0){warning("Negative incidences!")}
  }

  
  return(list(times.passive = tmp.passive$times,
                  deaths = tmp.passive$deaths,
                  times.active = tmp.active$time,
                  detectables = tmp.active[detectable.vars]))
}




#calculate the detection time in a surveillance system where detection also occurs when a threshold on dead animals is reached
detection.time.surveillance <- function(times.passive,  #vector with times for passive detection
                                        Deaths,  #vector with dead animals
                                        time.interval.pas,  #time interval for passive detection
                                        ints,  #number of days with subsequent increased mortality
                                        threshold,  #threshold for passive detection
                                        times.active,  #vector with times for active detection
                                        Detectables, #matrix with size length(times) and groups of detectable animals
                                        pfarm, #probability that this farm is being monitored
                                        panimal, #probability per detectable animal to being monitored,
                                        se, #sensitivity of test per type of detectable animal,
                                        seed = NULL,
                                        time.interval.ac = 1,
                                        roundTime = TRUE,
                                        ...
){
  set.seed(seed);
  #passive detection time
  pas.det.time <- pas.detection.time.threshold.subsequent(times,Deaths, time.interval.pas, threshold,ints);
  if(pfarm < runif(1)){
    #if this farm is not monitored return passive detection 
    return(data.frame(pas.det.time,ac.det.time = Inf, min.det.time = min(pas.det.time,Inf), ac.succes = FALSE))};
  
  #active detection time
  for(i in c(1:max(output$run)))
  {
    #temporary values of this run
    tmp <- tmp.output%>%filter(run ==i)
    #bin in time intervals. First round time to bin and than select the value at the latest value in the bin
    tmp$times <- round(tmp$times/time.interval.ac)*time.interval.ac
    #bin the detectables to time.interval.ac times
    tmp <- tibble(times = times,as_tibble(Detectables))
    tmp <- unique(tmp%>%group_by(times) %>%mutate_at(vars(!matches(c("times")) ),max))
    
    
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
    stop()
    
    #determine detection time by active surveillance
    if(length(se)!= length(Detectables[1,])&length(se)>1)stop("Not all sensitivities are defined");
    ac.det.time <- first(times[apply(X = Detectables,1, FUN = function(x){sum(mapply(d= x,s = se, FUN = function(d,s){rbinom(size =d, n =1, p = s)}))})>1]);
    if(roundTime)ac.det.time <- round(ac.det.time/time.interval.ac)*time.interval.ac;
    #replace NA by Inf as this one in undetected
    ac.det.time[ is.na(ac.det.time)] <- Inf
    #return
    tmp.det <- data.frame(pas.det.time,ac.det.time, min.det.time = min(pas.det.time,ac.det.time), ac.succes = pas.det.time>ac.det.time);
    return(tmp.det)
  }
  
#detection.time.surveillance both active and passive repetition
detection.times.surveillance.one.rep <- function(output,  #output
                                                   Deaths.vars,  #vector with dead animals
                                                   time.interval.pas,  #time interval for passive detection
                                                   ints,  #number of days with subsequent increased mortality
                                                   threshold,  #threshold for passive detection
                                                   Detectables.vars, #matrix with size length(times) and groups of detectable animals
                                                   pfarm, #probability that this farm is being monitored
                                                   panimal, #probability per detectable animal to being monitored,
                                                   se, #sensitivity of test per type of detectable animal,
                                                   seed = NULL,
                                                   roundTime = TRUE,
                                                   time.interval.ac = 1, 
                                                   init.ac =0,
                                                   Detectables.incidence = FALSE, #either a boolean or vector of length Detectables.vars with booleans indicating whether the values need to be transformed to incidence data. 
                                                   ...
  ){
    
    tmp.output <-output%>%select("time","run", Deaths.vars)%>%reframe(run = run,
                                                                      times =time,
                                                                      deaths = rowSums(.[Deaths.vars]));
    tmp.output <-cbind(tmp.output, output%>%select(Detectables.vars))
    
    
    
    
    
    if(exists("detection.times.output")){
      detection.times.output <- rbind(detection.times.output,
                                      data.frame(run = i, 
                                                 detection.time.surveillance(tmp$times, tmp$deaths,time.interval.pas,ints,threshold,as.matrix(tmp[,Detectables.vars]),pfarm,panimal, se, seed,time.interval.ac, roundTime, ...)))}
    else{
      detection.times.output <- data.frame(run = i, 
                                           detection.time.surveillance(tmp$times, tmp$deaths,time.interval.pas,ints,threshold,as.matrix(tmp[,Detectables.vars]),pfarm,panimal, se, seed,time.interval.ac, roundTime, ...))
    }
    
  }
  return(detection.times.output)
}


#detection.time.surveillance both active and passive 1 simulation multiple detection times
repeat.detection.time.surveillance <- function(output,  #output
                                               Deaths.vars,  #vector with dead animals
                                               time.interval.pas,  #time interval for passive detection
                                               ints,  #number of days with subsequent increased mortality
                                               threshold,  #threshold for passive detection
                                               Detectables.vars, #matrix with size length(times) and groups of detectable animals
                                               pfarm, #probability that this farm is being monitored
                                               panimal, #probability per detectable animal to being monitored,
                                               se, #sensitivity of test per type of detectable animal,
                                               seed = NULL,
                                               roundTime = TRUE, #rounding times
                                               time.interval.ac = 1,  #interval for binning active surveillance
                                               Detectables.incidence = FALSE, #either a boolean or vector of length Detectables.vars with booleans indicating whether the values need to be transformed to incidence data. 
                                               init.time = 0, #start active surveillance
                                               reps = 1, #default one repetiton
                                               ...
){
  for(i in c(1:reps)){
    detection.tmp <-detection.times.surveillance.one.rep(output,  #output
                                                         Deaths.vars,  #vector with dead animals
                                                         time.interval.pas,  #time interval for passive detection
                                                         ints,  #number of days with subsequent increased mortality
                                                         threshold,  #threshold for passive detection
                                                         Detectables.vars, #matrix with size length(times) and groups of detectable animals
                                                         pfarm, #probability that this farm is being monitored
                                                         panimal, #probability per detectable animal to being monitored,
                                                         se, #sensitivity of test per type of detectable animal,
                                                         seed,
                                                         roundTime,
                                                         time.interval.ac,
                                                         init.ac,
                                                         Detectables.incidence,
                                                         ...);
    if(!exists("sout.tmp")){sout.tmp <- cbind(detection.tmp,data.frame(rep = i))}
    else{sout.tmp<- rbind(sout.tmp,cbind(detection.tmp,data.frame(rep = i)))};
  }
  return(sout.tmp)
}

# det.times.interval <- data.frame(scenario =c(), run =c(), det.time =c())
# det.times.subseq <- data.frame(scenario =c(), run =c(), det.time =c())
# # add all death together
# output$D = output$DS.1 + output$DL.1+ output$DI.1 + output$DR.1+output$DS.2 + output$DL.2+ output$DI.2 + output$DR.2
# #
# for(i in c(1:4)){
#   for(j in c(1:10))  {
#     det.time<-detection.time.threshold.interval(times = output$time[output$scenario == i & output$run ==j],D = output$D[output$scenario == i & output$run ==j], time.interval.pas = 2,0.005*c(75000,45000,45000,45000)[i])
#     det.times.interval <- rbind(det.times,data.frame(scenario = i, run = j, det.time = det.time))
#     
#     det.time<-detection.time.threshold.subsequent(times = output$time[output$scenario == i & output$run ==j],D = output$D[output$scenario == i & output$run ==j], time.interval.pas = 1,n =2 , 0.005*c(75000,45000,45000,45000)[i])
#     det.times.subseq <- rbind(det.times.subseq,data.frame(scenario = i, run = j, det.time = det.time))
#   }
# }


# 
# det.time <- detection.time.threshold.interval(times = output$time[output$run ==4],D = output$DR.1[output$run ==4]+output$DR.2[output$run ==4], time.interval.pas = 1,0.05*75000)
# det.time
# 
# # #detect if a certain number of R animals are found 
# # detection.time.threshold <- function(times,cumD,threshold){
# #   return(head(times[cumD >= threshold],1))
# # }
# # 
# # detection.time.threshold(output$time[output$run ==1],output$DR.1[output$run ==1]+output$DR.2[output$run ==1], 0.05*20000)
# 