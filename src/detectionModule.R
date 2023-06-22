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
                                               detectables.vars, #names of columns with detectables animals
                                               se, #sensitivity of test per type of detectables animal,
                                               time.interval.ac = 1,  #time interval for binning of detectables animals
                                               init.ac = 0, #time when active monitoring starts
                                               detectables.incidence = FALSE, #either a boolean or vector of length detectables.vars with booleans indicating whether the values need to be transformed to incidence data. 
                                               pfarm, #probability that this farm is being monitored
                                               panimal, #probability per detectables animal to being monitored,
                                               seed = NULL, #seed
                                               roundTime = 1,# times will be rounded to multiple of roundTime. If roundTime == Null no rounding
                                               ...
                                               
    )
{
  #iterate over runs
  for(i in c(1:max(output$run))){
    
    #iterate over the number of reps
    for(j in c(1:reps))
    {
      
      #arrange data of this run for detection
      det.data <- arrange.data.detection(output%>%filter(run==i), 
                                         #passive detection
                                         deaths.vars,
                                         #active detection
                                         detectables.vars,
                                         time.interval.ac,
                                         ifelse(init.ac =="rand",round(runif(1, min = 0, max = time.interval.ac),roundTime),init.ac),
                                         detectables.incidence,
                                         #overall
                                         roundTime)
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
  return(surveillance.time)
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
  #no detectables animals
  if(max(detectables)==0){return(Inf)}
  
  #based on sensitivity the probability of those that are detectables at a certain moment are used to draw from a binomial distribution
  #return the first time for which the number of successes is larger than one (if none return Inf)
  ac.det.time <- first(c(times[apply(X = detectables,1, FUN = function(x){sum(mapply(d= x,s = panimal*se, FUN = function(d,s){rbinom(size =d, n =1, p = s)}))})>=1],Inf));
  return(ac.det.time);
}

#arrange data in such a way that it can be used for detection functions ####
arrange.data.detection <- function(output, 
                                   death.vars,
                                   detectables.vars,
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
  tmp.active <- output%>%select(all_of(c("time",detectables.vars)))
  tmp.active$time <- init.ac + ceiling((tmp.active$time-init.ac)/time.interval.ac)*time.interval.ac
  #bin the detectables to time.interval.ac times
  tmp.active <- unique(tmp.active%>%group_by(time) %>%mutate_at(vars(!matches(c("time")) ),max))
  #if incidence is required 
  if(detectables.incidence){
    tmp.active <-cbind(time = tail(tmp.active$time,-1), 
                       tail(tmp.active[detectables.vars],-1) - head(tmp.active[detectables.vars],-1))
    if(min(tmp.active[,detectables.vars])<0){
      warning("Negative incidences!")}
  }

  
  return(list(times.passive = tmp.passive$times,
                  deaths = tmp.passive$deaths,
                  times.active = tmp.active$time,
                  detectables = tmp.active[detectables.vars]))
}


