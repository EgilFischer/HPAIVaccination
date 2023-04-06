#########################################################
#                                                        
#                  MULTI-TYPE Sellke SEIR model
#                  Include detection events
#                                                        
#                  Author:          E.A.J.Fischer                     
#                  Contact:         e.a.j.fischer@uu.nl                     
#                  Creation date:   16-2-2023                         
#########################################################

#include libraries ####
packages <- c("ggplot2","dplyr")


## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#parameters####
animals <- 20000 #population size
beta <- 0.2 #transmission coefficient
dL <- function(U){return(-log(U)/0.5)} #Duration latency period as function of a random variable U
dI <- function(U){return(-log(U)/0.1)} #Duration infectious period as function of a random variable U
dLE <- function(U){return(-log(U)/0.0001)} #life expectancy as function of random variable U
L0 <- 5 #number of initially latently infected 
I0 <- 3 #number of initially infectious
R0 <- 0 # number of initially recoverd
runs <- 5 #number of runs



#initialize
output <- data.frame(time = 0, N = 0, S = 0 , L = 0, I = 0, R= 0, DS = 0,DL = 0, DI= 0, DR= 0, C =0, Th = 0, run = 0)
state <- output

seed <- as.numeric(Sys.time()); set.seed(seed)

#loop over number of runs
while(state$run < runs)
{
  #set state of the system with initial values
  state <- data.frame(time = 0, N = animals, S = animals - L0 , L = L0, I = I0, R= R0,DS = 0,DL =0, DI=0, DR=0, C=0, Th =0, run = last(state$run) + 1)
  #get for every animal the infection threshold Q, the length of the latent period and length of infectious period
  QIRtimes<- data.frame(Q = sort(rexp(state$S,1)),
                        E = sapply(FUN = dL, X = runif(state$S)), 
                        I = sapply(FUN = dI, X = runif(state$S)))
  #get for every animal the moment of death
  Mtimes <- data.frame(M = sapply(FUN = dLE, X = runif(state$N)));
  
  #initialize event list with id so that dead animals can be removed
  ids =  1:state$I;
  events <- data.frame(type = c("IR"), time = c(dI(runif(state$I))),id = ids);
  ids = (state$I+1) : (state$I+state$L)
  events <- rbind(events,data.frame(type = c("LI","IR"), time = c(replicate(cumsum(c(dL(runif(1)),dI(runif(1)))),n=state$L)),id = ids))
  #schedule mortality events
  events <- rbind(events,data.frame(type = c("M"), time = sapply(FUN = dLE, X = runif(state$N)),id = c(1:state$N)))
  #order events by time of execution
  events <- events[order(events$time),]
  ids = max(ids)+1;
  
  #set initial values 
  time <- 0 #time
  foi <- 0  #force of infection
  cumInf <- 0 #cumulative force of infection
  
  #counters for output/debugging not necessary for running the simulations
  cLI <- 0 #counter for latent to infectious transitions
  cIR <-0  #counter for infectious to recover transitions
  cSL <-0 #counter for susceptible to latent transitions
  cM <- 0 #counter for mortality events
  
  #set manual 'handbreak' that prevents endless loops
  handbreak =0
  while(length(events$time) > 0 & state$L + state$I > 0)
  {
    handbreak = handbreak + 1
    if(handbreak > 10*animals){stop}
    #process the first event in the list
    #update the cummulative infection pressure
    cumInf <- cumInf + foi * (first(events$time) - state$time)
    #set time to current time
    state$time <- first(events$time)
    #determine the next event
    #latent to infectious
    if(first(events$type) == "LI" ){
      cLI <- cLI + 1
      #add one to I type 
      state$I <- state$I + 1
      #subtract one from L type
      state$L <- state$L - 1
    }
    #infectious to recovered
    if(first(events$type) == "IR" )
    {
      cIR <- cIR + 1
      #add one to R
      state$R <- state$R + 1
      #subtract one from I
      if(state$I==0){
        print("error")
      }
      state$I <- state$I - 1
      
    }
    #susceptible to infectious
    if(first(events$type) == "SL" )
    {
      cSL <- cSL + 1
      #add one to L
      state$L <- state$L + 1
      #subtract one from S
      state$S <- state$S - 1
      #set transitions
      events<- rbind(events, data.frame(time =c(state$time + QIRtimes[1,]$E), type = c("LI"), id = ids)) #L -> I
      events<- rbind(events, data.frame(time =c(state$time + QIRtimes[1,]$E + QIRtimes[1,]$I), type = c("IR"), id = ids)) #I -> R
      ids = ids + 1
      
      #remove lowest resistance
      QIRtimes <-QIRtimes[-1,] 
      
    }
    #death
    if(first(events$type)== "M")
    {
      cM <- cM + 1
      
      if(min(state[c("S","L","I","R")])<0)
      {
        print(state);
      }
      
      #determine type
      typeM<- sample(c("S","L","I","R"), 1, prob = c(state$S,
                                                     state$L,
                                                     state$I,
                                                     state$R))
     
      if(state[typeM]==0)
      {
        print('error');
      }
      
      #remove threshold or events
      if(typeM == "S") { 
        #remove random resistance
        QIRtimes <-QIRtimes[-sample.int(length(QIRtimes$Q),1),] 
        #count
        state$DS =  state$DS +1}
      
      if(typeM == "L"){
        #remove LI and IR event
        animalid <- sample(subset(events, type == "LI")$id,1)
        #select and remove events of this animal
        events<- subset(events, id != animalid)
        #count
        state$DL =  state$DL +1
      }
      
      if(typeM == "I"){
        #remove  IR events
        animalid <- sample(events[ave(paste(events$id),events$id,FUN=function(x) length(x))==1,]$id,1)
        #select and remove events of this animal
        events<- subset(events, id != animalid)
        #count
        state$DI =  state$DI +1
      }
      if(typeM == "R")
      {
        #count
        state$DR =  state$DR +1
      }
     
      if(state[typeM] ==0){
        print("error")
      }
      #deal with the removal
      state[typeM] = state[typeM] - 1
      
    }
    
    
    
    #set force-of-infection
    foi <- beta * state$I /state$N
    
    #order events by time of execution
    events <- events[order(events$time),]
    
    
    #determine next infection event
    if(state$S * state$I > 0){
      infection <- state$time + (first(QIRtimes$Q) - cumInf)/foi
    } else {infection <- 10^10}
    
    #remove first event
    events <- events[-1,]
    
    #if this event is previous to other events schedule it
    if(length(events$time) > 0)
    {
      if(infection < first(events$time))
      {
        events <- rbind(data.frame(time = c(infection),type = "SL",id = NA),events)
        #order events by time of execution
        events <- events[order(events$time),]
      }
    }
    
    
    #record this moment
    output <- rbind(output, state)
  }
  output <- tail(output,-1)
}


##plot the output
ggplot(data = 
         output)+
 # geom_path(aes(x = time, y = S, group = run), color = "magenta")+
  geom_path(aes(x = time, y = L, group = run), color = "blue")+
  geom_path(aes(x = time, y = I, group = run), color = "red")+
  geom_path(aes(x = time, y = R, group = run), color = "green")+
#  geom_path(aes(x = time, y = DS + DL +DI+DR, group = run), color = "black")+
  ylab("L,I,R,D")+facet_grid(run~.)

