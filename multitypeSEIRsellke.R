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
packages <- c("ggplot2","dplyr","tidyverse")


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
itypes <- 2;

#these should be fitting to the number of types
N0 <- 20000
#proportion initially protected by vaccination
p.protect = 0.5;
initial.latent <- 10 
L0 <- round(c(1-p.protect,p.protect)*initial.latent,digits = 0) #number of initially latently infected 
I0 <- round(c(1-p.protect,p.protect)*0,digits = 0) #number of initially infectious
R0 <- c(0,0) # number of initially recoverd
S0 <- round(c(1-p.protect,p.protect)*N0,digits = 0)-L0-I0-R0
runs <- 5 #number of runs
N0 <- sum(S0+L0+I0+R0)

#parameters
#Type 1  = not protected by vaccination and type 2 = protected by vaccination
#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976
beta <- matrix(c(3.73, 3.73,0.058,0.058),ncol = itypes) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)
#choose very short latency period
latency.period <- c(0.001,0.001);
#Duration latency period as function of a random variable U
dL <- function(U,itype){return(-log(1 - U)*latency.period[itype])} #exponential distribution

#Duration infectious period as function of a random variable U
infectious.period <- c(1.47,1.47)
variance.infectious.period <- c(1.47,1.47)^2
#Duration infectious period as function of a random variable U
#dI <- function(U){return(-log(1 - U) * infectious.period)} #exponential function
dI <- function(U,itype){return(qgamma(U, shape = infectious.period[itype]^2/variance.infectious.period[itype] , rate = infectious.period[itype]/variance.infectious.period[itype]))} #exponential function

#probability of dying at end of infectious period
pdie <- c(0,0.7)

#
mortRate <- 0.0005/N0
dLE <- function(U){return(-log(1 -U)/mortRate)} #life expectancy as function of random variable U



#initialize
output <- list(time = c(0), N = c(N0),C=c(0), run =c(1),dt = c(0))
output$S<- S0
output$L<- L0
output$I<- I0
output$R <- R0
output$DS <- 0*S0
output$DL <- 0*L0
output$DI <- 0*I0
output$DR <- 0*R0


#seed <- as.numeric(Sys.time())
set.seed(1678980340)#(seed)
currun = 1;
#loop over number of runs
while(currun <= runs)
{
  #set state of the system with initial values
  state <- list(time = 0, N = N0,C=0, run = currun, dt = 0)
  state$S<- S0
  state$L<- L0
  state$I<- I0
  state$R <- R0
  state$DS <- 0*S0
  state$DL <- 0*L0
  state$DI <- 0*I0
  state$DR <- 0*R0
  
  ids.initI =  if(sum(state$I) > 0){1:sum(state$I)} else 0;
  ids.initL = (sum(state$I)+1) : (sum(state$I)+sum(state$L));
  ids.initS = c((sum(state$I)+sum(state$L)+1):(sum(state$S)+sum(state$I)+sum(state$L)))
  
  #individual states
  indiv.states <- rbind(if(sum(state$I)>0){data.frame(id = ids.initI,
                             state = "I",
                             itype = c(unlist(mapply(rep, times = state$I, x= c(1:length(state$I))))))}else data.frame(id =c(),state =c(),itype=c()),
                        data.frame(id = ids.initL,
                                         state = "L",
                                   itype = c(unlist(mapply(rep, times = state$L, x= c(1:length(state$L)))))),
                         data.frame(id = ids.initS,
                                    state = "S",
                                   itype = c(unlist(mapply(rep, times = state$S, x= c(1:length(state$S)))))))
  #die at recovery
  indiv.states$dieatrec <- FALSE;
  
  #debugging stuff
  indiv.states$state.or <- indiv.states$state
  indiv.states$itype.or <- indiv.states$itype
  
                                    
                        
  
  #initialize event list with id so that dead animals can be removed
  events <- if(sum(state$I)>0){data.frame(type = c("IR"), 
                       time =  mapply(dI,runif(sum(state$I)), indiv.states[indiv.states$state == "I",]$itype ),
                       id = ids.initI)} else {data.frame(type = c(), 
                                                         time =  c(),
                                                         id = c())};
  
  events <- rbind(events,data.frame(type = c("LI","IR"), 
                                    time = c(mapply(function(x,y){cumsum(c(x,y))},  c(mapply(dL,runif(sum(state$L)), indiv.states[indiv.states$state == "L",]$itype )),
                                                    c(mapply(dI,runif(sum(state$L)), indiv.states[indiv.states$state == "L",]$itype )))),
                                    id = rep(ids.initL, each = 2)))
  
 
  #schedule mortality events for all animals 
  events <- rbind(events,data.frame(type = c("M"), time = sapply(FUN = dLE, X = runif(state$N)),id = indiv.states$id))
  #order events by time of execution
  events <- events[order(events$time),]
 
 
  
  #for susceptible animal draw a threshold Q and latent and infectious period and set cumulative foi to 0
  QIRtimes<- data.frame(Q = rexp(sum(state$S),1),
                        M = sapply(FUN = dLE, X = runif(sum(state$S))),
                        cumfoi = rep(0,sum(state$S)),
                        id = ids.initS,
                        infectiontime = Inf);
  QIRtimes<- QIRtimes[order(QIRtimes$Q),]
  
  
  
  #set initial values 
  foi <- rep(0, itypes)  #force of infection
  
  #counters for output/debugging not necessary for running the simulations
  cLI <-  rep(0, itypes)  #counter for latent to infectious transitions
  cIR <- rep(0, itypes)   #counter for infectious to recover transitions
  cSL <- rep(0, itypes)  #counter for susceptible to latent transitions
  cM <- rep(0, itypes)  #counter for mortality events
  
  #set manual 'handbreak' that prevents endless loops
  handbreak =0
  while(length(events$time) > 0 & sum(state$L) + sum(state$I) > 0)
  {
    handbreak = handbreak + 1
    if(handbreak > 5*N0){stop}
    #process the first event in the list
    
    #update the cummulative infection pressure for all animals
    QIRtimes$cumfoi <- QIRtimes$cumfoi +foi[indiv.states[QIRtimes$id,]$itype] * (first(events$time) - state$time)
    
    #set time to current time and time step
    state$dt <- (first(events$time) - state$time)
    state$time <- first(events$time)
    
    #determine the next event
    #latent to infectious
    if(first(events$type) == "LI" ){
      cLI <- cLI + 1
      #add one to I type 
      state$I[indiv.states[first(events$id),]$itype] <- state$I[indiv.states[first(events$id),]$itype]  + 1
      #subtract one from L type
      state$L[indiv.states[first(events$id),]$itype]  <- state$L[indiv.states[first(events$id),]$itype]  - 1
      #set individual state
      indiv.states[first(events$id),]$state<- "I"
      if(min(state$L)<0){
        print("stop")
      }
      
    }
    #infectious to recovered
    if(first(events$type) == "IR" )
    {
      cIR <- cIR + 1
      #add one to R
      state$R[indiv.states[first(events$id),]$itype] <- state$R[indiv.states[first(events$id),]$itype] + 1
      #subtract one from I
      if(state$I[indiv.states[first(events$id),]$itype] ==0){
        print("error state$I" )
      }
      #
      state$I[indiv.states[first(events$id),]$itype]  <- state$I[indiv.states[first(events$id),]$itype]  - 1
      #set individual state
      indiv.states[first(events$id),]$state<- "R";
      #if the individual dies now deal with it
      if(indiv.states[first(events$id),]$dieatrec)
      {
        #subtract from R
        state$R[indiv.states[first(events$id),]$itype] <- state$R[indiv.states[first(events$id),]$itype] - 1
        state$N <- state$N - 1
        
        #set individual state
        indiv.states[first(events$id),]$state<- "DR";
        
        #remove events
        events<- subset(events, id != first(events$id))
      }
      
    }
    #susceptible to infectious
    if(first(events$type) == "SL" )
    { 
      #get individual type
      itype = indiv.states[first(events$id),]$itype
      cSL[itype]  <- cSL[itype] + 1
      #add one to L
      state$L[itype]  <- state$L[itype]  + 1
      #subtract one from S
      state$S[itype]  <- state$S[itype]  - 1
      #set transitions
      latPeriod <- dL(runif(1),itype);
      infPeriod <- dI(runif(1),itype);
      events<- rbind(events, data.frame(time =c(state$time + latPeriod), 
                                        type = c("LI"), 
                                        id = first(events$id))) #L -> I
      events<- rbind(events, data.frame(time =c(state$time + latPeriod+infPeriod),
                                        type = c("IR"), 
                                        id = first(events$id))) #I -> R
      indiv.states[first(events$id),]$dieatrec <- (rbinom(n =1, size = 1, p = pdie))>0
      
      #remove resistance of this particular animal
      QIRtimes <-subset(QIRtimes,QIRtimes$id !=first(events$id))
      #set individual state
      indiv.states[first(events$id),]$state<- "L"
      
    }
    #death
    if(first(events$type)== "M")
    {
      cM <- cM + 1
      state$N <- state$N - 1
    
      #get the animal id
      animalid  = first(events$id)
      typeM = indiv.states[indiv.states$id == animalid,]$state
      itype = indiv.states[indiv.states$id == animalid,]$itype
      
      #remove threshold 
      if(typeM == "S") { 
        #remove resistance of this anomal
        QIRtimes <-QIRtimes[QIRtimes$id != animalid,]
        } #or events
      else if(typeM =="L" || typeM == "I"){
        #select and remove events of this animal
        events<- subset(events, id != animalid)
      }
      
     
      if(state[typeM][[1]][itype] <= 0){
        print("error mortality event");
        stop();
      }
      
      
      #deal with the removal
      indiv.states[indiv.states$id == animalid,]$state = paste0("D",typeM); #set individual state to be dead in certain state
      state[typeM][[1]][itype] = state[typeM][[1]][itype] - 1; #remove one from the counters
      state[paste0("D",typeM)][[1]][itype] = state[paste0("D",typeM)][[1]][itype] + 1; #add one to the death counters
    }
    
    if(min(state$L)<0){
      print("stop")
    }
    #order events by time of execution
    events <- events[order(events$time),]
    
    #set force-of-infection for each type
    foi <- c(beta %*% state$I /state$N)
    
    #order events by time of execution
    events <- events[order(events$time),]
    
    
    #determine next infection event
    if(sum(state$S) * sum(state$I) > 0){
      QIRtimes$infectiontime <- state$time + (QIRtimes$Q - QIRtimes$cumfoi)/foi[indiv.states$itype[QIRtimes$id]]
    } else {infection <- 10^10}
    
    #remove first event
    events <- events[-1,]
    
    #if this event is previous to other events schedule it
    if(length(events$time) > 0)
    {
      if(min(QIRtimes$infectiontime)  < first(events$time))
      {
        #get the first infection
        QIRtimes <-QIRtimes[order(QIRtimes$infectiontime),]
        #
        events <- rbind(data.frame(time = first(QIRtimes$infectiontime),type = "SL",
                                   id = first(QIRtimes$id) ),events)
        #order events by time of execution
        events <- events[order(events$time),]
      }
    }
    
    #record this moment
    output$time <- rbind(output$time, state$time);
    output$N <- rbind(output$N, state$N);
    output$C <- rbind(output$C, state$C);
    output$run <- rbind(output$run, state$run);
    output$S <- rbind(output$S, state$S);
    output$L <- rbind(output$L, state$L);
    output$I <- rbind(output$I, state$I);
    output$R <- rbind(output$R, state$R);
    output$DS <- rbind(output$DS, state$DS);
    output$DL <- rbind(output$DL, state$DL);
    output$DI <- rbind(output$DI, state$DI);
    output$DR   <- rbind(output$DR, state$DR);
    output$dt   <- rbind(output$dt, state$dt);
    
  }
  currun <- currun +1
}
##plot the output
ggplot(data = 
         data.frame(output))+
  #geom_path(aes(x = time, y = S.1, group = run), color = "magenta", linetype = 1)+
  #geom_path(aes(x = time, y = S.2, group = run), color = "magenta", linetype = 2)+
  # geom_step(aes(x = time, y = L.1, group = run), color = "blue", linetype = 1)+
  # geom_step(aes(x = time, y = L.2, group = run), color = "blue", linetype = 2)+
  geom_step(aes(x = time, y = I.1, group = run), color = "red")+
  geom_step(aes(x = time, y = I.2, group = run), color = "red", linetype = 2)+
  #geom_step(aes(x = time, y = R.1, group = run), color = "darkgreen")+
  # geom_step(aes(x = time, y = R.2, group = run), color = "darkgreen", linetype = 2)+
  # geom_step(aes(x = time, y = DS.1 + DL.1 +DI.1+DR.1, group = run), color = "black")+
  # geom_step(aes(x = time, y = DS.2 + DL.2 +DI.2+DR.2, group = run), color = "black", linetype = 2)+
  ylab("L,I,R,D")+facet_grid(run~.)

#post process to exposure
#Assumption of exposure is that the rate at which humans get exposed to infection is proportionate to the infectivity towards chickens
a <-1; #factor scaling towards human exposure
beta.human <- a*beta

library(tidyverse)
output.df<- data.frame(output)
exposure.human <- output.df%>%group_by(run)%>%summarise(time = time,
                                      I.1 = I.1,
                                      I.2= I.2,psurvdt = exp(-beta.human[1,1]*I.1*dt-beta.human[1,2]* I.2*dt))
exposure.human <-exposure.human%>%group_by(run) %>% mutate(pinf = 1 - cumprod(psurvdt))
ggplot(data = exposure.human) + geom_path(aes(x = time, y = psurvdt, group = run))
ggplot(data = exposure.human) + geom_path(aes(x = time, y = pinf, group = run))
