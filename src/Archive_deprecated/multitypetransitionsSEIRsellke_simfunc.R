#########################################################
#                                                        
#                  MULTI-TYPE Sellke SEIR model
#                  
#                                                        
#                  Author:          E.A.J.Fischer                     
#                  Contact:         e.a.j.fischer@uu.nl                     
#                  Creation date:   16-2-2023                         
#########################################################

#Notes
  #dL and dT are still exponential 
  #still using lots of rbind 

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

#define number of types
itypes = 2;

#create a list of all parameters
param.list <- list(
  scenario = "Test", #scenario
  runs = 10, #number of runs
  max.time = 19*30,#length of the run
  itypes = itypes, #types####
  N0 = 46000, #population size####
  initial= 10 , #initially infected
  p.protect = 0.001, #1 - 6/26,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  latency.period = c(1,1),
  k.latency =2,
  infectious.period = c(3.0,4.0),#Duration infectious period
  k.infectious = 20,
  #variance.infectious.period = c(3.0,4.0)^2, #Variance infectious period
  transRate = matrix(c(0,0,0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(1.0,0.01),#probability of dying at end of infectious period
  mortRate = 0.0 #per capita death rate #Mortality events
)

inits <- with(param.list,{list(
  L0 = round(c(1-p.protect,p.protect)*initial,digits = 0), #number of initially latently infected 
  I0 = round(c(1-p.protect,p.protect)*0,digits = 0), #number of initially infectious
  R0 = c(0,0), # number of initially recovered
  S0 = round(c(1-p.protect,p.protect)*(N0-initial),digits = 0)#initially susceptible
)})

functions <- with(param.list, {
  list(
    dL = function(U, itype) {
      shape <- k.latency
      scale <- latency.period[itype] / k.latency
      return(rgamma(n = 1, shape = shape, rate = 1 / scale))
    },
    dI = function(U, itype) {
      shape <- k.infectious  # Shape parameter for Erlang distribution
      scale <- infectious.period[itype] / k.infectious  
      return(rgamma(n = 1, shape = shape, rate = 1 / scale))  
    },
    dT = function(U, itype1, itype2) {
      if (transRate[itype1, itype2] == 0) {
        return(Inf)  # No transition if rate is zero
      } else {
        return(-log(1 - U) / transRate[itype1, itype2])
      }
    },
    dLE = function(U) {
      if (mortRate == 0) {
        return(Inf)  # No mortality if mortality is zero
      } else {
        return(-log(1 - U) / mortRate)
      }
    }
  )
})


#function to run the model
sim.multitypeSEIR <- function(param.list,init, functions, seed = NULL)
{
  with(c(param.list,init,functions),{
    set.seed(seed)
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
    currun = 1;
    
    
    #loop over number of runs
    while(currun <= runs){
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
      
      start.time.run <- Sys.time()
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
      
      #schedule transition events for all animals
      for(i in c(1:itypes)){#type from which to transition
        for(j in c(1:itypes)) #type transition to
        {
          if(i!=j){
            events <- rbind(events,data.frame(type = c(paste0("T",i,j)), time = sapply(FUN = dT, X = runif(state$S[i]), 
                                                                                       itype1 = i, itype2 = j),
                                              id = indiv.states[indiv.states$state=="S" & indiv.states$itype == i,]$id))
          }
        }
      }
      #remove events that will never happen
      events <- filter(events,time != Inf)
      events <- filter(events,time < max.time)
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
      while(length(events$time) > 0 & sum(state$L) + sum(state$I) > 0 & state$time < max.time)
      {
        handbreak = handbreak + 1
        if(handbreak > 5*N0){stop}
        #process the first event in the list
        
        #update the cumulative infection pressure for all animals
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
            
            #count death
            state$DR[indiv.states[first(events$id),]$itype] <- state$DR[indiv.states[first(events$id),]$itype]  +1 
            
            #set individual state
            indiv.states[first(events$id),]$state<- "DR";
            
            #remove events of this animal but keep the current event
            events<- rbind(first(events),subset(events, id != first(events$id)))
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
          if(exists("maxInfperiod")){ maxInfperiod<- max(maxInfperiod, infPeriod,na.rm =T) }else maxInfperiod <- infPeriod;
          indiv.states[first(events$id),]$dieatrec <- (rbinom(n =1, size = 1, p = pdie[indiv.states[first(events$id),]$itype]))>0
          
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
            events<- rbind(first(events),subset(events, id != animalid))
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
        
        #type transition event
        if(grepl("T",first(events$type)))
        {
          #if uninfected transfer the individual to the other type
          if(indiv.states[indiv.states$id == first(events$id),]$state == "S"){
            itype1 <- as.numeric(substring(first(events$type),2,2))
            itype2 <- as.numeric(substring(first(events$type),3,3))
            
            indiv.states[indiv.states$id == first(events$id),]$itype = itype2;
            state$S[itype1]= state$S[itype1] -1 
            state$S[itype2]= state$S[itype2] +1 
            
          }
          
          
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
      #record output and parameters
      op <- list(out = output, pars = param.list)
      saveRDS(op, file = paste0(format(Sys.Date(),"%Y%m%d"),"output",scenario,"",currun,".RDS"))
      
      #add on to the current run counter
      currun <- currun +1
      
    }
    return(output);
  })}; timing.overall <-Sys.time() - start.time.overal;



output <- sim.multitypeSEIR(param.list, inits, functions, seed = 1440)


setwd("C:/Users/plett/OneDrive/Desktop/HPAI_vaccination")
files <- list.files(pattern = "outputTest.*.RDS")
all_outputs <- lapply(files, readRDS)

combined_output <- bind_rows(
  lapply(seq_along(all_outputs), function(run) {
    # Extract the `out` component and convert to a data frame
    out <- as.data.frame(all_outputs[[run]]$out)
    out$run <- run  # Add a column for the run number
    return(out)
  })
)
# PLOTTING
#Dynamics
 ggplot(data = combined_output)+
   geom_path(aes(x = time, y = S.1 + S.2, group = run, colour = "S")) +
   #geom_path(aes(x = time, y = S.1, group = run, colour = "S.1")) +
   #geom_path(aes(x = time, y = S.2, group = run, colour = "S.2")) +
   #geom_path(aes(x = time, y = L.1 + L.2,group = run, colour = "L"))+
   #geom_path(aes(x = time, y = I.2,group = run, colour = "I"))+
   #geom_path(aes(x = time, y = R.2,group = run, colour = "R"))+
   #geom_path(aes(x = time, y = DR.1 + DR.2, group = run, colour = "DR")) +
   #geom_path(aes(x = time, y = DI.1 + DI.2, group = run, colour = "DI")) +
   #geom_path(aes(x = time, y = DL.1 + DL.2, group = run, colour = "DL")) +
   labs(
     title = "SEIR Chickens",
     subtitle = paste("Population Size (N0):", param.list$N0, "| Introduction Time:", param.list$intro.time),
     x = "Time",
     y = "Number of Chickens",
     colour = "Compartment"
  )


#Final size and number of small outbreaks
 fs <- combined_output %>%
   group_by(run) %>%
   summarize(fs = max(R.1 + R.2 + DI.1 + DI.2 + DR.1 + DR.2 + DL.1 + DL.2, na.rm = TRUE))
 
# Define small outbreaks as outbreaks that are less than 10% of the population size
 small_outbreaks <- sum(fs$fs < 0.1 * param.list$N0)
 total_runs <- param.list$runs
 proportion_small_outbreaks <- small_outbreaks / total_runs
 print(proportion_small_outbreaks)

 ggplot(fs, aes(x = fs)) +
   geom_histogram(binwidth = 2000, fill = "salmon", color = "black", boundary = 0) +
   labs(
     title = "Distribution of Outbreak Sizes",
     subtitle = paste(
       "Number of small outbreaks (<10% of population):", small_outbreaks, "out of", total_runs, "runs\n",
       "Population Size (N0):", param.list$N0
     ),
     x = "Outbreak Size",
     y = "Frequency"
   ) +
   theme_minimal(base_size = 14) +
   theme(
     plot.title = element_text(hjust = 0.5, face = "bold"),
     plot.subtitle = element_text(hjust = 0.5),
     axis.title = element_text(face = "bold")
   )
