#########################################################
#                                                        
#                  MULTI-TYPE Sellke SIR model
#                  
#                                                        
#                  Author:          E.A.J.Fischer                     
#                  Contact:         e.a.j.fischer@uu.nl                     
#                  Creation date:   16-2-2023                         
#########################################################

#load libraries
source("./src/loadLibraries.R") 


#define number of types
itypes = 2;

#create a list of all parameters
param.list <- list(
  scenario = "TestSIR", #scenario
  runs =2, #number of runs
  max.time = 17*30,#length of the run
  itypes = itypes, #types####
  N0 = 450, #population size####
  initial= 10 , #initially infected
  p.hightitre = 1 - 6/26,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2, #Variance infectious period
  transRate = matrix(c(0,0.012,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.01,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005/7 #per capita death rate #Mortality events - based on performance reports and pers. com. mieke matthijs 0.5% per week
)

inits <- with(param.list,{list(
  I0 = round(c(1-p.hightitre,p.hightitre)*initial,digits = 0), #number of initially infectious
  R0 = c(0,0), # number of initially recovered
  S0 = round(c(1-p.hightitre,p.hightitre)*(N0-initial),digits = 0)#initially susceptible
)})

functions <- with(param.list, {list(
 dI = function(U,itype){return(qgamma(U, shape = infectious.period[itype]^2/variance.infectious.period[itype] , rate = infectious.period[itype]/variance.infectious.period[itype]))}, #exponential function#Duration infectious period as function of a random variable U
 dT = function(U,itype1,itype2){return(-log(1 -U)/transRate[itype1,itype2])}, #exponential function
 dLE = function(U){return(-log(1 -U)/mortRate)} #life expectancy as function of random variable U
)})



#function to run the model ####
sim.multitypeSIR <- function(param.list,init, functions, seed = NULL)
{
  with(c(param.list,init,functions),{
  set.seed(seed); 
  #initialize
  output <- list(time = c(0), N = c(N0),C=c(0), run =c(1),dt = c(0))
  output$S<- S0
  output$I<- I0
  output$R <- R0
  output$DS <- 0*S0
  output$DI <- 0*I0
  output$DR <- 0*R0
  currun = 1;
  
  #create output folder for this scenario if it does not yet exists
  path <- paste0("./output/",gsub(pattern= "[.]",replacement = "", scenario));
  if(!file.exists(path)){dir.create(path)}
  
  #loop over number of runs
  while(currun <= runs){
    #initialize
    output <- list(time = c(0), N = c(N0),C=c(0), run =c(currun),dt = c(0))
    output$S<- S0
    output$I<- I0
    output$R <- R0
    output$DS <- 0*S0
    output$DI <- 0*I0
    output$DR <- 0*R0
    
    start.time.run <- Sys.time()
    #set state of the system with initial values
    state <- list(time = 0, N = N0,C=0, run = currun, dt = 0)
    state$S<- S0
    state$I<- I0
    state$R <- R0
    state$DS <- 0*S0
    state$DI <- 0*I0
    state$DR <- 0*R0
    
    ids.initI =  if(sum(state$I) > 0){1:sum(state$I)} else 0;
    ids.initS = c((sum(state$I)+1):(sum(state$S)+sum(state$I)))
    
    #individual states
    indiv.states <- rbind(if(sum(state$I)>0){data.frame(id = ids.initI,
                                                        state = "I",
                                                        itype = c(unlist(mapply(rep, times = state$I, x= c(1:length(state$I))))))}else data.frame(id =c(),state =c(),itype=c()),
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


    #schedule mortality events for all animals 
    events <- rbind(events,data.frame(type = c("M"), time = sapply(FUN = dLE, X = runif(state$N)),id = indiv.states$id))
    
    #schedule transition events for all animals
    for(i in c(1:itypes)){#type from which to transition
      for(j in c(1:itypes)) #type transition to
      {
        if(i!=j){
          if(length(indiv.states[indiv.states$state=="S" & indiv.states$itype == i,]$id)>0)
          {
          events <- rbind(events,data.frame(type = c(paste0("T",i,j)), time = sapply(FUN = dT, X = runif(state$S[i]), 
                                                                                     itype1 = i, itype2 = j),
                                            id = indiv.states[indiv.states$state=="S" & indiv.states$itype == i,]$id))
          }
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
    cIR <- rep(0, itypes)   #counter for infectious to recover transitions
    cSI <- rep(0, itypes)  #counter for susceptible to latent transitions
    cM <- rep(0, itypes)  #counter for mortality events
    
    #set manual 'handbreak' that prevents endless loops
    handbreak =0; 
    while(length(events$time) > 0 &  sum(state$I) > 0 & state$time < max.time)
    { 
      handbreak = handbreak + 1
      if(handbreak > 5*N0){stop("Too many events processed > 5 N0 ")}
      #process the first event in the list
      
      #update the cumulative infection pressure for all animals
      QIRtimes$cumfoi <- QIRtimes$cumfoi +foi[indiv.states[QIRtimes$id,]$itype] * (first(events$time) - state$time)
      
      #set time to current time and time step
      state$dt <- (first(events$time) - state$time)
      state$time <- first(events$time)
      
      #determine the next event
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
      if(first(events$type) == "SI" )
      { 
        #get individual type
        itype = indiv.states[first(events$id),]$itype
        cSI[itype]  <- cSI[itype] + 1
        #add one to L
        state$I[itype]  <- state$I[itype]  + 1
        #subtract one from S
        state$S[itype]  <- state$S[itype]  - 1
        #set transitions
        infPeriod <- dI(runif(1),itype);
        events<- rbind(events, data.frame(time =c(state$time+infPeriod),
                                          type = c("IR"), 
                                          id = first(events$id))) #I -> R
        if(exists("maxInfperiod")){ maxInfperiod<- max(maxInfperiod, infPeriod,na.rm =T) }else maxInfperiod <- infPeriod;
        indiv.states[first(events$id),]$dieatrec <- (rbinom(n =1, size = 1, p = pdie[indiv.states[first(events$id),]$itype]))>0
        
        #remove resistance of this particular animal
        QIRtimes <-subset(QIRtimes,QIRtimes$id !=first(events$id))
        #set individual state
        indiv.states[first(events$id),]$state<- "I"
        
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
        else if(typeM == "I"){
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
          events <- rbind(data.frame(time = first(QIRtimes$infectiontime),type = "SI",
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
      output$I <- rbind(output$I, state$I);
      output$R <- rbind(output$R, state$R);
      output$DS <- rbind(output$DS, state$DS);
      output$DI <- rbind(output$DI, state$DI);
      output$DR   <- rbind(output$DR, state$DR);
      output$dt   <- rbind(output$dt, state$dt);
    }
    #record output and parameters
    op <- list(out = output, pars = param.list)
    
    saveRDS(op, file = paste0(path,"/",format(Sys.Date(),"%Y%m%d"),gsub(pattern = "[.]",replacement = "",x = scenario),"",currun,".RDS"))
    #if(exists("all.output")){all.output <- rbind(all.output,data.frame(output))}else{all.output <- data.frame(output)}
    #add on to the current run counter
    currun <- currun +1
    
  }
  #return(all.output);
})};

# wrapper function so that only requires param.list ####
simulate.multitypeSIR <- function(param.list, seed = NULL){
  #functions to produce initial numbers and functions
  inits <- with(param.list,{list(
    I0 = round(c(1-p.hightitre,p.hightitre)*initial,digits = 0), #number of initially infectious
    R0 = c(0,0) # number of initially recovered
    )});
  inits$S0 <- with(param.list,{round(c(1-p.hightitre,p.hightitre)*N0,digits = 0)})-inits$I0 - inits$R0
  funcs <- with(param.list, {list(
    dI = function(U,itype){return(qgamma(U, shape = infectious.period[itype]^2/variance.infectious.period[itype] , rate = infectious.period[itype]/variance.infectious.period[itype]))}, #exponential function#Duration infectious period as function of a random variable U
    dT = function(U,itype1,itype2){return(-log(1 -U)/transRate[itype1,itype2])}, #exponential function
    dLE = function(U){return(-log(1 -U)/mortRate)} #life expectancy as function of random variable U
  )}) ;
  
  return(sim.multitypeSIR(param.list , init = inits, functions = funcs))
}
