#########################################################
#                                                        
#                  Stochastic multitype transitions model                                 
#                  
#                                                        
#                  Author: Egil Fischer                              
#                  Contact: e.a.j.fischer@uu.nl                             
#                  Creation date                         
#########################################################
#include libraries ####
packages <- c("ggplot2","deSolve","tidyverse","bbmle")


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
  runs =5, #number of runs
  max.time = 19*30,#length of the run
  itypes = itypes, #types####
  deltat = 1, #time step#
  N0 = 75000, #population size####
  initial= 1 , #initially infected
  p.protect = 1, #1 - 6/26,#proportion initially protected by vaccination
  beta =  matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
  latency.period = c(.1,.1),#choose very short latency period
  k.latency =2,##k parameter Erlang distribution
  infectious.period = c(3.0,4.0),#Duration infectious period 
  #variance.infectious.period = c(3.0,4.0)^2, #Variance infectious period
  k.infectious = 20, #k parameter Erlang distribution
  transRate = matrix(c(0,0.012,0.0,0), nrow = itypes), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  trans.mean.wane = 514, #only considers transition from type 2 to 1
  trans.var.wane = 85^2, #only considers transition from type 2 to 1
  trans.mean.buildup = 14, #only considers transition from type 2 to 1
  trans.var.buildup = 7, #only considers transition from type 2 to 1
  pdie = c(1.0,0.01),#probability of dying at end of infectious period
  mortRate = 0,#0.0005 #per capita death rate #Mortality events
  intro.time = 0 #time at which the infection is introduced
)

#gamma distributed time until build-up of immunity - notice that if trans.var == (trans.mean^2) this distribution is an exponential distribution
gamma.buildup.distribution <- function(t, param.list){

  #gamma distribution giving the probability of waning at that time point
 with(param.list,{
    shape_buildup = (trans.mean.buildup^2)/trans.var.buildup;
    scale_buildup =trans.var.buildup/trans.mean.buildup;
    if(is.nan(pgamma(t+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE))) {
      return(0);
    }
    if(pgamma(t+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)-(pgamma(t, shape = shape_buildup, scale = scale_buildup, lower.tail = TRUE)) == 0){return(0)}
    #probability of transition is probability changed at time t + delta - probability at time t given that you did not change yet 
    round((pgamma(t+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)-(pgamma(t, shape = shape_buildup, scale = scale_buildup, lower.tail = TRUE)))/(1.-pgamma(t, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)),digits = 25)
  })
  
}

#initial number of high titre given the time to build-up of immunity
inits.gamma.buildup <- with(param.list,{
  #adjust p.protect with waning
  shape_buildup = (trans.mean.buildup^2)/trans.var.buildup;
  scale_buildup =trans.var.buildup/trans.mean.buildup;
  #fraction of those that will be protected p.protect that is already protected
  p.protect.adjust <- p.protect*(pgamma(intro.time, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE));
  
  return(list(
    L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected 
    I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
    R0 = c(0,0), # number of initially recovered
    S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
  )})


#gamma distributed time until waning - notice that if trans.var == (trans.mean^2) this distribution is an exponential distribution
gamma.waning.distribution <- function(t, param.list){
  #gamma distribution giving the probability of waning at that time point
  with(param.list,{
    shape_wane = (trans.mean.wane^2)/trans.var.wane;
    scale_wane =trans.var.wane/trans.mean.wane;
    round((pgamma(t, shape = shape_wane , scale = scale_wane, lower.tail = FALSE)-(pgamma(t+deltat, shape = shape_wane , scale = scale_wane, lower.tail = FALSE)))/(pgamma(t, shape = shape_wane, scale = scale_wane, lower.tail = FALSE)),digits = 15)
    
  })
}

#initial number of high titre given the time to transition
inits.gamma <- with(param.list,{
  #adjust p.protect with waning
  p.protect.adjust <- p.protect*pgamma(intro.time, shape = (trans.mean.wane^2)/trans.var.wane , scale = trans.var.wane/trans.mean.wane, lower.tail = FALSE);
  
  return(list(
  L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected 
  I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
  R0 = c(0,0), # number of initially recovered
  S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
)})


#initial number of high titre given time of introduction taking into account build-up and waning of immunity
inits.gamma.buildup.wane <-  with(param.list,{
  
  #adjust p.protect with waning
  shape_buildup = (trans.mean.buildup^2)/trans.var.buildup;
  scale_buildup =trans.var.buildup/trans.mean.buildup;
  #fraction of those that will be protected p.protect that is already protected
  p.protect.adjust <- p.protect*(pgamma(intro.time, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE));
  
  #adjust p.protect with waning
  shape_wane = (trans.mean.wane^2)/trans.var.wane;
  scale_wane =trans.var.wane/trans.mean.wane;
  p.protect.adjust <- p.protect.adjust*pgamma(intro.time, shape = shape_wane , scale = scale_wane, lower.tail = FALSE);
  
  return(list(
    L0 = round(c(1-p.protect.adjust,p.protect.adjust)*initial,digits = 0), #number of initially latently infected 
    I0 = round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0), #number of initially infectious
    R0 = c(0,0), # number of initially recovered
    S0 = round(c(1-p.protect.adjust,p.protect.adjust)*(N0-initial),digits = 0))#initially susceptible
  )});inits.gamma.buildup.wane 

#TODO initial values should be within simulator

#function of simulation of spread with an exponentially distributed transition between states
sim.multitypeSEIR_tleap_expwaning <- function(param.list,init, seed = NULL){
  with(c(param.list,init),{
    set.seed(seed)
    #initialize
    currun = 1;
   #calculate the transition rates in the latency and infectious states
    lLat <- k.latency/latency.period
    lInf <- k.infectious/infectious.period
    
    #run the simulations
    while(currun <= runs)
    {
    #initialize the state of the system with k.latency and k.infectious states
    #set state of the system with initial values
    state <- list(time = intro.time, N = N0,C=0, run = currun, dt = 0)
    state$S<- S0
    state$L<- matrix(0, itypes,k.latency)
    state$L[,1] <- L0
    state$I<- matrix(0, itypes,k.infectious)
    state$I[,1] <- I0
    state$R<- R0
    state$DS <- 0*S0
    state$DL <- 0*L0
    state$DI <- 0*I0
    state$DR <- 0*R0
    
   
     
      #one run until maximum time or no infected animals left
      while(state$time < max.time & sum(state$L)+sum(state$I)>0)
      {
        #record this moment
        if(!exists("output")) output <- NULL;
        output$time <- rbind(output$time, state$time);
        output$run <- rbind(output$run, run = currun);
        output$S <- rbind(output$S, state$S);
        output$L <- rbind(output$L, rowSums(state$L));
        output$I <- rbind(output$I, rowSums(state$I));
        #output$Lstates <- rbind(output$Lstates, state$L);
        #output$Istates <- rbind(output$Istates, state$I);
        output$R <- rbind(output$R, state$R);
        output$DS <- rbind(output$DS, c(state$DS));
        output$DL <- rbind(output$DL, c(state$DL));
        output$DI <- rbind(output$DI, c(state$DI));
        output$DR   <- rbind(output$DR, c(state$DR));
        output$N <- rbind(output$N, c(state$N))
         
        #first transitions from one vaccination to the other vaccination state
        #first transitions from high titre to low titre
        transProb <-c(1-exp(-transRate*deltat))
        vacState <- matrix(mapply(rbinom, MoreArgs = list(n = 1), size = c(state$S,state$S),prob = transProb), ncol = itypes)
        #only off diagonals are of interest 
        diag(vacState) <- 0
        diag(vacState) <- state$S - rowSums(vacState) #S remaining in the current state
        #add all transition up 
        state$S <- colSums(vacState)
        #only off diagonals are of interest 
        
        
        #calculate the force of infection
        foi <- (beta %*% rowSums(state$I))/state$N
        
        #determine transitions to next state or death
        transS <- rbinom(rep(1, itypes), state$S,1-exp(-(foi+mortRate)*deltat))
        transL <- matrix(rbinom(rep(1,itypes*k.latency), state$L, 1-exp(-(lLat+mortRate)*deltat)), nrow = itypes)
        transI <- matrix(rbinom(rep(1,itypes * k.infectious), state$I, 1-exp(-(lInf+mortRate)*deltat)), nrow = itypes)
        mortR <- rbinom(rep(1, itypes), state$R, 1-exp(-mortRate*deltat))
        
        #divide between death and transition
        mortS <- round(transS*mortRate/(foi+10^-100 + mortRate)) #prevent dividing by 0
        StoL <- transS - mortS
        mortL <- round(transL*mortRate/(lLat + mortRate))
        LtoNext <- cbind(StoL, transL - mortL)
        mortI <- round(transI*mortRate/(lInf + mortRate)) 
        ItoNext <- cbind(LtoNext[,ncol(LtoNext)], transI - mortI)
        
        #deal with changes
        state$S <- state$S - transS
        state$DS <-state$DS + mortS
        state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
        state$DL <- state$DL + rowSums(mortL)
        state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)] 
        #those dying at recovery
        mortIdisease <- rbinom(rep(1, itypes),ItoNext[,dim(ItoNext)[2]] ,pdie)
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        #update time
        state$time = state$time + deltat
        
       
     
      }
      currun = currun +1
    }
    return(output)
    
  })
  
  
}

# rm(x);x <-sim.multitypeSEIR_tleap_expwaning(param.list,inits, seed = NULL)
# 
# 
# ggplot(data = data.frame(x))+ 
#   geom_path(aes(x = time, y = S.1 + S.2,group = run, colour = "S"))+
#   geom_path(aes(x = time, y = L.1 + L.2,group = run, colour = "L"))+
#   geom_path(aes(x = time, y = I.1 + I.2,group = run, colour = "I"))+
#   geom_path(aes(x = time, y = R.1 + R.2,group = run, colour = "R"))
# 
# ggplot(data = data.frame(x))+ 
#   geom_path(aes(x = time, y = S.1,group = run, colour = "S.1"))+
#   geom_path(aes(x = time, y = S.2,group = run, colour = "S.2"))
# 
# ggplot(data = data.frame(x))+ 
#   geom_path(aes(x = time, y = I.1,group = run, colour = "I.1"))+
#   geom_path(aes(x = time, y = I.2,group = run, colour = "I.2"))

#function of simulation of spread with an arbitrary distributed transition between states
sim.multitypeSEIR_tleap_distwaning <- function(param.list,init,waning.distribution, seed = NULL){
  with(c(param.list,init),{
    set.seed(seed)
    #initialize
     currun = 1;
    #calculate the transition rates in the latency and infectious states
    lLat <- k.latency/latency.period
    lInf <- k.infectious/infectious.period
    
    #run the simulations
    while(currun <= runs)
    {
      #initialize the state of the system with k.latency and k.infectious states
      #set state of the system with initial values
      state <- list(time = intro.time, N = N0,C=0, run = currun, dt = 0)
      state$S<- S0
      state$L<- matrix(0, itypes,k.latency)
      state$L[,1] <- L0
      state$I<- matrix(0, itypes,k.infectious)
      state$I[,1] <- I0
      state$R<- R0
      state$DS <- 0*S0
      state$DL <- 0*L0
      state$DI <- 0*I0
      state$DR <- 0*R0
      
      
      
      #one run until maximum time or no infected animals left
      while(state$time < max.time & (sum(state$L)+sum(state$I)>0 | state$time < intro.time) )
      {
        #record this moment
        if(!exists("output")) output <- NULL;
        output$time <- rbind(output$time, state$time);
        output$run <- rbind(output$run, run = currun);
        output$S <- rbind(output$S, state$S);
        output$L <- rbind(output$L, rowSums(state$L));
        output$I <- rbind(output$I, rowSums(state$I));
        #output$Lstates <- rbind(output$Lstates, state$L);
        #output$Istates <- rbind(output$Istates, state$I);
        output$R <- rbind(output$R, state$R);
        output$DS <- rbind(output$DS, c(state$DS));
        output$DL <- rbind(output$DL, c(state$DL));
        output$DI <- rbind(output$DI, c(state$DI));
        output$DR   <- rbind(output$DR, c(state$DR));
        output$N <- rbind(output$N, c(state$N));
        output$dt <- deltat;
        
        #first transitions from high titre to low titre
        transProb <-c(0,waning.distribution(state$time, param.list),0,0)
        vacState <- matrix(mapply(rbinom, MoreArgs = list(n = 1), size = c(state$S,state$S),prob = transProb), ncol = itypes)
        #only off diagonals are of interest 
        diag(vacState) <- state$S - rowSums(vacState) #S remaining in the current state
        #add all transition up 
        state$S <- colSums(vacState)
        
        
        #calculate the force of infection
        foi <- (beta %*% rowSums(state$I))/state$N
      
        #determine transitions to next state or death
        transS <- rbinom(rep(1, itypes), state$S,1-exp(-(foi+mortRate)*deltat))
        transL <- matrix(rbinom(rep(1,itypes*k.latency), state$L, 1-exp(-(lLat+mortRate)*deltat)), nrow = itypes)
        transI <- matrix(rbinom(rep(1,itypes * k.infectious), state$I, 1-exp(-(lInf+mortRate)*deltat)), nrow = itypes)
        mortR <- rbinom(rep(1, itypes), state$R, 1-exp(-mortRate*deltat))
        
        #divide between death and transition
        mortS <- round(transS*mortRate/(foi+10^-100 + mortRate)) #prevent dividing by 0
        StoL <- transS - mortS
        mortL <- round(transL*mortRate/(lLat + mortRate))
        LtoNext <- cbind(StoL, transL - mortL)
        mortI <- round(transI*mortRate/(lInf + mortRate)) 
        ItoNext <- cbind(LtoNext[,ncol(LtoNext)], transI - mortI)
        
        #deal with changes
        state$S <- state$S - transS
        state$DS <-state$DS + mortS
        state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
        state$DL <- state$DL + rowSums(mortL)
        state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)] 
        #those dying at recovery
        mortIdisease <- rbinom(rep(1, itypes),ItoNext[,dim(ItoNext)[2]] ,pdie)
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        #update time
        state$time = state$time + deltat
        
        
        
      }
      currun = currun +1
    }
    return(output)
    
  })
  
  
}
# 
# rm(x);x <-sim.multitypeSEIR_tleap_distwaning(param.list,
#                                              inits.gamma,
#                                              gamma.waning.distribution, seed = NULL)
# 
# 
# ## #####
# ggplot(data = data.frame(x))+ 
#   geom_path(aes(x = time, y = S.1 + S.2,group = run, colour = "S"))+
#   geom_path(aes(x = time, y = L.1 + L.2,group = run, colour = "L"))+
#   geom_path(aes(x = time, y = I.1 + I.2,group = run, colour = "I"))+
#   geom_path(aes(x = time, y = R.1 + R.2,group = run, colour = "R"))
# 
# with(param.list,{
# ggplot(data = data.frame(x))+ 
#   geom_point(aes(x = x$time, y = N0*pgamma(x$time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean )))+
#   geom_path(aes(x = time, y = S.1,group = run, colour = "S.1"))+
#   geom_path(aes(x = time, y = S.2,group = run, colour = "S.2"))
# })
# 
# ggplot(data = data.frame(x))+ 
#   geom_path(aes(x = time, y = I.1,group = run, colour = "I.1"))+
#   geom_path(aes(x = time, y = I.2,group = run, colour = "I.2"))


#function of simulation of spread with an arbitrary distributed transition between states
sim.multitypeSEIR_tleap_buildup_distwaning <- function(param.list,init,buildup_distribution, waning.distribution, seed = NULL){
  with(c(param.list,init),{
    set.seed(seed)
    #initialize
    currun = 1;
    #calculate the transition rates in the latency and infectious states
    lLat <- k.latency/latency.period
    lInf <- k.infectious/infectious.period
    #output
    output = NULL;
    #run the simulations
    while(currun <= runs)
    {
      #initialize the state of the system with k.latency and k.infectious states
      #set state of the system with initial values
      state <- list(time = intro.time, N = N0,C=0, run = currun, dt = 0)
      state$S<- S0
      state$L<- matrix(0, itypes,k.latency)
      state$L[,1] <- L0
      state$I<- matrix(0, itypes,k.infectious)
      state$I[,1] <- I0
      state$R<- R0
      state$DS <- 0*S0
      state$DL <- 0*L0
      state$DI <- 0*I0
      state$DR <- 0*R0
      
      
      
      #one run until maximum time or no infected animals left
      while(state$time < max.time & (sum(state$L)+sum(state$I)>0 ) )
      {
        #record this moment
        if(!exists("output")) output <- NULL;
        output$time <- rbind(output$time, state$time);
        output$run <- rbind(output$run, run = currun);
        output$S <- rbind(output$S, state$S);
        output$L <- rbind(output$L, rowSums(state$L));
        output$I <- rbind(output$I, rowSums(state$I));
        #output$Lstates <- rbind(output$Lstates, state$L);
        #output$Istates <- rbind(output$Istates, state$I);
        output$R <- rbind(output$R, state$R);
        output$DS <- rbind(output$DS, c(state$DS));
        output$DL <- rbind(output$DL, c(state$DL));
        output$DI <- rbind(output$DI, c(state$DI));
        output$DR   <- rbind(output$DR, c(state$DR));
        output$N <- rbind(output$N, c(state$N));
        output$dt <- deltat;
       
        #Transitions between titre status
        transProb <-c(0,waning.distribution(state$time, param.list),buildup_distribution(state$time,param.list),0)
        vacState <- matrix(mapply(rbinom, MoreArgs = list(n = 1), size = c(state$S,state$S),prob = transProb), ncol = itypes)
        #use expected number of transitions
        #vacState <- matrix(ceiling(c(state$S,state$S)*transProb), ncol = itypes)
        
        #off diagonals are transition form high to low or low to high titre 
        diag(vacState) <- state$S - rowSums(vacState) #S remaining in the current state
        #add all transition up 
        state$S <- colSums(vacState)
        
        #calculate the force of infection
        foi <- (beta %*% rowSums(state$I))/state$N
        
        #determine transitions to next state or death
        transS <- rbinom(rep(1, itypes), state$S,1-exp(-(foi+mortRate)*deltat))
        transL <- matrix(rbinom(rep(1,itypes*k.latency), state$L, 1-exp(-(lLat+mortRate)*deltat)), nrow = itypes)
        transI <- matrix(rbinom(rep(1,itypes * k.infectious), state$I, 1-exp(-(lInf+mortRate)*deltat)), nrow = itypes)
        mortR <- rbinom(rep(1, itypes), state$R, 1-exp(-mortRate*deltat))
        
        #divide between death and transition
        mortS <- round(transS*mortRate/(foi+10^-100 + mortRate)) #prevent dividing by 0
        StoL <- transS - mortS
        mortL <- round(transL*mortRate/(lLat + mortRate))
        LtoNext <- cbind(StoL, transL - mortL)
        mortI <- round(transI*mortRate/(lInf + mortRate)) 
        ItoNext <- cbind(LtoNext[,ncol(LtoNext)], transI - mortI)
        
        #deal with changes
        state$S <- state$S - transS
        state$DS <-state$DS + mortS
        state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
        state$DL <- state$DL + rowSums(mortL)
        state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)] 
        #those dying at recovery
        mortIdisease <- rbinom(rep(1, itypes),ItoNext[,dim(ItoNext)[2]] ,pdie)
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        #update time
        state$time = state$time + deltat
        
        
        
      }
      currun = currun +1
    }
    return(output)
    
  })
  
  
}
# 
# rm(x);x <-sim.multitypeSEIR_tleap_distwaning(param.list,
#                                              inits.gamma,
#                                              gamma.waning.distribution, seed = NULL)
# 
# 
# ## #####
# ggplot(data = data.frame(x))+ 
#   geom_path(aes(x = time, y = S.1 + S.2,group = run, colour = "S"))+
#   geom_path(aes(x = time, y = L.1 + L.2,group = run, colour = "L"))+
#   geom_path(aes(x = time, y = I.1 + I.2,group = run, colour = "I"))+
#   geom_path(aes(x = time, y = R.1 + R.2,group = run, colour = "R"))
# 
# with(param.list,{
# ggplot(data = data.frame(x))+ 
#   geom_point(aes(x = x$time, y = N0*pgamma(x$time, shape = (trans.mean^2)/trans.var , scale = trans.var/trans.mean )))+
#   geom_path(aes(x = time, y = S.1,group = run, colour = "S.1"))+
#   geom_path(aes(x = time, y = S.2,group = run, colour = "S.2"))
# })
# 
# ggplot(data = data.frame(x))+ 
#   geom_path(aes(x = time, y = I.1,group = run, colour = "I.1"))+
#   geom_path(aes(x = time, y = I.2,group = run, colour = "I.2"))
