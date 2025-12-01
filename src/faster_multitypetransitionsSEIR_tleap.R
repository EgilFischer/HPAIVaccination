#############################################################################################################
#                                                        
#                  Stochastic multitype transitions model with egg production                                 
#                  
#                                                        
#                  Author: Egil Fischer   & Sarah Pletts                           
#                  Contact: e.a.j.fischer@uu.nl                             
#                                          
#############################################################################################################
# 
# #include libraries ####
# packages <- c("ggplot2","deSolve","tidyverse","bbmle")
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
#include libraries ####
# packages <- c("ggplot2","deSolve","tidyverse","bbmle")
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




#function of simulation of spread with transition between states based on gamma distributed time until transition ####
sim.multitypeSEIR_tleap <- function(param.list, seed = NULL){
  if(is.finite(param.list$trans.mean.wane) & is.finite(param.list$trans.mean.buildup)){
    #Both build up and waning
    return(sim.multitypeSEIR_tleap_buildupandwaning(param.list, seed))}
  if(is.finite(param.list$trans.mean.buildup) & !is.finite(param.list$trans.mean.wane)){
    #if build-up is finite, use the sim.multitypeSEIR_tleap_onlybuilup function
    return(sim.multitypeSEIR_tleap_onlybuildup(param.list, seed))}
  if(!is.finite(param.list$trans.mean.buildup) & is.finite(param.list$trans.mean.wane)){
    #if waning is finite, use the sim.multitypeSEIR_tleap_onlywaning function
    return(sim.multitypeSEIR_tleap_onlywaning(param.list, seed))}
  if(!is.finite(param.list$trans.mean.buildup) & !is.finite(param.list$trans.mean.wane)){
    #if both waning and build-up are infinite, use the sim.multitypeSEIR_tleap_noimmunity function
    return(sim.multitypeSEIR_tleap_constantimmunity(param.list, seed))
  }
  
} 


#function to determine the distribution of the number of introduction of titre groups ####
StoLintro_fun <- function(S, N, introduction, random){
  #if length introduction >1 this means it is preset
  if(length(introduction)> 1){
    #error if not all itypes are included
    if(length(introduction)!= length(S)){stop("no.introduction should either be a single number or a vector with length of itypes")}
    #take the minimum of S and number of introduction per group
    return(pmin(S,introduction)) 
  }else if(random)
  {
    # draw from a multinomial distbrition with S/N probabilities. Also take the minimum of S and number of introduction per group
    return(pmin(S,c(rmultinom(1, introduction, S/N))))
  }else
  {
    #distributed evenly over itypes and take minimum to prevent more introductions than S
    return(pmin(S,round(introduction * state$S/sum(state$N)))) 
  }
    
}


#simulation function with build-up and waning of immunity
sim.multitypeSEIR_tleap_buildupandwaning <- function(param.list, seed = NULL){
  #check if both build-up and waning are finite  
  with(c(param.list),{
    set.seed(seed)
    
    steps <- ceiling(length.round / deltat) + 1
    n_steps <- steps * runs
    
    #pre-allocate storage for each output variable (faster than using rbind)
    output <- list(
      time = numeric(n_steps),
      run = integer(n_steps),
      S = matrix(0, nrow = n_steps, ncol = itypes),
      L = matrix(0, nrow = n_steps, ncol = itypes),
      I = matrix(0, nrow = n_steps, ncol = itypes),
      R = matrix(0, nrow = n_steps, ncol = itypes),
      DS = matrix(0, nrow = n_steps, ncol = itypes),
      DL = matrix(0, nrow = n_steps, ncol = itypes),
      DI = matrix(0, nrow = n_steps, ncol = itypes),
      DR = matrix(0, nrow = n_steps, ncol = itypes),
      N = numeric(n_steps),
      Egg_H_farm = numeric(n_steps),
      Egg_I_farm = numeric(n_steps),
      Egg_Total_farm = numeric(n_steps),
      Total_Egg_H_laid = numeric(n_steps),
      Total_Egg_I_laid = numeric(n_steps),
      Egg_H_shipped = numeric(n_steps),
      Egg_I_shipped = numeric(n_steps)
    )
    
    #Egg : Flock production function governing percentage of birds laying eggs per day at time t
    flock_production <- function(t, age_at_d0) {
      t_shift <- (t + age_at_d0)/7 #function assumes week 0 = birth, model uses day 0 = first day of maturity
      a <- 103
      b <- 0.0016
      c <- 1.16
      d <- 20.75
      ((a * exp(-b * t_shift)) / (1 + exp(-c * (t_shift - d))))/100
    }
    
    #gamma distributed time until build-up of immunity function - notice that if trans.var == (trans.mean^2) this distribution is an exponential distribution
    buildup.distribution <- function(t){
      #gamma distribution giving the probability of waning at that time point
      shape_buildup  <- (trans.mean.buildup^2)/trans.var.buildup
      scale_buildup  <- trans.var.buildup/trans.mean.buildup
      
      biological_age <- t + age_at_d0 
      
      if(is.nan(pgamma(biological_age+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE))) {
        return(0)
      }
      if(pgamma(biological_age+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)-(pgamma(biological_age, shape = shape_buildup, scale = scale_buildup, lower.tail = TRUE)) == 0){return(0)}
      #probability of transition is probability changed at time t + delta - probability at time t given that you did not change yet 
      round((pgamma(biological_age+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)-(pgamma(biological_age, shape = shape_buildup, scale = scale_buildup, lower.tail = TRUE)))/(1.-pgamma(biological_age, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)),digits = 25)
    }
    
    #gamma distributed time until waning function - notice that if trans.var == (trans.mean^2) this distribution is an exponential distribution
    waning.distribution <- function(t){
      #gamma distribution giving the probability of waning at that time point
      shape_wane  <- (trans.mean.wane^2)/trans.var.wane
      scale_wane  <- trans.var.wane/trans.mean.wane
      
      biological_age <- t + age_at_d0 
      
      round((pgamma(biological_age, shape = shape_wane , scale = scale_wane, lower.tail = FALSE)-(pgamma(biological_age+deltat, shape = shape_wane , scale = scale_wane, lower.tail = FALSE)))/(pgamma(biological_age, shape = shape_wane, scale = scale_wane, lower.tail = FALSE)),digits = 15)
      
    }
    
    #initial number of high titre given time of introduction taking into account build-up and waning of immunity
    #adjust p.protect with waning
    shape_buildup  <- (trans.mean.buildup^2)/trans.var.buildup
    scale_buildup <- trans.var.buildup/trans.mean.buildup
    #fraction of those that will be protected p.protect that is already protected
    p.protect.adjust <- p.protect*(pgamma(age_at_d0, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE));
    
    #adjust p.protect with waning
    shape_wane  <- (trans.mean.wane^2)/trans.var.wane
    scale_wane  <- trans.var.wane/trans.mean.wane
    p.protect.adjust <- p.protect.adjust*pgamma(age_at_d0, shape = shape_wane , scale = scale_wane, lower.tail = FALSE)
    
    #Initial size
    L0 <- round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0) #number of initially latently infected 
    I0 <- round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0) #number of initially infectious
    R0 <- c(0,0) # number of initially recovered
    S0 <- round(c(1-p.protect.adjust,p.protect.adjust)*(N0),digits = 0)#initially susceptible
    Egg_H0_farm <- 0
    Egg_I0_farm <- 0
    Egg_Total0_farm <- 0
    Total_Egg_H_laid0 <- 0
    Total_Egg_I_laid0 <- 0
    Egg_H_shipped0 <- 0
    Egg_I_shipped0 <- 0
    
    # start_time <- Sys.time()
    # print(paste("Simulation started at:", start_time))
    currun <- 1;
    step <- 1
    
    #calculate the transition rates in the latency and infectious states
    lLat <- k.latency/latency.period
    lInf <- k.infectious/infectious.period
    
    #run the simulations
    while(currun <= runs)  {
      #initialize the state of the system with k.latency and k.infectious states
      #set state of the system with initial values
      state <- list(time = 0, N = N0,C=0, run = currun, dt = 0)
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
      #Egg
      state$Egg_H_farm <- Egg_H0_farm
      state$Egg_I_farm <- Egg_I0_farm
      state$Egg_Total_farm <- Egg_Total0_farm
      state$Total_Egg_H_laid <- Total_Egg_H_laid0
      state$Total_Egg_I_laid <- Total_Egg_I_laid0
      state$Egg_H_shipped <- Egg_H_shipped0
      state$Egg_I_shipped <- Egg_I_shipped0
      
      #marker that the infection has not yet been introduced. This is to prevent that rounding does not trigger the introduction
      introduction.occurred <- FALSE; 
      #one run until maximum time or no infected animals left or the introduction has not yet occurred:
      #while((state$time < length.round & (sum(state$L)+sum(state$I)>0 || !introduction.occurred))) 
      #one run until maximum time or the introduction has not yet occurred:
      while((state$time < length.round || !introduction.occurred) && state$N > 0)
      {
        #check if the infection is introduced and if the time of introduction has passed
        if(!introduction.occurred & state$time >= intro.time)
        {
          #set no.introduction S animals to L
          introStoL <- StoLintro_fun(state$S, state$N, no.introduction, introduction.at.random)
          
          #subtract the introductions from the susceptibles
          state$S <- state$S - introStoL;
          #add to latent
          state$L[,1]<- introStoL;
          #set introduction as done.
          introduction.occurred = TRUE
        }
        
        output$time[step] <- state$time
        output$run[step] <- currun
        output$S[step, ] <- state$S
        output$L[step, ] <- rowSums(state$L)
        output$I[step, ] <- rowSums(state$I)
        output$R[step, ] <- state$R
        output$DS[step, ] <- state$DS
        output$DL[step, ] <- state$DL
        output$DI[step, ] <- state$DI
        output$DR[step, ] <- state$DR
        output$N[step] <- state$N
        output$Egg_H_farm[step] <- state$Egg_H_farm
        output$Egg_I_farm[step] <- state$Egg_I_farm
        output$Egg_Total_farm[step] <- state$Egg_Total_farm
        output$Total_Egg_H_laid[step] <- state$Total_Egg_H_laid
        output$Total_Egg_I_laid[step] <- state$Total_Egg_I_laid
        output$Egg_H_shipped[step] <- state$Egg_H_shipped
        output$Egg_I_shipped[step] <- state$Egg_I_shipped
        
        
        #Transitions between titre status ####
        transProb <-c(0,waning.distribution(state$time),
                      buildup.distribution(state$time),0)
        vacState <- matrix(mapply(rbinom, MoreArgs = list(n = 1), 
                                  size = c(state$S,state$S),prob = transProb), 
                           ncol = itypes)
        #off diagonals are transition form high to low or low to high titre
        diag(vacState) <- state$S - rowSums(vacState) #S remaining in the current state
        #add all transition up
        state$S <- colSums(vacState)
        
        #calculate the force of infection ####
        foi <- (beta %*% rowSums(state$I))/state$N
        
        #Susceptible: 
        #total number of S that leave (from infection or death)
        transS <- rbinom(n = itypes, size = state$S, prob = 1 - exp(-(foi + mortRate) * deltat))
        #of those determine how many leave because of death:
        mortS <- rbinom(n = itypes, size = transS, prob = mortRate / (foi + mortRate + 10^-100))
        #remaining to latent:
        StoL  <- transS - mortS
        #Latent:
        #for each of the k.latency stages draw the number of transitions
        transL <- matrix(rbinom(n = itypes * k.latency, size = as.vector(state$L), prob = 1 - exp(-(lLat + mortRate) * deltat)), nrow = itypes)
        #of these transitions, determine how many are due to mortality:
        mortL <- matrix(rbinom(n = length(transL), size = as.vector(transL), prob = mortRate/(lLat + mortRate)), nrow = itypes)
        #combine the new latent infections coming from susceptibles with those progressing within L:
        LtoNext <- cbind(StoL, transL - mortL)
        #Infectious:
        transI <- matrix(rbinom(n = itypes * k.infectious, size = as.vector(state$I), prob = 1 - exp(-(lInf + mortRate) * deltat)), nrow = itypes)
        mortI <- matrix(rbinom(n = length(transI), size = as.vector(transI), prob = mortRate/(lInf + mortRate)), nrow = itypes)
        #individuals that finish the infectious period (and survive) move on:
        ItoNext <- cbind(LtoNext[, ncol(LtoNext)], transI - mortI)
        #Recovered:
        mortR <- rbinom(n = itypes, size = state$R, prob = 1 - exp(-mortRate * deltat))
        
        #deal with changes
        state$S <- state$S - transS
        state$DS <-state$DS + mortS
        state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
        state$DL <- state$DL + rowSums(mortL)
        state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)]
        
        #those dying at recovery
        mortIdisease <- rbinom(rep(1, itypes),ItoNext[,dim(ItoNext)[2]] ,pdie)
        state$DR <- state$DR + mortR
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        
        # Eggs
        prop_laying <- flock_production(state$time, param.list$age_at_d0)
        new_healthy_eggs <- (sum(state$S) + sum(state$L) + sum(state$R)) * (eh * prop_laying * deltat)
        new_infected_eggs <- sum(state$I) * (ei * prop_laying * deltat) * (1-disfigured)
        
        # Update cumulative eggs laid (does not reset every pickup_day)
        state$Total_Egg_H_laid <- state$Total_Egg_H_laid + new_healthy_eggs
        state$Total_Egg_I_laid <- state$Total_Egg_I_laid + new_infected_eggs
        
        # Update eggs on the farm (i.e. eggs currently held on-farm)
        state$Egg_H_farm <- state$Egg_H_farm + new_healthy_eggs
        state$Egg_I_farm <- (state$Egg_I_farm + new_infected_eggs) 
        state$Egg_Total_farm <- state$Egg_H_farm + state$Egg_I_farm
        
        # At pickup time: update shipments and reset farm eggs simultaneously
        if (state$time > 0 && abs(state$time %% pickup_time) < deltat) {
          # Update shipped eggs to include all eggs laid so far
          state$Egg_H_shipped <- state$Total_Egg_H_laid
          state$Egg_I_shipped <- state$Total_Egg_I_laid
          # Reset farm egg counts
          state$Egg_H_farm <- 0
          state$Egg_I_farm <- 0
          state$Egg_Total_farm <- 0
        } else {
          # Otherwise, calculate shipments as before
          state$Egg_H_shipped <- state$Total_Egg_H_laid - state$Egg_H_farm
          state$Egg_I_shipped <- state$Total_Egg_I_laid - state$Egg_I_farm
        }
        
        
        #update time
        state$time = state$time + deltat
        step <- step + 1
        
      }
      currun <- currun +1
      
    }
    
    # end_time <- Sys.time()
    # print(paste("Simulation ended at:", end_time))
    # 
    # runtime <- end_time - start_time
    # runtime_secs <- as.numeric(runtime, units = "secs")
    # runtime_hours <- floor(runtime_secs / 3600) 
    # runtime_mins <- floor((runtime_secs %% 3600) / 60) 
    # 
    # print(paste("Total runtime:", runtime_hours, "hours and", runtime_mins, "minutes"))
    
    output <- lapply(output, function(x) {
      if (is.matrix(x)) {
        x[1:(step - 1), , drop = FALSE]
      } else {
        x[1:(step - 1)]
      }
    })
    
    return(output)
  })
} 
#simulation function only build-up of immunity
sim.multitypeSEIR_tleap_onlybuildup <- function(param.list, seed = NULL){
  with(c(param.list),{
    set.seed(seed)
    
    steps <- ceiling(length.round / deltat) + 1
    n_steps <- steps * runs
    
    #pre-allocate storage for each output variable (faster than using rbind)
    output <- list(
      time = numeric(n_steps),
      run = integer(n_steps),
      S = matrix(0, nrow = n_steps, ncol = itypes),
      L = matrix(0, nrow = n_steps, ncol = itypes),
      I = matrix(0, nrow = n_steps, ncol = itypes),
      R = matrix(0, nrow = n_steps, ncol = itypes),
      DS = matrix(0, nrow = n_steps, ncol = itypes),
      DL = matrix(0, nrow = n_steps, ncol = itypes),
      DI = matrix(0, nrow = n_steps, ncol = itypes),
      DR = matrix(0, nrow = n_steps, ncol = itypes),
      N = numeric(n_steps),
      Egg_H_farm = numeric(n_steps),
      Egg_I_farm = numeric(n_steps),
      Egg_Total_farm = numeric(n_steps),
      Total_Egg_H_laid = numeric(n_steps),
      Total_Egg_I_laid = numeric(n_steps),
      Egg_H_shipped = numeric(n_steps),
      Egg_I_shipped = numeric(n_steps)
    )
    
    #Egg : Flock production function governing percentage of birds laying eggs per day at time t
    flock_production <- function(t, age_at_d0) {
      t_shift <- (t + age_at_d0)/7 #function assumes week 0 = birth, model uses day 0 = first day of maturity
      a <- 103
      b <- 0.0016
      c <- 1.16
      d <- 20.75
      ((a * exp(-b * t_shift)) / (1 + exp(-c * (t_shift - d))))/100
    }
    
    #gamma distributed time until build-up of immunity function - notice that if trans.var == (trans.mean^2) this distribution is an exponential distribution
    buildup.distribution <- function(t){
      #gamma distribution giving the probability of waning at that time point
      shape_buildup  <- (trans.mean.buildup^2)/trans.var.buildup
      scale_buildup  <- trans.var.buildup/trans.mean.buildup
      
      biological_age <- t + age_at_d0 
      
      if(is.nan(pgamma(biological_age+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE))) {
        return(0)
      }
      if(pgamma(biological_age+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)-(pgamma(biological_age, shape = shape_buildup, scale = scale_buildup, lower.tail = TRUE)) == 0){return(0)}
      #probability of transition is probability changed at time t + delta - probability at time t given that you did not change yet 
      round((pgamma(biological_age+deltat, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)-(pgamma(biological_age, shape = shape_buildup, scale = scale_buildup, lower.tail = TRUE)))/(1.-pgamma(biological_age, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE)),digits = 25)
    }
    
    
    
    #initial number of high titre given time of introduction taking into account build-up and waning of immunity
    #adjust p.protect with waning
    shape_buildup  <- (trans.mean.buildup^2)/trans.var.buildup
    scale_buildup <- trans.var.buildup/trans.mean.buildup
    #fraction of those that will be protected p.protect that is already protected
    p.protect.adjust <- p.protect*(pgamma(age_at_d0, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE));
    
    
    #Initial size
    L0 <- round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0) #number of initially latently infected 
    I0 <- round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0) #number of initially infectious
    R0 <- c(0,0) # number of initially recovered
    S0 <- round(c(1-p.protect.adjust,p.protect.adjust)*(N0),digits = 0)#initially susceptible
    Egg_H0_farm <- 0
    Egg_I0_farm <- 0
    Egg_Total0_farm <- 0
    Total_Egg_H_laid0 <- 0
    Total_Egg_I_laid0 <- 0
    Egg_H_shipped0 <- 0
    Egg_I_shipped0 <- 0
    
    # start_time <- Sys.time()
    # print(paste("Simulation started at:", start_time))
    currun <- 1;
    step <- 1
    
    #calculate the transition rates in the latency and infectious states
    lLat <- k.latency/latency.period
    lInf <- k.infectious/infectious.period
    
    #run the simulations
    while(currun <= runs)  {
      #initialize the state of the system with k.latency and k.infectious states
      #set state of the system with initial values
      state <- list(time = 0, N = N0,C=0, run = currun, dt = 0)
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
      #Egg
      state$Egg_H_farm <- Egg_H0_farm
      state$Egg_I_farm <- Egg_I0_farm
      state$Egg_Total_farm <- Egg_Total0_farm
      state$Total_Egg_H_laid <- Total_Egg_H_laid0
      state$Total_Egg_I_laid <- Total_Egg_I_laid0
      state$Egg_H_shipped <- Egg_H_shipped0
      state$Egg_I_shipped <- Egg_I_shipped0
      
      #marker that the infection has not yet been introduced. This is to prevent that rounding does not trigger the introduction
      introduction.occurred <- FALSE; 
      #one run until maximum time or no infected animals left or the introduction has not yet occurred:
      #while((state$time < length.round & (sum(state$L)+sum(state$I)>0 || !introduction.occurred))) 
      #one run until maximum time or the introduction has not yet occurred:
      while((state$time < length.round || !introduction.occurred) && state$N > 0)
      {
        #check if the infection is introduced and if the time of introduction has passed
        if(!introduction.occurred & state$time >= intro.time)
        {
          #set no.introduction S animals to L
          introStoL <- StoLintro_fun(state$S, state$N, no.introduction, introduction.at.random)
          
          #subtract the introductions from the susceptibles
          state$S <- state$S - introStoL;
          #add to latent
          state$L[,1]<- introStoL;
          #set introduction as done.
          introduction.occurred = TRUE
        }
        
        output$time[step] <- state$time
        output$run[step] <- currun
        output$S[step, ] <- state$S
        output$L[step, ] <- rowSums(state$L)
        output$I[step, ] <- rowSums(state$I)
        output$R[step, ] <- state$R
        output$DS[step, ] <- state$DS
        output$DL[step, ] <- state$DL
        output$DI[step, ] <- state$DI
        output$DR[step, ] <- state$DR
        output$N[step] <- state$N
        output$Egg_H_farm[step] <- state$Egg_H_farm
        output$Egg_I_farm[step] <- state$Egg_I_farm
        output$Egg_Total_farm[step] <- state$Egg_Total_farm
        output$Total_Egg_H_laid[step] <- state$Total_Egg_H_laid
        output$Total_Egg_I_laid[step] <- state$Total_Egg_I_laid
        output$Egg_H_shipped[step] <- state$Egg_H_shipped
        output$Egg_I_shipped[step] <- state$Egg_I_shipped
        
        
        #Transitions between titre status ####
        transProb <-c(0,0,
                      buildup.distribution(state$time),0)
        vacState <- matrix(mapply(rbinom, MoreArgs = list(n = 1), 
                                  size = c(state$S,state$S),prob = transProb), 
                           ncol = itypes)
        #off diagonals are transition form high to low or low to high titre
        diag(vacState) <- state$S - rowSums(vacState) #S remaining in the current state
        #add all transition up
        state$S <- colSums(vacState)
        
        #calculate the force of infection ####
        foi <- (beta %*% rowSums(state$I))/state$N
        
        #Susceptible: 
        #total number of S that leave (from infection or death)
        transS <- rbinom(n = itypes, size = state$S, prob = 1 - exp(-(foi + mortRate) * deltat))
        #of those determine how many leave because of death:
        mortS <- rbinom(n = itypes, size = transS, prob = mortRate / (foi + mortRate + 10^-100))
        #remaining to latent:
        StoL  <- transS - mortS
        #Latent:
        #for each of the k.latency stages draw the number of transitions
        transL <- matrix(rbinom(n = itypes * k.latency, size = as.vector(state$L), prob = 1 - exp(-(lLat + mortRate) * deltat)), nrow = itypes)
        #of these transitions, determine how many are due to mortality:
        mortL <- matrix(rbinom(n = length(transL), size = as.vector(transL), prob = mortRate/(lLat + mortRate)), nrow = itypes)
        #combine the new latent infections coming from susceptibles with those progressing within L:
        LtoNext <- cbind(StoL, transL - mortL)
        #Infectious:
        transI <- matrix(rbinom(n = itypes * k.infectious, size = as.vector(state$I), prob = 1 - exp(-(lInf + mortRate) * deltat)), nrow = itypes)
        mortI <- matrix(rbinom(n = length(transI), size = as.vector(transI), prob = mortRate/(lInf + mortRate)), nrow = itypes)
        #individuals that finish the infectious period (and survive) move on:
        ItoNext <- cbind(LtoNext[, ncol(LtoNext)], transI - mortI)
        #Recovered:
        mortR <- rbinom(n = itypes, size = state$R, prob = 1 - exp(-mortRate * deltat))
        
        #deal with changes
        state$S <- state$S - transS
        state$DS <-state$DS + mortS
        state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
        state$DL <- state$DL + rowSums(mortL)
        state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)]
        
        #those dying at recovery
        mortIdisease <- rbinom(rep(1, itypes),ItoNext[,dim(ItoNext)[2]] ,pdie)
        state$DR <- state$DR + mortR
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        
        # Eggs
        prop_laying <- flock_production(state$time, param.list$age_at_d0)
        new_healthy_eggs <- (sum(state$S) + sum(state$L) + sum(state$R)) * (eh * prop_laying * deltat)
        new_infected_eggs <- sum(state$I) * (ei * prop_laying * deltat) * (1-disfigured)
        
        # Update cumulative eggs laid (does not reset every pickup_day)
        state$Total_Egg_H_laid <- state$Total_Egg_H_laid + new_healthy_eggs
        state$Total_Egg_I_laid <- state$Total_Egg_I_laid + new_infected_eggs
        
        # Update eggs on the farm (i.e. eggs currently held on-farm)
        state$Egg_H_farm <- state$Egg_H_farm + new_healthy_eggs
        state$Egg_I_farm <- (state$Egg_I_farm + new_infected_eggs) 
        state$Egg_Total_farm <- state$Egg_H_farm + state$Egg_I_farm
        
        # At pickup time: update shipments and reset farm eggs simultaneously
        if (state$time > 0 && abs(state$time %% pickup_time) < deltat) {
          # Update shipped eggs to include all eggs laid so far
          state$Egg_H_shipped <- state$Total_Egg_H_laid
          state$Egg_I_shipped <- state$Total_Egg_I_laid
          # Reset farm egg counts
          state$Egg_H_farm <- 0
          state$Egg_I_farm <- 0
          state$Egg_Total_farm <- 0
        } else {
          # Otherwise, calculate shipments as before
          state$Egg_H_shipped <- state$Total_Egg_H_laid - state$Egg_H_farm
          state$Egg_I_shipped <- state$Total_Egg_I_laid - state$Egg_I_farm
        }
        
        
        #update time
        state$time = state$time + deltat
        step <- step + 1
        
      }
      currun <- currun +1
      
    }
    
    # end_time <- Sys.time()
    # print(paste("Simulation ended at:", end_time))
    # 
    # runtime <- end_time - start_time
    # runtime_secs <- as.numeric(runtime, units = "secs")
    # runtime_hours <- floor(runtime_secs / 3600) 
    # runtime_mins <- floor((runtime_secs %% 3600) / 60) 
    # 
    # print(paste("Total runtime:", runtime_hours, "hours and", runtime_mins, "minutes"))
    
    output <- lapply(output, function(x) {
      if (is.matrix(x)) {
        x[1:(step - 1), , drop = FALSE]
      } else {
        x[1:(step - 1)]
      }
    })
    
    return(output)
  })
} 

#simulation function only waning of immunity
sim.multitypeSEIR_tleap_onlywaning <- function(param.list, seed = NULL){
  with(c(param.list),{
    set.seed(seed)
    
    steps <- ceiling(length.round / deltat) + 1
    n_steps <- steps * runs
    
    #pre-allocate storage for each output variable (faster than using rbind)
    output <- list(
      time = numeric(n_steps),
      run = integer(n_steps),
      S = matrix(0, nrow = n_steps, ncol = itypes),
      L = matrix(0, nrow = n_steps, ncol = itypes),
      I = matrix(0, nrow = n_steps, ncol = itypes),
      R = matrix(0, nrow = n_steps, ncol = itypes),
      DS = matrix(0, nrow = n_steps, ncol = itypes),
      DL = matrix(0, nrow = n_steps, ncol = itypes),
      DI = matrix(0, nrow = n_steps, ncol = itypes),
      DR = matrix(0, nrow = n_steps, ncol = itypes),
      N = numeric(n_steps),
      Egg_H_farm = numeric(n_steps),
      Egg_I_farm = numeric(n_steps),
      Egg_Total_farm = numeric(n_steps),
      Total_Egg_H_laid = numeric(n_steps),
      Total_Egg_I_laid = numeric(n_steps),
      Egg_H_shipped = numeric(n_steps),
      Egg_I_shipped = numeric(n_steps)
    )
    
    #Egg : Flock production function governing percentage of birds laying eggs per day at time t
    flock_production <- function(t, age_at_d0) {
      t_shift <- (t + age_at_d0)/7 #function assumes week 0 = birth, model uses day 0 = first day of maturity
      a <- 103
      b <- 0.0016
      c <- 1.16
      d <- 20.75
      ((a * exp(-b * t_shift)) / (1 + exp(-c * (t_shift - d))))/100
    }
    
    
    #gamma distributed time until waning function - notice that if trans.var == (trans.mean^2) this distribution is an exponential distribution
    waning.distribution <- function(t){
      #gamma distribution giving the probability of waning at that time point
      shape_wane  <- (trans.mean.wane^2)/trans.var.wane
      scale_wane  <- trans.var.wane/trans.mean.wane
      
      biological_age <- t + age_at_d0 
      
      round((pgamma(biological_age, 
                    shape = shape_wane, 
                    scale = scale_wane, lower.tail = FALSE)-(pgamma(biological_age+deltat, shape = shape_wane , scale = scale_wane, lower.tail = FALSE)))/(pgamma(biological_age, shape = shape_wane, scale = scale_wane, lower.tail = FALSE)),digits = 15)
      
    }
    
    #initial number of high titre given time of introduction taking into account build-up and waning of immunity
    #adjust p.protect with waning
    shape_wane  <- (trans.mean.wane^2)/trans.var.wane
    scale_wane  <- trans.var.wane/trans.mean.wane
    p.protect.adjust <- p.protect*pgamma(age_at_d0, shape = shape_wane , scale = scale_wane, lower.tail = FALSE)
    
    #Initial size
    L0 <- round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0) #number of initially latently infected 
    I0 <- round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0) #number of initially infectious
    R0 <- c(0,0) # number of initially recovered
    S0 <- round(c(1-p.protect.adjust,p.protect.adjust)*(N0),digits = 0)#initially susceptible
    Egg_H0_farm <- 0
    Egg_I0_farm <- 0
    Egg_Total0_farm <- 0
    Total_Egg_H_laid0 <- 0
    Total_Egg_I_laid0 <- 0
    Egg_H_shipped0 <- 0
    Egg_I_shipped0 <- 0
    
    # start_time <- Sys.time()
    # print(paste("Simulation started at:", start_time))
    currun <- 1;
    step <- 1
    
    #calculate the transition rates in the latency and infectious states
    lLat <- k.latency/latency.period
    lInf <- k.infectious/infectious.period
    
    #run the simulations
    while(currun <= runs)  {
      #initialize the state of the system with k.latency and k.infectious states
      #set state of the system with initial values
      state <- list(time = 0, N = N0,C=0, run = currun, dt = 0)
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
      #Egg
      state$Egg_H_farm <- Egg_H0_farm
      state$Egg_I_farm <- Egg_I0_farm
      state$Egg_Total_farm <- Egg_Total0_farm
      state$Total_Egg_H_laid <- Total_Egg_H_laid0
      state$Total_Egg_I_laid <- Total_Egg_I_laid0
      state$Egg_H_shipped <- Egg_H_shipped0
      state$Egg_I_shipped <- Egg_I_shipped0
      
      #marker that the infection has not yet been introduced. This is to prevent that rounding does not trigger the introduction
      introduction.occurred <- FALSE; 
      #one run until maximum time or no infected animals left or the introduction has not yet occurred:
      #while((state$time < length.round & (sum(state$L)+sum(state$I)>0 || !introduction.occurred))) 
      #one run until maximum time or the introduction has not yet occurred:
      while((state$time < length.round || !introduction.occurred) && state$N > 0)
      {
        #check if the infection is introduced and if the time of introduction has passed
        if(!introduction.occurred & state$time >= intro.time)
        {
          #set no.introduction S animals to L
          introStoL <- StoLintro_fun(state$S, state$N, no.introduction, introduction.at.random)
          
          #subtract the introductions from the susceptibles
          state$S <- state$S - introStoL;
          #add to latent
          state$L[,1]<- introStoL;
          #set introduction as done.
          introduction.occurred = TRUE
        }
        
        output$time[step] <- state$time
        output$run[step] <- currun
        output$S[step, ] <- state$S
        output$L[step, ] <- rowSums(state$L)
        output$I[step, ] <- rowSums(state$I)
        output$R[step, ] <- state$R
        output$DS[step, ] <- state$DS
        output$DL[step, ] <- state$DL
        output$DI[step, ] <- state$DI
        output$DR[step, ] <- state$DR
        output$N[step] <- state$N
        output$Egg_H_farm[step] <- state$Egg_H_farm
        output$Egg_I_farm[step] <- state$Egg_I_farm
        output$Egg_Total_farm[step] <- state$Egg_Total_farm
        output$Total_Egg_H_laid[step] <- state$Total_Egg_H_laid
        output$Total_Egg_I_laid[step] <- state$Total_Egg_I_laid
        output$Egg_H_shipped[step] <- state$Egg_H_shipped
        output$Egg_I_shipped[step] <- state$Egg_I_shipped
        
        
        #Transitions between titre status ####
        transProb <-c(0,waning.distribution(state$time),
                      0,0)
        vacState <- matrix(mapply(rbinom, MoreArgs = list(n = 1), 
                                  size = c(state$S,state$S),prob = transProb), 
                           ncol = itypes)
        #off diagonals are transition form high to low or low to high titre
        diag(vacState) <- state$S - rowSums(vacState) #S remaining in the current state
        #add all transition up
        state$S <- colSums(vacState)
        
        #calculate the force of infection ####
        foi <- (beta %*% rowSums(state$I))/state$N
        
        #Susceptible: 
        #total number of S that leave (from infection or death)
        transS <- rbinom(n = itypes, size = state$S, prob = 1 - exp(-(foi + mortRate) * deltat))
        #of those determine how many leave because of death:
        mortS <- rbinom(n = itypes, size = transS, prob = mortRate / (foi + mortRate + 10^-100))
        #remaining to latent:
        StoL  <- transS - mortS
        #Latent:
        #for each of the k.latency stages draw the number of transitions
        transL <- matrix(rbinom(n = itypes * k.latency, size = as.vector(state$L), prob = 1 - exp(-(lLat + mortRate) * deltat)), nrow = itypes)
        #of these transitions, determine how many are due to mortality:
        mortL <- matrix(rbinom(n = length(transL), size = as.vector(transL), prob = mortRate/(lLat + mortRate)), nrow = itypes)
        #combine the new latent infections coming from susceptibles with those progressing within L:
        LtoNext <- cbind(StoL, transL - mortL)
        #Infectious:
        transI <- matrix(rbinom(n = itypes * k.infectious, size = as.vector(state$I), prob = 1 - exp(-(lInf + mortRate) * deltat)), nrow = itypes)
        mortI <- matrix(rbinom(n = length(transI), size = as.vector(transI), prob = mortRate/(lInf + mortRate)), nrow = itypes)
        #individuals that finish the infectious period (and survive) move on:
        ItoNext <- cbind(LtoNext[, ncol(LtoNext)], transI - mortI)
        #Recovered:
        mortR <- rbinom(n = itypes, size = state$R, prob = 1 - exp(-mortRate * deltat))
        
        #deal with changes
        state$S <- state$S - transS
        state$DS <-state$DS + mortS
        state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
        state$DL <- state$DL + rowSums(mortL)
        state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)]
        
        #those dying at recovery
        mortIdisease <- rbinom(rep(1, itypes),ItoNext[,dim(ItoNext)[2]] ,pdie)
        state$DR <- state$DR + mortR
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        
        # Eggs
        prop_laying <- flock_production(state$time, param.list$age_at_d0)
        new_healthy_eggs <- (sum(state$S) + sum(state$L) + sum(state$R)) * (eh * prop_laying * deltat)
        new_infected_eggs <- sum(state$I) * (ei * prop_laying * deltat) * (1-disfigured)
        
        # Update cumulative eggs laid (does not reset every pickup_day)
        state$Total_Egg_H_laid <- state$Total_Egg_H_laid + new_healthy_eggs
        state$Total_Egg_I_laid <- state$Total_Egg_I_laid + new_infected_eggs
        
        # Update eggs on the farm (i.e. eggs currently held on-farm)
        state$Egg_H_farm <- state$Egg_H_farm + new_healthy_eggs
        state$Egg_I_farm <- (state$Egg_I_farm + new_infected_eggs) 
        state$Egg_Total_farm <- state$Egg_H_farm + state$Egg_I_farm
        
        # At pickup time: update shipments and reset farm eggs simultaneously
        if (state$time > 0 && abs(state$time %% pickup_time) < deltat) {
          # Update shipped eggs to include all eggs laid so far
          state$Egg_H_shipped <- state$Total_Egg_H_laid
          state$Egg_I_shipped <- state$Total_Egg_I_laid
          # Reset farm egg counts
          state$Egg_H_farm <- 0
          state$Egg_I_farm <- 0
          state$Egg_Total_farm <- 0
        } else {
          # Otherwise, calculate shipments as before
          state$Egg_H_shipped <- state$Total_Egg_H_laid - state$Egg_H_farm
          state$Egg_I_shipped <- state$Total_Egg_I_laid - state$Egg_I_farm
        }
        
        
        #update time
        state$time = state$time + deltat
        step <- step + 1
        
      }
      currun <- currun +1
      
    }
    
    # end_time <- Sys.time()
    # print(paste("Simulation ended at:", end_time))
    # 
    # runtime <- end_time - start_time
    # runtime_secs <- as.numeric(runtime, units = "secs")
    # runtime_hours <- floor(runtime_secs / 3600) 
    # runtime_mins <- floor((runtime_secs %% 3600) / 60) 
    # 
    # print(paste("Total runtime:", runtime_hours, "hours and", runtime_mins, "minutes"))
    
    output <- lapply(output, function(x) {
      if (is.matrix(x)) {
        x[1:(step - 1), , drop = FALSE]
      } else {
        x[1:(step - 1)]
      }
    })
    
    return(output)
  })
} 

# simulation function without waning and build-up of immunity
sim.multitypeSEIR_tleap_constantimmunity <- function(param.list, seed = NULL){
  with(c(param.list),{
    set.seed(seed)
    
    steps <- ceiling(length.round / deltat) + 1
    n_steps <- steps * runs
    
    #pre-allocate storage for each output variable (faster than using rbind)
    output <- list(
      time = numeric(n_steps),
      run = integer(n_steps),
      S = matrix(0, nrow = n_steps, ncol = itypes),
      L = matrix(0, nrow = n_steps, ncol = itypes),
      I = matrix(0, nrow = n_steps, ncol = itypes),
      R = matrix(0, nrow = n_steps, ncol = itypes),
      DS = matrix(0, nrow = n_steps, ncol = itypes),
      DL = matrix(0, nrow = n_steps, ncol = itypes),
      DI = matrix(0, nrow = n_steps, ncol = itypes),
      DR = matrix(0, nrow = n_steps, ncol = itypes),
      N = numeric(n_steps),
      Egg_H_farm = numeric(n_steps),
      Egg_I_farm = numeric(n_steps),
      Egg_Total_farm = numeric(n_steps),
      Total_Egg_H_laid = numeric(n_steps),
      Total_Egg_I_laid = numeric(n_steps),
      Egg_H_shipped = numeric(n_steps),
      Egg_I_shipped = numeric(n_steps)
    )
    
    #Egg : Flock production function governing percentage of birds laying eggs per day at time t
    flock_production <- function(t, age_at_d0) {
      t_shift <- (t + age_at_d0)/7 #function assumes week 0 = birth, model uses day 0 = first day of maturity
      a <- 103
      b <- 0.0016
      c <- 1.16
      d <- 20.75
      ((a * exp(-b * t_shift)) / (1 + exp(-c * (t_shift - d))))/100
    }
    
    #Initial size
    L0 <- round(c(1-p.protect,p.protect)*0,digits = 0) #number of initially latently infected 
    I0 <- round(c(1-p.protect,p.protect)*0,digits = 0) #number of initially infectious
    R0 <- c(0,0) # number of initially recovered
    S0 <- round(c(1-p.protect,p.protect)*(N0),digits = 0)#initially susceptible
    Egg_H0_farm <- 0
    Egg_I0_farm <- 0
    Egg_Total0_farm <- 0
    Total_Egg_H_laid0 <- 0
    Total_Egg_I_laid0 <- 0
    Egg_H_shipped0 <- 0
    Egg_I_shipped0 <- 0
    
    # start_time <- Sys.time()
    # print(paste("Simulation started at:", start_time))
    currun <- 1;
    step <- 1
    
    #calculate the transition rates in the latency and infectious states
    lLat <- k.latency/latency.period
    lInf <- k.infectious/infectious.period
    
    #run the simulations
    while(currun <= runs)  {
      #initialize the state of the system with k.latency and k.infectious states
      #set state of the system with initial values
      state <- list(time = 0, N = N0,C=0, run = currun, dt = 0)
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
      #Egg
      state$Egg_H_farm <- Egg_H0_farm
      state$Egg_I_farm <- Egg_I0_farm
      state$Egg_Total_farm <- Egg_Total0_farm
      state$Total_Egg_H_laid <- Total_Egg_H_laid0
      state$Total_Egg_I_laid <- Total_Egg_I_laid0
      state$Egg_H_shipped <- Egg_H_shipped0
      state$Egg_I_shipped <- Egg_I_shipped0
      
      #marker that the infection has not yet been introduced. This is to prevent that rounding does not trigger the introduction
      introduction.occurred <- FALSE; 
      #one run until maximum time or no infected animals left or the introduction has not yet occurred:
      #while((state$time < length.round & (sum(state$L)+sum(state$I)>0 || !introduction.occurred))) 
      #one run until maximum time or the introduction has not yet occurred:
      while((state$time < length.round || !introduction.occurred) && state$N > 0) 
      {
        #check if the infection is introduced and if the time of introduction has passed
        if(!introduction.occurred & state$time >= intro.time)
        {
          #set no.introduction S animals to L
          introStoL <- StoLintro_fun(state$S, state$N, no.introduction, introduction.at.random)
          
          #subtract the introductions from the susceptibles
          state$S <- state$S - introStoL;
          #add to latent
          state$L[,1]<- introStoL;
          #set introduction as done.
          introduction.occurred = TRUE
        }
        
        output$time[step] <- state$time
        output$run[step] <- currun
        output$S[step, ] <- state$S
        output$L[step, ] <- rowSums(state$L)
        output$I[step, ] <- rowSums(state$I)
        output$R[step, ] <- state$R
        output$DS[step, ] <- state$DS
        output$DL[step, ] <- state$DL
        output$DI[step, ] <- state$DI
        output$DR[step, ] <- state$DR
        output$N[step] <- state$N
        output$Egg_H_farm[step] <- state$Egg_H_farm
        output$Egg_I_farm[step] <- state$Egg_I_farm
        output$Egg_Total_farm[step] <- state$Egg_Total_farm
        output$Total_Egg_H_laid[step] <- state$Total_Egg_H_laid
        output$Total_Egg_I_laid[step] <- state$Total_Egg_I_laid
        output$Egg_H_shipped[step] <- state$Egg_H_shipped
        output$Egg_I_shipped[step] <- state$Egg_I_shipped
        
        
        #calculate the force of infection ####
        foi <- (beta %*% rowSums(state$I))/state$N
        
        #Susceptible: 
        #total number of S that leave (from infection or death)
        transS <- rbinom(n = itypes, size = state$S, prob = 1 - exp(-(foi + mortRate) * deltat))
        
        #of those determine how many leave because of death:
        mortS <- rbinom(n = itypes, size = transS, prob = mortRate / (foi + mortRate + 10^-100))
        #remaining to latent:
        StoL  <- transS - mortS
        #Latent:
        #for each of the k.latency stages draw the number of transitions
        transL <- matrix(rbinom(n = itypes * k.latency, size = as.vector(state$L), prob = 1 - exp(-(lLat + mortRate) * deltat)), nrow = itypes)
        #of these transitions, determine how many are due to mortality:
        mortL <- matrix(rbinom(n = length(transL), size = as.vector(transL), prob = mortRate/(lLat + mortRate)), nrow = itypes)
        #combine the new latent infections coming from susceptibles with those progressing within L:
        LtoNext <- cbind(StoL, transL - mortL)
        #Infectious:
        transI <- matrix(rbinom(n = itypes * k.infectious, size = as.vector(state$I), prob = 1 - exp(-(lInf + mortRate) * deltat)), nrow = itypes)
        mortI <- matrix(rbinom(n = length(transI), size = as.vector(transI), prob = mortRate/(lInf + mortRate)), nrow = itypes)
        #individuals that finish the infectious period (and survive) move on:
        ItoNext <- cbind(LtoNext[, ncol(LtoNext)], transI - mortI)
        #Recovered:
        mortR <- rbinom(n = itypes, size = state$R, prob = 1 - exp(-mortRate * deltat))
        
        #deal with changes
        state$S <- state$S - transS
        state$DS <-state$DS + mortS
        state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
        state$DL <- state$DL + rowSums(mortL)
        state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)]
        
        #those dying at recovery
        mortIdisease <- rbinom(rep(1, itypes),ItoNext[,dim(ItoNext)[2]] ,pdie)
        state$DR <- state$DR + mortR
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        
        # Eggs
        prop_laying <- flock_production(state$time, param.list$age_at_d0)
        new_healthy_eggs <- (sum(state$S) + sum(state$L) + sum(state$R)) * (eh * prop_laying * deltat)
        new_infected_eggs <- sum(state$I) * (ei * prop_laying * deltat) * (1-disfigured)
        
        # Update cumulative eggs laid (does not reset every pickup_day)
        state$Total_Egg_H_laid <- state$Total_Egg_H_laid + new_healthy_eggs
        state$Total_Egg_I_laid <- state$Total_Egg_I_laid + new_infected_eggs
        
        # Update eggs on the farm (i.e. eggs currently held on-farm)
        state$Egg_H_farm <- state$Egg_H_farm + new_healthy_eggs
        state$Egg_I_farm <- (state$Egg_I_farm + new_infected_eggs) 
        state$Egg_Total_farm <- state$Egg_H_farm + state$Egg_I_farm
        
        # At pickup time: update shipments and reset farm eggs simultaneously
        if (state$time > 0 && abs(state$time %% pickup_time) < deltat) {
          # Update shipped eggs to include all eggs laid so far
          state$Egg_H_shipped <- state$Total_Egg_H_laid
          state$Egg_I_shipped <- state$Total_Egg_I_laid
          # Reset farm egg counts
          state$Egg_H_farm <- 0
          state$Egg_I_farm <- 0
          state$Egg_Total_farm <- 0
        } else {
          # Otherwise, calculate shipments as before
          state$Egg_H_shipped <- state$Total_Egg_H_laid - state$Egg_H_farm
          state$Egg_I_shipped <- state$Total_Egg_I_laid - state$Egg_I_farm
        }
        
        
        #update time
        state$time = state$time + deltat
        step <- step + 1
        
      }
      currun <- currun +1
      
    }
    
    # end_time <- Sys.time()
    # print(paste("Simulation ended at:", end_time))
    # 
    # runtime <- end_time - start_time
    # runtime_secs <- as.numeric(runtime, units = "secs")
    # runtime_hours <- floor(runtime_secs / 3600) 
    # runtime_mins <- floor((runtime_secs %% 3600) / 60) 
    # 
    # print(paste("Total runtime:", runtime_hours, "hours and", runtime_mins, "minutes"))
    
    output <- lapply(output, function(x) {
      if (is.matrix(x)) {
        x[1:(step - 1), , drop = FALSE]
      } else {
        x[1:(step - 1)]
      }
    })
    
    return(output)
  })
} 



# simulation function with table for high and low titre fractions
sim.multitypeSEIR_tleap_titretable <- function(param.list,
                                               titre_table,
                                               seed = NULL){
  #check if both build-up and waning are finite  
  with(c(param.list),{
    set.seed(seed)
    
    steps <- ceiling(length.round / deltat) + 1
    n_steps <- steps * runs
    
    #pre-allocate storage for each output variable (faster than using rbind)
    output <- list(
      time = numeric(n_steps),
      run = integer(n_steps),
      S = matrix(0, nrow = n_steps, ncol = itypes),
      L = matrix(0, nrow = n_steps, ncol = itypes),
      I = matrix(0, nrow = n_steps, ncol = itypes),
      R = matrix(0, nrow = n_steps, ncol = itypes),
      DS = matrix(0, nrow = n_steps, ncol = itypes),
      DL = matrix(0, nrow = n_steps, ncol = itypes),
      DI = matrix(0, nrow = n_steps, ncol = itypes),
      DR = matrix(0, nrow = n_steps, ncol = itypes),
      N = numeric(n_steps),
      Egg_H_farm = numeric(n_steps),
      Egg_I_farm = numeric(n_steps),
      Egg_Total_farm = numeric(n_steps),
      Total_Egg_H_laid = numeric(n_steps),
      Total_Egg_I_laid = numeric(n_steps),
      Egg_H_shipped = numeric(n_steps),
      Egg_I_shipped = numeric(n_steps)
    )
    
    #Egg : Flock production function governing percentage of birds laying eggs per day at time t
    flock_production <- function(t, age_at_d0) {
      t_shift <- (t + age_at_d0)/7 #function assumes week 0 = birth, model uses day 0 = first day of maturity
      a <- 103
      b <- 0.0016
      c <- 1.16
      d <- 20.75
      ((a * exp(-b * t_shift)) / (1 + exp(-c * (t_shift - d))))/100
    }
    
    #create a table with titre for each times step 
    step_titre_table <- merge(data.frame(time = seq(0, length.round, deltat)), titre_table, by = "time", all.x = TRUE) %>%
      tidyr::fill(value) 
    
    
    #initial number of high titre given time of introduction taking into account build-up and waning of immunity
    #fraction of those that will be protected p.protect that is already protected
    p.protect.adjust <- p.protect*(pgamma(age_at_d0, shape = shape_buildup , scale = scale_buildup, lower.tail = TRUE));
    
    #adjust p.protect with waning
    shape_wane  <- (trans.mean.wane^2)/trans.var.wane
    scale_wane  <- trans.var.wane/trans.mean.wane
    p.protect.adjust <- p.protect.adjust*pgamma(age_at_d0, shape = shape_wane , scale = scale_wane, lower.tail = FALSE)
    
    #Initial size
    L0 <- round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0) #number of initially latently infected 
    I0 <- round(c(1-p.protect.adjust,p.protect.adjust)*0,digits = 0) #number of initially infectious
    R0 <- c(0,0) # number of initially recovered
    S0 <- round(c(1-p.protect.adjust,p.protect.adjust)*(N0),digits = 0)#initially susceptible
    Egg_H0_farm <- 0
    Egg_I0_farm <- 0
    Egg_Total0_farm <- 0
    Total_Egg_H_laid0 <- 0
    Total_Egg_I_laid0 <- 0
    Egg_H_shipped0 <- 0
    Egg_I_shipped0 <- 0
    
    # start_time <- Sys.time()
    # print(paste("Simulation started at:", start_time))
    currun <- 1;
    step <- 1
    
    #calculate the transition rates in the latency and infectious states
    lLat <- k.latency/latency.period
    lInf <- k.infectious/infectious.period
    
    #run the simulations
    while(currun <= runs)  {
      #initialize the state of the system with k.latency and k.infectious states
      #set state of the system with initial values
      state <- list(time = 0, N = N0,C=0, run = currun, dt = 0)
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
      #Egg
      state$Egg_H_farm <- Egg_H0_farm
      state$Egg_I_farm <- Egg_I0_farm
      state$Egg_Total_farm <- Egg_Total0_farm
      state$Total_Egg_H_laid <- Total_Egg_H_laid0
      state$Total_Egg_I_laid <- Total_Egg_I_laid0
      state$Egg_H_shipped <- Egg_H_shipped0
      state$Egg_I_shipped <- Egg_I_shipped0
      
      #marker that the infection has not yet been introduced. This is to prevent that rounding does not trigger the introduction
      introduction.occurred <- FALSE; 
      #one run until maximum time or no infected animals left or the introduction has not yet occurred:
      #while((state$time < length.round & (sum(state$L)+sum(state$I)>0 || !introduction.occurred))) 
      #one run until maximum time or the introduction has not yet occurred:
      while((state$time < length.round || !introduction.occurred) && state$N > 0)
      {
        #check if the infection is introduced and if the time of introduction has passed
        if(!introduction.occurred & state$time >= intro.time)
        {
          #set no.introduction S animals to L
          introStoL <- StoLintro_fun(state$S, state$N, no.introduction, introduction.at.random)
          
          #subtract the introductions from the susceptibles
          state$S <- state$S - introStoL;
          #add to latent
          state$L[,1]<- introStoL;
          #set introduction as done.
          introduction.occurred = TRUE
        }
        
        output$time[step] <- state$time
        output$run[step] <- currun
        output$S[step, ] <- state$S
        output$L[step, ] <- rowSums(state$L)
        output$I[step, ] <- rowSums(state$I)
        output$R[step, ] <- state$R
        output$DS[step, ] <- state$DS
        output$DL[step, ] <- state$DL
        output$DI[step, ] <- state$DI
        output$DR[step, ] <- state$DR
        output$N[step] <- state$N
        output$Egg_H_farm[step] <- state$Egg_H_farm
        output$Egg_I_farm[step] <- state$Egg_I_farm
        output$Egg_Total_farm[step] <- state$Egg_Total_farm
        output$Total_Egg_H_laid[step] <- state$Total_Egg_H_laid
        output$Total_Egg_I_laid[step] <- state$Total_Egg_I_laid
        output$Egg_H_shipped[step] <- state$Egg_H_shipped
        output$Egg_I_shipped[step] <- state$Egg_I_shipped
        
        
        #Transitions between titre status ####
        transProb <-c(0,waning.distribution(state$time),
                      buildup.distribution(state$time),0)
        vacState <- matrix(mapply(rbinom, MoreArgs = list(n = 1), 
                                  size = c(state$S,state$S),prob = transProb), 
                           ncol = itypes)
        #off diagonals are transition form high to low or low to high titre
        diag(vacState) <- state$S - rowSums(vacState) #S remaining in the current state
        #add all transition up
        state$S <- colSums(vacState)
        
        #calculate the force of infection ####
        foi <- (beta %*% rowSums(state$I))/state$N
        
        #Susceptible: 
        #total number of S that leave (from infection or death)
        transS <- rbinom(n = itypes, size = state$S, prob = 1 - exp(-(foi + mortRate) * deltat))
        #of those determine how many leave because of death:
        mortS <- rbinom(n = itypes, size = transS, prob = mortRate / (foi + mortRate + 10^-100))
        #remaining to latent:
        StoL  <- transS - mortS
        #Latent:
        #for each of the k.latency stages draw the number of transitions
        transL <- matrix(rbinom(n = itypes * k.latency, size = as.vector(state$L), prob = 1 - exp(-(lLat + mortRate) * deltat)), nrow = itypes)
        #of these transitions, determine how many are due to mortality:
        mortL <- matrix(rbinom(n = length(transL), size = as.vector(transL), prob = mortRate/(lLat + mortRate)), nrow = itypes)
        #combine the new latent infections coming from susceptibles with those progressing within L:
        LtoNext <- cbind(StoL, transL - mortL)
        #Infectious:
        transI <- matrix(rbinom(n = itypes * k.infectious, size = as.vector(state$I), prob = 1 - exp(-(lInf + mortRate) * deltat)), nrow = itypes)
        mortI <- matrix(rbinom(n = length(transI), size = as.vector(transI), prob = mortRate/(lInf + mortRate)), nrow = itypes)
        #individuals that finish the infectious period (and survive) move on:
        ItoNext <- cbind(LtoNext[, ncol(LtoNext)], transI - mortI)
        #Recovered:
        mortR <- rbinom(n = itypes, size = state$R, prob = 1 - exp(-mortRate * deltat))
        
        #deal with changes
        state$S <- state$S - transS
        state$DS <-state$DS + mortS
        state$L <- state$L - transL + LtoNext[,1:(dim(LtoNext)[2]-1)]
        state$DL <- state$DL + rowSums(mortL)
        state$I <- state$I - transI + ItoNext[,1:(dim(ItoNext)[2]-1)]
        
        #those dying at recovery
        mortIdisease <- rbinom(rep(1, itypes),ItoNext[,dim(ItoNext)[2]] ,pdie)
        state$DR <- state$DR + mortR
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        
        # Eggs
        prop_laying <- flock_production(state$time, param.list$age_at_d0)
        new_healthy_eggs <- (sum(state$S) + sum(state$L) + sum(state$R)) * (eh * prop_laying * deltat)
        new_infected_eggs <- sum(state$I) * (ei * prop_laying * deltat) * (1-disfigured)
        
        # Update cumulative eggs laid (does not reset every pickup_day)
        state$Total_Egg_H_laid <- state$Total_Egg_H_laid + new_healthy_eggs
        state$Total_Egg_I_laid <- state$Total_Egg_I_laid + new_infected_eggs
        
        # Update eggs on the farm (i.e. eggs currently held on-farm)
        state$Egg_H_farm <- state$Egg_H_farm + new_healthy_eggs
        state$Egg_I_farm <- (state$Egg_I_farm + new_infected_eggs) 
        state$Egg_Total_farm <- state$Egg_H_farm + state$Egg_I_farm
        
        # At pickup time: update shipments and reset farm eggs simultaneously
        if (state$time > 0 && abs(state$time %% pickup_time) < deltat) {
          # Update shipped eggs to include all eggs laid so far
          state$Egg_H_shipped <- state$Total_Egg_H_laid
          state$Egg_I_shipped <- state$Total_Egg_I_laid
          # Reset farm egg counts
          state$Egg_H_farm <- 0
          state$Egg_I_farm <- 0
          state$Egg_Total_farm <- 0
        } else {
          # Otherwise, calculate shipments as before
          state$Egg_H_shipped <- state$Total_Egg_H_laid - state$Egg_H_farm
          state$Egg_I_shipped <- state$Total_Egg_I_laid - state$Egg_I_farm
        }
        
        
        #update time
        state$time = state$time + deltat
        step <- step + 1
        
      }
      currun <- currun +1
      
    }
    
    # end_time <- Sys.time()
    # print(paste("Simulation ended at:", end_time))
    # 
    # runtime <- end_time - start_time
    # runtime_secs <- as.numeric(runtime, units = "secs")
    # runtime_hours <- floor(runtime_secs / 3600) 
    # runtime_mins <- floor((runtime_secs %% 3600) / 60) 
    # 
    # print(paste("Total runtime:", runtime_hours, "hours and", runtime_mins, "minutes"))
    
    output <- lapply(output, function(x) {
      if (is.matrix(x)) {
        x[1:(step - 1), , drop = FALSE]
      } else {
        x[1:(step - 1)]
      }
    })
    
    return(output)
  })
} 
