

#outdated: please use faster_multitypetransitionsSEIR_tleap.R




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

#Define number of types
itypes = 2; #Type 1  = not protected by vaccination; Type 2 = protected by vaccination

#List of all parameters
param.list <- list(
  #Simulation settings
  scenario = "Test", 
  runs = 100, #number of runs
  deltat = 0.005, #time step; nb: 0.01 is faster but not accurate, 0.001 is very slow and only mild improvement in accuracy
  #Flock settings
  length.round = 19*30,#length of the round. Simulation starts at time 0
  age_at_d0 = 120, #unit: days; model starts at t = 0 when chickens are 120 days/17 wks old (important for egg production curve and vx buildup)
  N0 = 46000, #population size#### 
  #Vaccination settings
  itypes = itypes, #types####
  p.protect = 1, #proportion of the flock that is vaccinated
  #Introduction settings
  no.introduction= 1, #initially infected at introduction
  intro.time = 0, #time during the round at which the infection is introduced
  #Infection settings
  beta = matrix(c(1.13, 1.13,0.05,0.05), ncol = itypes), #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
  latency.period = c(1,1), #latency period of 1 day
  k.latency =2, #k parameter Erlang distribution
  infectious.period = c(3.0,4.0), #Duration infectious period 
  k.infectious = 20, #k parameter Erlang distribution, k=1 gives sim results comparable to pmaj predictions for exponential 
  trans.mean.wane = 514, #only considers transistion from type 2 to 1; for waning off set to 20000
  trans.var.wane = 185^2, #only considers transistion from type 2 to 1
  trans.mean.buildup = 14, #mean time to build-up immunity
  trans.var.buildup = 2, #variance of time to build-up immunity
  pdie = c(1.0,0.01),#probability of dying at end of infectious period
  mortRate = 0, #0.0005, #per capita death rate for broilers (0.0005 too high for layers) 
  #Egg settings
  eh = 0.57, # daily egg-laying rate for healthy chickens (S, E, R)
  ei = 0.57 * 0.5,# daily egg-laying rate for infected chickens (I) (50% reduction in egg laying)
  disfigured = 0.1, # daily removal rate for infected eggs (10% do not pass inspection)
  pickup_time = 3 #how often eggs are picked up 
)



#function of simulation of spread with transition between states based on gamma distributed time until transition ####
sim.multitypeSEIR_tleap <- function(param.list, seed = NULL){
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
      Egg_H = numeric(n_steps),
      Egg_I = numeric(n_steps),
      Egg_Total = numeric(n_steps),
      FOI = matrix(0, nrow = n_steps, ncol = itypes)
    )
    
    #Egg : Flock production function governing percentage of birds laying eggs per day at time t
    flock_production <- function(t, age_at_d0) {
      t_shift <- (t + age_at_d0)/7 #function assumes week 0 = birth, model uses day 0 = first day of maturity
      a <- 126.65
      b <- 0.01251
      c <- 1.164
      d <- 23.82
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
    Egg_H0 <- c(0,0)
    Egg_I0 <- c(0,0)
    Egg_Total0 <- c(0,0) 
    
    start_time <- Sys.time()
    print(paste("Simulation started at:", start_time))
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
      state$Egg_H <- Egg_H0
      state$Egg_I <- Egg_I0
      state$Egg_Total <- Egg_Total0
      
      #marker that the infection has not yet been introduced. This is to prevent that rounding does not trigger the introduction
      introduction.occurred <- FALSE; 
      #one run until maximum time or no infected animals left or the introduction has not yet occurred
      #while((state$time < length.round & (sum(state$L)+sum(state$I)>0 || !introduction.occurred))) 
      #one run until maximum time or the introduction has not yet occurred
      while(state$time < length.round || !introduction.occurred) 
      {
        #check if the infection is introduced and if the time of introduction has passed
        if(!introduction.occurred & state$time >= intro.time)
        {
          #set no.introduction S animals to L
          #if no.introduction is one number this is done proportional to the types 
          if(length(no.introduction)==1)
          {
            introStoL <- round(no.introduction * state$S/sum(state$N));
            #correct for when the sum of introStoL is less than 1 (ie when p.protect is 0.5 and rounded down to 0) so that the single infection is type 1. 
            if(sum(introStoL)<1)introStoL <- c(1,0)
          }else
            #else use the 
          {
            #get the maximum possible transitions from StoL
            introStoL <- pmin(state$S,no.introduction);
            
          }
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
        output$Egg_H[step] <- sum(state$Egg_H)
        output$Egg_I[step] <- sum(state$Egg_I)
        output$Egg_Total[step] <- sum(state$Egg_Total)
        
        
        #Transitions between titre status ####
        transProb <-c(0,waning.distribution(state$time),buildup.distribution(state$time),0)
        vacState <- matrix(mapply(rbinom, MoreArgs = list(n = 1), size = c(state$S,state$S),prob = transProb), ncol = itypes)
        
        #off diagonals are transition form high to low or low to high titre 
        diag(vacState) <- state$S - rowSums(vacState) #S remaining in the current state
        #add all transition up 
        state$S <- colSums(vacState)
        
        #calculate the force of infection ####
        foi <- (beta %*% rowSums(state$I))/state$N
        output$FOI[step, ] <- foi
        #determine transitions to next state or death ####
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
        state$DR <- state$R - mortR
        state$R <- state$R + ItoNext[,dim(ItoNext)[2]] - mortIdisease - mortR
        state$DI<- state$DI + rowSums(mortI) + mortIdisease
        state$N <- state$N - sum(mortS + rowSums(mortL)  +rowSums(mortI)+mortR + mortIdisease)
        
        #Egg
        prop_laying <- flock_production(state$time, param.list$age_at_d0) * deltat
        new_healthy_eggs <- (sum(state$S) + sum(state$L) + sum(state$R)) * (eh * prop_laying *deltat)
        #new_healthy_eggs <- (state$S[,1] + state$S[,2] + rowSums(state$L)[1] + rowSums(state$L)[2] + state$R[1] + state$R[2]) * eh * prop_laying #tracking by type
        new_infected_eggs <- rowSums(state$I) * (ei * prop_laying *deltat)
        #new_infected_eggs <- (rowSums(state$I)[1] + rowSums(state$I)[2]) * ei * prop_laying #tracking by type
        state$Egg_H <- state$Egg_H + new_healthy_eggs
        state$Egg_I <- (state$Egg_I + new_infected_eggs) * exp(-disfigured)
        #state$Egg_Total <- state$Egg_Total + sum(new_healthy_eggs) + sum(new_infected_eggs)
        state$Egg_Total <- state$Egg_Total + new_healthy_eggs + new_infected_eggs
        
        #Egg : clear every 3 days from pickup
        if ((state$time %% pickup_time) < deltat) {
          state$Egg_H <- c(0,0)
          state$Egg_I <- c(0,0)
          state$Egg_Total <- c(0,0)
        } 
        
        #update time
        state$time = state$time + deltat
        step <- step + 1
        
      }
      currun <- currun +1
      
    }
    
    end_time <- Sys.time()
    print(paste("Simulation ended at:", end_time))
    
    runtime <- end_time - start_time
    runtime_secs <- as.numeric(runtime, units = "secs")
    runtime_hours <- floor(runtime_secs / 3600) 
    runtime_mins <- floor((runtime_secs %% 3600) / 60) 
    
    print(paste("Total runtime:", runtime_hours, "hours and", runtime_mins, "minutes"))
    
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

 
#if(exists("x")){rm(x)};x <-sim.multitypeSEIR_tleap(param.list, seed = NULL)
#if(exists("x")){rm(x)};x <-sim.multitypeSEIR_tleap(param.list, seed = 1440)


# PLOTTING
##Dynamics
# ggplot(data = data.frame(x))+
#   geom_path(aes(x = time, y = S.1 + S.2, group = run, colour = "S")) +
   #geom_path(aes(x = time, y = S.1, group = run, colour = "S.1")) +
   #geom_path(aes(x = time, y = S.2, group = run, colour = "S.2")) +
   #geom_path(aes(x = time, y = L.1 + L.2,group = run, colour = "L"))+
   #geom_path(aes(x = time, y = I.2,group = run, colour = "I"))+
   #geom_path(aes(x = time, y = R.2,group = run, colour = "R"))+
   #geom_path(aes(x = time, y = DR.1 + DR.2, group = run, colour = "DR")) +
   #geom_path(aes(x = time, y = DI.1 + DI.2, group = run, colour = "DI")) +
   #geom_path(aes(x = time, y = DL.1 + DL.2, group = run, colour = "DL")) +
  #  labs(
  #    title = "SEIR Chickens",
  #    subtitle = paste("Population Size (N0):", param.list$N0, "| Introduction Time:", param.list$intro.time),
  #    x = "Time",
  #    y = "Number of Chickens",
  #    colour = "Compartment"
  # )
# 
# 
# ##Final size and number of small outbreaks
# all_runs <- bind_rows(x)
# fs <- data.frame(cbind(rowSums(all_runs$R + all_runs$DI + all_runs$DR + all_runs$DL), all_runs$run)) %>%
#   group_by(X2) %>%
#   reframe(fs = max(X1))
# 
# # Define small outbreaks as outbreaks that are less than 10% of the population size
# small_outbreaks <- sum(fs$fs < 0.1 * param.list$N0)
# total_runs <- param.list$runs
# proportion_small_outbreaks <- small_outbreaks / total_runs
# print(proportion_small_outbreaks)
# 
# ggplot(fs, aes(x = fs)) +
#   geom_histogram(binwidth = 2000, fill = "salmon", color = "black", boundary = 0) +
#   labs(
#     title = "Distribution of Outbreak Sizes",
#     subtitle = paste(
#       "Number of small outbreaks (<10% of population):", small_outbreaks, "out of", total_runs, "runs\n",
#       "Population Size (N0):", param.list$N0, "| Introduction Time:", param.list$intro.time
#     ),
#     x = "Outbreak Size",
#     y = "Frequency"
#   ) +
#   theme_minimal(base_size = 14) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5),
#     axis.title = element_text(face = "bold")
#   )
# 
# 
# 
# ## Can skip: Plot vaccination status incl. key time points ##
#  simulation_data <- data.frame(x)
#  simulation_data$fraction_vaccinated <- (simulation_data$S.2 + simulation_data$L.2 + simulation_data$I.2 + simulation_data$R.2) /
#    (simulation_data$S.1 + simulation_data$S.2 + simulation_data$L.1 + simulation_data$L.2 +
#       simulation_data$I.1 + simulation_data$I.2 + simulation_data$R.1 + simulation_data$R.2)
# 
#  # Define target fractions
#  target_fractions <- c(0.70,0.65)
# 
#  # Find the closest points for each target fraction
#  closest_points <- lapply(target_fractions, function(target) {
#    simulation_data[which.min(abs(simulation_data$fraction_vaccinated - target)), ]
#  })
# # Combine the closest points into a data frame
#  closest_points_df <- do.call(rbind, closest_points)
# 
#  # Plot with points and labels for each target fraction
#  ggplot(data = simulation_data) +
#    geom_path(aes(x = time, y = fraction_vaccinated, group = run)) +
#    ylim(0, 1) +
#    geom_hline(yintercept = target_fractions, linetype = "dashed", color = "red") +
#    geom_point(data = closest_points_df, aes(x = time, y = fraction_vaccinated), color = "black", size = 3) +
#    geom_text(data = closest_points_df, aes(x = time, y = fraction_vaccinated, label = round(time, 2)),
#              vjust = -1, color = "red") +
#    labs(
#      title = "Vaccination Status",
#      subtitle = paste("Population Size (N0):", param.list$N0, "| Introduction Time:", param.list$intro.time),
#      x = "Time",
#     y = "Fraction Vaccinated"
#    )
# 
