#########################################################
#                                                        
#                 Detection module                            
#                                                        
#                  Author: Egil Fischer                               
#                  Contact: e.a.j.fischer@uu.nl                             
#                  Creation date: 28-3-2023                         
#########################################################
source("./src/loadLibraries.R") 


#Surveillance protocols (and notes):

#1. Input: Use the "pre_processing_detModule.R" on simulation/scenario outputs; this "detectionModule" assumes daily data, daily.death.incidence, daily.mortality.rate etc are all available columns 
  #Output: detection times for each protocol, and the number of contaminated eggs that have been shipped and still in holding (on farm) 

#2.Passive detection v3: detection based on a mortality ratio (if the sum of mortality in 2 days is 3x higher than previous week *and/or* 2 day mortality is greater than or equal to 0.5% and confirmation testing)
    #previous passive detection v2: detection based on a mortality ratio (if the sum of weekly mortality is 3x higher than previous week) (see earlier version for code: 18 February 2025)
    #previous passive detection v1: based on 2 days of mortality greater than or equal to 0.5% N (see earlier versions for code)
#Passive detection starting time can be changed by altering "current.init.pas" to start at 0 or to sync with current.init.ac
    
#3. Active surveillance: detection based on active sampling, followed by a confirmation test on all dead birds & a sample of live  birds 
    #easily changed to sample of all live birds by including I.2 in section 5. Scenario Processing 
#Active detection starting time can be changed by altering "current.init.ac" 

#############################
## 0. Detection parameters ##
#############################
picked_sick_input <- c(0,0,0,0,1,1,0,0) #indicates the chance that during a survey an specific type is added to the selection based on clinical signs.
se_input = 1; #sensitivity for dead infected birds
se.survey.live_input = 0.98; #sensitivity for dead infected birds see EfSA report
se.confirm.live_input = 1; #sensitivity for dead infected birds see EfSA report
threshold_input = 0.005; #warning in passive when more than this fraction dies
pfarm_input = 1; #probability that this farm is being monitored
panimal_input = 1; #probability per detectables animal to being monitored,
tlag_input = 1; #time to confirmation
roundTime_input = 1; #rounds to nearest day 
survey.live.sample.size_input = 60;
confirm.live.sample.size_input = 20;
ints_input = 2;#consecutive days for threshold 

#############################
## 1. Passive Detection  ####
#############################

#Ratio approach: detect if a certain number of R dead are found within a certain time interval
passive.detection.ratio <- function(ratio.data, init.pas) {
  days <- ratio.data$day
  detection.time.ratio <- Inf #initialize with Inf
  
  #need at least 7 previous days so start at (init.pas + 7) and end at the day before the last day
  for (d in (init.pas + 7):(max(days) - 1)) {
    index <- which(ratio.data$day == d)
    if (length(index) == 0 || (index + 1) > nrow(ratio.data) || (index - 7) < 1) next
    
    current.2day <- sum(ratio.data$daily.death.incidence[index:(index + 1)], na.rm = TRUE)
    previous.week <- sum(ratio.data$daily.death.incidence[(index - 7):(index - 1)], na.rm = TRUE)
    
    if (!is.na(previous.week) && previous.week > 0 &&
        !is.na(current.2day) && current.2day >= 3 * previous.week) {
      detection.time.ratio <- d
      break 
    }
  }
  
  return(detection.time.ratio)
}

#Passive detection: Threshold approach
passive.detection.time.threshold <- function(threshold.data, threshold, ints, init.pas) {
  # Only consider days >= init.pas
  threshold.data <- threshold.data %>% filter(day >= init.pas) %>%
    #trigger for days where deaths are at or above the threshold
    mutate(exceed = as.numeric(mortality.rate >= threshold))
  
  detection.time.threshold <- Inf #initialize with Inf
  #look for the first occurrence of 'ints' consecutive days with deaths > the threshold
  for (i in 1:(nrow(threshold.data) - ints + 1)) {
    if (sum(threshold.data$exceed[i:(i + ints - 1)]) == ints) {
      detection.time.threshold <- threshold.data$day[i]
      break
    }
  }
  
  return(detection.time.threshold)
}


#############################
## 2. Active Detection   ####
#############################

#Active detection  
active.detection.time<- function(times, survey.times, lumped.detectables, survey.live.detectables, survey.positive.live.vars,confirm.live.detectables, confirm.positive.live.vars,confirm.dead.data,
                      se, se.confirm.live,se.survey.live, pfarm, panimal, tlag, confirm.live.sample.size, survey.live.sample.size,
                      init.ac, roundTime, picked_sick = 0){
  #when there is no active detection by lumped.detectables
  if(nrow(lumped.detectables)==0) lumped.detectables <- tibble(detectables = 0);
  # check that sensitivity values are defined for each detectable
  if (length(se) != ncol(lumped.detectables) && length(se) > 1)
    stop("Not all sensitivities are defined")
  if (length(se.confirm.live) == 1) {
    se.confirm.live <- rep(se.confirm.live, ncol(confirm.live.detectables)-1)
  }
  if (length(se.confirm.live) != ncol(confirm.live.detectables)-1)
    stop("Not all sensitivities are defined for live detectables of confirmation")
  
  if (length(se.survey.live) == 1) {
    se.survey.live <- rep(se.survey.live, length(survey.positive.live.vars))
  }
  if (length(se.survey.live) != length(survey.positive.live.vars))
    stop("Not all sensitivities are defined for live detectables for survey")
  
  # farm may escape detection by active sampling
  if (runif(1) > pfarm)
    return(list(
      #bucket sampling results
      first.ac.detection.time = Inf,
      confirm.result = "NA",
      live.confirmation.result = "NA",
      total.detectables = 0,
      live.positives = 0,
      #live bird survey results
      first.ac.detection.survey.time = Inf,
      confirm.survey.result = "NA",
      live.confirmation.survey.result = "NA",
      total.detectables.survey = 0,
      live.positives.survey = 0
      
      
    ))
  
  # if no detectable animals are or can be found, return Inf detection time
  if (max(as.matrix(lumped.detectables)) == 0 & survey.live.sample.size == 0)
    return(list(#bucket sampling results
      first.ac.detection.time = Inf,
      confirm.result = "NA",
      live.confirmation.result = "NA",
      total.detectables = 0,
      live.positives = 0,
      #live bird survey results
      first.ac.detection.survey.time = Inf,
      confirm.survey.result = "NA",
      live.confirmation.survey.result = "NA",
      total.detectables.survey = 0,
      live.positives.survey = 0    ))
  
  #determine first active detection time using a binomial draw per row ####
  #panimal = 1 assumes that all dead birds are sampled
  #return the first time for which the number of successes is larger than one (if none return Inf)
  detection.flags <- apply(lumped.detectables, 1, function(x) {
    sum(mapply(function(d, s) {
      rbinom(n = 1, size = d, p = panimal * s)
    }, d = x, s = panimal * se))
  })
  
  valid.times <- times[detection.flags >= 1]
  first.ac.detection.time <- if (length(valid.times) > 0) min(valid.times) else Inf
  
  #
  if (is.infinite(first.ac.detection.time)) {
    # return(list(first.ac.detection.time = Inf, confirm.result = "NA", 
    #             live.confirmation.result = "NA", total.detectables = 0, live.positives = NA))
    first.ac.detection.time = Inf; 
    confirm.result = "NA";
    live.confirmation.result = "NA"; 
    total.detectables = 0; 
    live.positives = NA
  }else{
  
    #repeat active surveillance until positive confirmation or no detection
    while(length(valid.times)>=1){
  #Active dead confirmation
  #confirmation test after a lag period to account for test delays: all dead birds
  confirm.test.day <- first.ac.detection.time + tlag
  confirm.dead <- confirm.dead.data %>%
    filter(day > first.ac.detection.time & day <= confirm.test.day) %>%
    summarize(total.detectables = sum(daily.death.incidence.detectables, na.rm = TRUE))
  
  total.detectables <- confirm.dead$total.detectables %||% 0
  confirm.positive <- rbinom(n = 1, size = total.detectables, p = panimal * se) >= 1
  confirm.result <- ifelse(confirm.positive, "POSITIVE", "NEGATIVE")
  
  #Active live confirmation
  #secondary confirmation (independent to dead confirmation): live bird sampling on test day
  day.before <- confirm.test.day - 1
  live.samples <- confirm.live.detectables %>%
    filter(day > day.before & day <= confirm.test.day)
  
  if (nrow(live.samples) == 0) {
    live.confirmation.result <- "NEGATIVE"
    total.live.detectables <- 0
  } else {
    total.live.detectables <- rowSums(select(live.samples, where(is.numeric)), na.rm = TRUE)
    total.live.detectables <- sum(total.live.detectables)
    
    if (total.live.detectables == 0) {
      live.confirmation.result <- "NEGATIVE"
    } else {
      positive.birds <- sum(live.samples %>% select(all_of(confirm.positive.live.vars)) %>% unlist(), na.rm = TRUE)
      live.panimal <- positive.birds / total.live.detectables
      sample.size <- min(confirm.live.sample.size, total.live.detectables)
      live.positives <- rbinom(n = 1, size = sample.size, p = live.panimal * se.confirm.live)
      live.confirmation.result <- ifelse(live.positives >= 1, "POSITIVE", "NEGATIVE")
    }
  }
  if(confirm.result != "POSITIVE" & live.confirmation.result != "POSITIVE" & length(valid.times) >0)
  {
    valid.times = valid.times[-1];
  first.ac.detection.time <- if (length(valid.times) > 0) min(valid.times) else Inf}else{valid.times = c()} 
    }
  }
  
  
  #use survey of live birds ####
  if(survey.live.sample.size>0){
    #
    selection_function = function(x){tabulate(sample(rep(1:length(x), times = x), 
                                                     size = min(survey.live.sample.size,sum(x)), replace = FALSE), nbins = length(x))}
    #get the number of live birds to be sampled during survey per day
    survey_live_detectables <- survey.live.detectables %>%select(-day)%>%
      rowwise() %>%
      reframe(sampled_row = list(selection_function(c_across()))) %>%
      unnest_wider(sampled_row, names_sep = ".") %>%
      setNames(names(survey.live.detectables)[-1]) %>%
      select(any_of(survey.positive.live.vars))
    
    #assume infectious birds will always be found as these are really sick. use pick_sick
    if(length(picked_sick)==1)picked_sick = rep(0, ncol(survey.live.detectables)-1)
    if(sum(picked_sick)>0){
      #add sick birds with a certain probability
      add_birds <- t(sapply(1:nrow(survey.live.detectables),function(i){mapply(function(s,p){rbinom(1,s,p)},survey.live.detectables[i,], c(0,picked_sick))}))
      survey_live_detectables <- survey_live_detectables+add_birds[,survey.positive.live.vars]
      
    }
    
    #get detection by survey flags
    detection.flags.survey <-apply(survey_live_detectables, 1, function(x) {
      sum(mapply(function(d, s) {
        rbinom(n = 1, size = d, p =  s)
      }, d = x, s =  se.survey.live[names(survey_live_detectables)%in% survey.positive.live.vars]))
    })
    #
    valid.times.survey <- survey.times[detection.flags.survey >= 1]
    first.ac.detection.survey.time <- if (length(valid.times.survey) > 0) min(valid.times.survey) else Inf
    
  }else first.ac.detection.survey.time <- Inf
  
  #
  if (is.infinite(first.ac.detection.survey.time )) {
    first.ac.detection.survey.time  = Inf; confirm.survey.result = "NA"; 
    live.confirmation.survey.result = "NA"; total.detectables.survey = 0; live.positives.survey  = NA;
  }else{
  while(length(valid.times.survey)>=1){
  #Active dead confirmation
  #confirmation test after a lag period to account for test delays: all dead birds
  confirm.test.day <- first.ac.detection.survey.time  + tlag
  confirm.dead <- confirm.dead.data %>%
    filter(day > first.ac.detection.survey.time  & day <= confirm.test.day) %>%
    summarize(total.detectables.survey = sum(daily.death.incidence.detectables, na.rm = TRUE))
  
  total.detectables.survey <- confirm.dead$total.detectables.survey %||% 0
  confirm.positive.survey <- rbinom(n = 1, size = total.detectables.survey, p = panimal * se) >= 1
  confirm.survey.result <- ifelse(confirm.positive.survey, "POSITIVE", "NEGATIVE")
  
  #Active live confirmation
  #secondary confirmation (independent to dead confirmation): live bird sampling on test day
  day.before <- confirm.test.day - 1
  live.samples <- confirm.live.detectables %>%
    filter(day > day.before & day <= confirm.test.day)
  
  if (nrow(live.samples) == 0) {
    live.confirmation.survey.result <- "NEGATIVE"
    total.live.detectables.survey <- 0
  } else {
    total.live.detectables.survey <- rowSums(select(live.samples, where(is.numeric), -day), na.rm = TRUE)
    total.live.detectables.survey <- sum(total.live.detectables.survey)
    
    if (total.live.detectables.survey == 0) {
      live.confirmation.survey.result <- "NEGATIVE"
    } else {
      positive.birds <- sum(live.samples %>% select(all_of(confirm.positive.live.vars)) %>% unlist(), na.rm = TRUE)
      live.panimal <- positive.birds / total.live.detectables.survey
      sample.size <- min(confirm.live.sample.size, total.live.detectables.survey)
      live.positives.survey  <- rbinom(n = 1, size = sample.size, p = live.panimal * se.confirm.live)
      live.confirmation.survey.result <- ifelse(live.positives.survey  >= 1, "POSITIVE", "NEGATIVE")
    }
  }
  if(confirm.survey.result != "POSITIVE" & live.confirmation.survey.result != "POSITIVE" & length(valid.times.survey) >0)
  {
    valid.times.survey = valid.times.survey[-1];
    first.ac.detection.survey.time <- if (length(valid.times.survey) > 0) min(valid.times.survey) else Inf}else{valid.times.survey = c()} 
  }
  }
  
  
  return(list(
    #bucket sampling results
    first.ac.detection.time = first.ac.detection.time,
    confirm.result = confirm.result,
    live.confirmation.result = live.confirmation.result,
    total.detectables = total.detectables,
    live.positives = ifelse(exists("live.positives"), live.positives, 0),
    #live bird survey results
    first.ac.detection.survey.time = first.ac.detection.survey.time,
    confirm.survey.result = confirm.survey.result,
    live.confirmation.survey.result = live.confirmation.survey.result,
    total.detectables.survey = total.detectables.survey,
    live.positives.survey = ifelse(exists("live.positives.survey"), live.positives.survey, 0)
    
    
  ))
}



#############################
## 3. Data Arrangement   ####
#############################
arrange.data.detection <- function(output, confirm.live.detectables.vars,survey.live.detectables.vars, time.interval.ac,time.interval.survey, init.ac, init.survey, dispose.interval = time.interval.ac) {
  #Passive detection dataset
  #For ratio (using daily.death.incidence):
  ratio.data <- output %>%
  arrange(day) %>%
  select(day, 
         deaths = N.dead, 
         daily.death.incidence) 
  
#For threshold (using daily.mortality.rate):
threshold.data <- output %>%
  arrange(day) %>%
  select(day, 
         deaths = N.dead,
         daily.death.incidence,
         mortality.rate = daily.mortality.rate)
  
  # Active detection dataset (using daily.death.incidence.detectables) and sampling interval assuming birds are disposed in between without testing:
  active.data <- output %>%
    arrange(day) %>%
    mutate(active.day = init.ac + ceiling((day - init.ac) / time.interval.ac) * time.interval.ac) %>% 
    mutate(non_disposed_day = (active.day - day)< dispose.interval)%>%filter(non_disposed_day ==TRUE)%>%
    group_by(active.day) %>%
    summarize(
      detectables = sum(daily.death.incidence.detectables, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(day = active.day)
  
  #Dead confirmation dataset (using daily.death.incidence.detectables): 
  confirm.dead.data <- output %>%
    arrange(day) %>%
    select(day, daily.death.incidence.detectables)
  
  # Live confirmation dataset:
  survey.live.data <- output %>%
    arrange(day) %>%
    filter(day %in% unique(init.survey + ceiling((day - init.survey) / time.interval.survey) * time.interval.survey))%>% #filter out sampling days
    select(day, all_of(survey.live.detectables.vars))
  
  # Live confirmation dataset:
  confirm.live.data <- output %>%
    arrange(day) %>%
    select(day, all_of(confirm.live.detectables.vars))
  
  return(list(
    ratio.data = ratio.data,
    threshold.data = threshold.data,
    active.data = active.data,
    survey.live.data = survey.live.data,
    confirm.live.data = confirm.live.data,
    confirm.dead.data = confirm.dead.data
  ))
}


###############################
## 4. Repeat Detection Module #
###############################

# function to run the detection process (active and passive) for each simulation run and repetition.
repeat.detection.time.surveillance <- function(output, reps, deaths.vars,
                                               detectables.vars, survey.live.detectables.vars, survey.positive.live.vars,confirm.live.detectables.vars, confirm.positive.live.vars,
                                               se, se.survey.live, se.confirm.live,picked_sick,
                                               time.interval.ac, dispose.interval = time.interval.ac, time.interval.survey,
                                               #init.ac,
                                               pfarm, panimal, tlag, survey.live.sample.size,confirm.live.sample.size,
                                               threshold, ints, intro.time, farm.size, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  
  #create a data frame from simulation output
  output.df <- output
  
  results.list <- list()
  
  for (i in seq_len(max(output.df$run))) {
    for (j in seq_len(reps)) {
      #generate separate random starting times:
      current.init.ac <- if(is.finite(time.interval.ac)){intro.time + round(runif(1, 0, time.interval.ac-1))} else { 10^6}
      current.init.survey <- intro.time +round(runif(1, 0, time.interval.survey-1))
      current.init.survey <- min(current.init.ac,current.init.survey)
      current.init.pas <- intro.time + round(runif(1, 0, 0)) #currently uses intro.time
      
      #arrange simulation data using arrange data functions
      det.data <- arrange.data.detection(
        output = output.df %>% filter(run == i),
        confirm.live.detectables.vars = confirm.live.detectables.vars,
        survey.live.detectables.vars = survey.live.detectables.vars,
        time.interval.ac = time.interval.ac,
        time.interval.survey = time.interval.survey,
        init.ac = current.init.ac,
        init.survey = current.init.survey
      )
      
      #Passive detection 
      pas.det.time.ratio <- passive.detection.ratio(
        ratio.data = det.data$ratio.data,
        init.pas = current.init.pas
      )
      pas.det.time.threshold <- passive.detection.time.threshold(
        threshold.data = det.data$threshold.data,
        threshold = threshold,
        ints = ints,
        init.pas = current.init.pas
      )
      
      pas.det.time <- min(pas.det.time.ratio, pas.det.time.threshold)
      
      # Active detection
      ac.det <- active.detection.time(
        times = det.data$active.data$day,
        survey.times = det.data$survey.live.data$day,
        lumped.detectables = det.data$active.data %>% select(-day),
        survey.live.detectables = det.data$survey.live.data[,c("day",survey.live.detectables.vars)],
        survey.positive.live.vars = survey.positive.live.vars,
        confirm.live.detectables = det.data$confirm.live.data[,c("day",confirm.live.detectables.vars)],
        confirm.dead.data = det.data$confirm.dead.data, 
        confirm.positive.live.vars = confirm.positive.live.vars,
        se = se,
        se.confirm.live = se.confirm.live,
        se.survey.live = se.survey.live,
        pfarm = pfarm,
        panimal = panimal,
        tlag = tlag,
        confirm.live.sample.size = confirm.live.sample.size,
        survey.live.sample.size = survey.live.sample.size,
        init.ac = current.init.ac,
        roundTime = 1,
        picked_sick = picked_sick 
      )
      #collect results. Note that the active surveillance on dead and live birds is taken together.
      result.row <- data.frame(
        run = i,
        rep = j,
        pas.det.time.ratio = pas.det.time.ratio,
        pas.det.time.threshold = pas.det.time.threshold,
        pas.det.time = pas.det.time,
        #bucket
        first.ac.detection.time = ac.det$first.ac.detection.time,
        confirm.result = ac.det$confirm.result,
        live.confirmation.result = ac.det$live.confirmation.result,
        total.detectables = ac.det$total.detectables,
        live.positives = ac.det$live.positives,
        #live bird survey
        first.ac.detection.survey.time = ac.det$first.ac.detection.survey.time,
        confirm.survey.result = ac.det$confirm.survey.result,
        live.confirmation.survey.result = ac.det$live.confirmation.survey.result,
        total.detectables.survey = ac.det$total.detectables.survey,
        live.positives.survey = ac.det$live.positives.survey,
        ac.success = pas.det.time > ac.det$first.ac.detection.time,
        ac.success.survey = pas.det.time > ac.det$first.ac.detection.survey.time
      )
      if(nrow(result.row)>1)stop("MOre than one entry!")
      results.list[[length(results.list) + 1]] <- result.row
    }
  }
  
  surveillance.time <- do.call(rbind, results.list)
  return(surveillance.time)
}

#############################
## 5. Scenario Processing  ##
#############################

# calculate the final outbreak size for each simulation run
calculate.final.size <- function(sim.out) {
  bind_rows(sim.out) %>%
    mutate(
      R = R.1 + R.2,
      DI = DI.1 + DI.2,
      DR = DR.1 + DR.2,
      DL = DL.1 + DL.2,
      final.size = R + DI + DR + DL
    ) %>%filter(!is.na(final.size)) %>%
    group_by(run) %>%
    summarize(final.size = max(final.size), .groups = "drop")
}

#Contaminated eggs at detection:
#get the detectables, dead, and shipped/on farm eggs at time of detection 
calculate.detectables.at.min.time <- function(sim.out, min.det.times) {
  sim.df <- sim.out %>%
    distinct(run, day, .keep_all = TRUE)
  
  #split detection times into finite and Inf rows
  finite.dets <- min.det.times %>% filter(is.finite(min.det.time))
  inf.dets <- min.det.times %>% filter(!is.finite(min.det.time))
  
  #join finite detection times with matching day 
  finite.df <- left_join(finite.dets, sim.df, by = c("run", "absolute.min.det.time" = "day")) %>%
    mutate(
      detect.N.dead = N.dead,
      detect.N.dead.detectables = N.dead.detectables,
      detect.N.live.detectables = N.live.detectables,
      detect.N.live.positives.survey = live.positives.survey,
      detect.N = N,
      detect.Egg_I_farm = Egg_I_farm,
      detect.Egg_I_shipped = Egg_I_shipped
    )
  
  #for Inf detection times: use the last (max day) value per run ie. the total number of shipped eggs since never detected 
  inf.df <- sim.df %>%
    group_by(run) %>%
    filter(day == max(day)) %>%
    ungroup() %>%
    right_join(inf.dets, by = "run") %>%
    mutate(
      detect.N.dead = N.dead,
      detect.N.dead.detectables = N.dead.detectables,
      detect.N.live.detectables = N.live.detectables,
      detect.N.live.positives.survey = live.positives.survey,
      detect.N = N,
      detect.Egg_I_farm = Egg_I_farm,
      detect.Egg_I_shipped = Egg_I_shipped
    )
  
  #combine
  df.out <- bind_rows(finite.df, inf.df) %>%
    select(run, absolute.min.det.time, time.interval.ac,
           detect.N.dead, detect.N.dead.detectables, detect.N.live.detectables,
           detect.Egg_I_farm, detect.Egg_I_shipped) %>%
    distinct()
  
  return(df.out)
}


# run the detection module for one scenario and a given active sampling frequency
run.detection.module <- function(sim.out, freq.sampling,freq.survey, dispose.interval, intro.time, farm.size, scenario.name, fs,...) {
  detection.result <- repeat.detection.time.surveillance(
    output = sim.out,
    reps = 1,
    #Passive
    deaths.vars = c("DI.1", "DI.2", "DR.1", "DR.2", "DS.1", "DS.2", "DL.1", "DL.2"),
    detectables.vars = c("DI.1", "DI.2", "DR.1", "DR.2", "DL.1", "DL.2"),
    survey.live.detectables.vars = c("S.1","S.2","L.1","L.2","I.1","I.2","R.1","R.2"),
    survey.positive.live.vars = c("I.1","I.2","R.1","R.2"),
    confirm.live.detectables.vars = c("I.1" ,"I.2"), #Only consider I birds as diseased
    confirm.positive.live.vars =  c("I.1", "I.2"), #I.1 excluded as sampling of healthy birds. If protocols change, only need to add in here.
    se = se_input,
    se.survey.live = se.survey.live_input, #in case different sensitivity for live bird test
    se.confirm.live = se.confirm.live_input, #in case different sensitivity for live bird test
    picked_sick = picked_sick_input,
    threshold = threshold_input, #detection threshold (don't need to multiply by farm size as uses daily.mortality.rates that already account for farm size)
    #Active
    time.interval.ac = freq.sampling,
    time.interval.survey = freq.survey,
    dispose.interval = dispose.interval,
    #init.ac = 0,
    pfarm = pfarm_input, #probability that this farm is being monitored
    panimal = panimal_input, #probability per detectables animal to being monitored,
    tlag = tlag_input, #time to confirmation
    roundTime = 1, #rounds to nearest day 
    survey.live.sample.size = survey.live.sample.size_input,
    confirm.live.sample.size = confirm.live.sample.size_input,
    ints = ints_input, #consecutive days for threshold 
    intro.time = intro.time,
    farm.size = farm.size,
    ...
  )
  
  detection.result %>%
    mutate(
      time.interval.ac = freq.sampling,
      time.interval.survey = freq.survey,
      scenario = scenario.name,
      intro.time = intro.time,
      N0 = farm.size,
      det.since.intro.pas = pas.det.time - intro.time,
      det.since.intro.ac = first.ac.detection.time - intro.time,
      min.det.time = pmin(det.since.intro.pas, det.since.intro.ac),
      absolute.min.det.time = round(min.det.time + intro.time)
    ) %>%
    left_join(fs, by = "run") %>%
    mutate(outbreak.size = ifelse(final.size < 0.1 * N0, "MINOR", "MAJOR"))
}



#############################
## 6. Main Execution     ####
#############################

output.dir <- "./HPC_output/DailyData/"
#change the pattern below to load different scenarios.
rds.files <- list.files(output.dir, pattern = "DailyData.*\\.RDS$", full.names = TRUE)
#rds.files <- list.files(output.dir, pattern = "Daily\\Data.*S1.*\\.RDS$", full.names = TRUE)
#rds.files <- list.files(output.dir, pattern = "Daily\\Data.*S2.1_HomoWaning.*\\.RDS$", full.names = TRUE)
#rds.files <- list.files(output.dir, pattern = "Daily\\Data.*S2.2_HeteroWaning.*\\.RDS$", full.names = TRUE)

# Bucket sampling as active surveillance ####
# 
# sampling.frequency <- c(2, 7, 14, 30)
# dispose.interval <- 1
# # process each simulation file (scenario)
# detection.results <- map(rds.files, function(file) {
#   output <- readRDS(file)
#   sim.out <- output$out
#   intro.time <- output$pars$intro.time
#   farm.size <- output$pars$N0
#   scenario.name <- output$pars$scenario
#   
#   message("Processing scenario: ", scenario.name)
#   
#   fs <- calculate.final.size(sim.out)
#   
#   # run detection for each active sampling frequency and store results
#   results <- map_dfr(sampling.frequency, function(freq) {
#     message("Active Detection Interval: ", freq)
#     run.detection.module(sim.out, freq, dispose.interval,intro.time, farm.size, scenario.name, fs)
#   })
#   detect.at.min <- calculate.detectables.at.min.time(sim.out, results)
#   
#   #join and return
#   left_join(results, detect.at.min, by = c("run", "absolute.min.det.time", "time.interval.ac"))
#   
# })
# 
# # Combine all results
# combined.results <- bind_rows(detection.results)%>%
#   select(-total.detectables, -live.positives)
# 
# #change pattern for correct scenario name
# #create directory HPC_output/detection_results if it does not exist
# if (!dir.exists("./HPC_output/Detection/")) {
#   dir.create("./HPC_output/Detection/")
# }
# #save results
# saveRDS(combined.results, file = paste0("HPC_output/Detection/", format(Sys.Date(), "%Y_%m_%d_"), "Detection_results_complete.RDS"))

# EU protocol as active surveillance - the minimum detection is based on passive and bucket sampling ####
sampling.frequency <- c(2,7,14,30,Inf)
survey.frequency <- c(30)
dispose.interval <- 1
# process each simulation file (scenario)

detection.results <- map(rds.files, function(file) {
  output <- readRDS(file)
  sim.out <- output$out
  intro.time <- output$pars$intro.time
  farm.size <- output$pars$N0
  scenario.name <- output$pars$scenario
  
  message("Processing scenario: ", scenario.name)
  
  fs <- calculate.final.size(sim.out)
  
  # run detection for each active sampling frequency and store results
  results <- map2_dfr(sampling.frequency,survey.frequency, function(freq.sampling, freq.survey) {
    message("Active Detection Interval: ", freq.sampling," And live birds surveillance: ", freq.survey)
    run.detection.module(sim.out, freq.sampling, freq.survey, dispose.interval,intro.time, farm.size, scenario.name, fs, seed = 1440)
  })
  detect.at.min <- calculate.detectables.at.min.time(sim.out, results)
  
  #join and return
  left_join(results, detect.at.min, by = c("run", "absolute.min.det.time", "time.interval.ac"))
  
})

# Combine all results
combined.results <- bind_rows(detection.results)%>%
  select(-total.detectables, -live.positives)

#change pattern for correct scenario name
#create directory HPC_output/detection_results if it does not exist
if (!dir.exists("./HPC_output/Detection/")) {
  dir.create("./HPC_output/Detection/")
}
#save results
saveRDS(combined.results, file = paste0("HPC_output/Detection/", format(Sys.Date(), "%Y_%m_%d_"), "Detection_results_complete_EU.RDS"))

