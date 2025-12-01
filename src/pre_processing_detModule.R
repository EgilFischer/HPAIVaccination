#########################################################
#                                                        
#                 Data pre-processing                            
#                   
#########################################################

#Input: saved simulation/scenario output with multiple time steps per day
#Output: simulation/scenarios with aggregated daily data, and new columns for daily.death.incidence (total and for detectables), daily.mortality.rate, etc, which are required by "detectionModule.R" 

source("./src/loadLibraries.R") 

output.dir <- "./HPC_output/"# location to find simulations
#change pattern as needed
file.list <- list.files(output.dir, pattern = "*\\.RDS$", full.names = TRUE)

for (f in file.list) {
  sim <- readRDS(f)
  
  #extract simulation output as a data frame
  df <- as.data.frame(sim$out)
  
  #aggregate to daily data by taking the first "snapshot" of each day per run
  df.daily <- df %>%
    mutate(day = floor(time)) %>%
    group_by(run, day) %>%
    summarize(across(everything(), first), .groups = "drop") %>%
    arrange(run, day)
  
  #add new columns on the daily aggregated data (so easier to check the detection results later, and make sure pop size is stable):
  df.daily <- df.daily %>%
    mutate(
      N.dead = DS.1 + DS.2 + DL.1 + DL.2 + DI.1 + DI.2 + DR.1 + DR.2,
      N.dead.detectables = DL.1 + DL.2 + DI.1 + DI.2 + DR.1 + DR.2,
      N.live.detectables = L.1 + L.2 + I.1 + I.2 + R.1 + R.2,
      N.total = N + N.dead
    )
  
  #add new columns for daily death incidence and mortality rate
  df.daily <- df.daily %>%
    arrange(run, day) %>%
    group_by(run) %>%
    mutate(
      daily.death.incidence = c(N.dead[1], diff(N.dead)),
      daily.mortality.rate = daily.death.incidence / N, 
      daily.death.incidence.detectables = c(N.dead.detectables[1], diff(N.dead.detectables))
    ) %>%
    ungroup()
  
  #reorder columns for nicer reading
  df.daily <- df.daily %>%
    select(
      run, day, time, 
      S.1, S.2, L.1, L.2, I.1, I.2, R.1, R.2,
      DS.1, DS.2, DL.1, DL.2, DI.1, DI.2, DR.1, DR.2,
      N, N.dead, N.total, 
      daily.death.incidence, daily.death.incidence.detectables, daily.mortality.rate,
      N.dead.detectables, N.live.detectables,
      #remaining columns
      everything()
    )
  
  #new sim output:'pars' unchanged and 'out' replaced by the daily data
  sim.daily <- list(out = df.daily, pars = sim$pars)
  
  #save the daily pre-aggregated data with new name
  #create a subdirectory "DailyData" if it does not exist
  if (!dir.exists(file.path(dirname(f), "DailyData"))) {
    dir.create(file.path(dirname(f), "DailyData"))
  }
  #create a new file name in subdirectory "DailyData" under the original file's directory
  new.file <- file.path(dirname(f), paste0("/DailyData/DailyData.", basename(f)))
  saveRDS(sim.daily, new.file)
  message("Created daily aggregate: ", new.file)
}

#Optional: 
# #view the daily data file (change file name as necessary)
# output.dir <- "./HPC_output/"
# file.list <- list.files(output.dir, pattern = "DailyData.*\\.RDS$", full.names = TRUE)
# choose_file <- 4; file.list[choose_file]
# S.sim <- readRDS(file.list[choose_file])
# S.out <- S.sim$out
# View(as.data.frame(S.out))
# #check runs by visualizing daily data
# ggplot(data =S.out, aes(x = day, y = daily.death.incidence, group = as.factor(run))) +
#   geom_line() +
#   labs(title = "Daily Death Incidence by Run", x = "Day", y = "Daily Death Incidence") +
#   theme_minimal() +
#   scale_color_discrete(name = "Run")
# 
# head(S.out)
