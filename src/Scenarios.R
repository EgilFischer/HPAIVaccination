#########################################################
#                                                        
#                  Simulate scenarios: tleap                              
#                                                        
#                  Author:Egil Fischer                              
#                  Contact:           e.a.j.fischer@uu.nl                   
#                  Creation date: 2-5-2023                         
#########################################################

#load libraries
source("./src/loadLibraries.R") 
source("./src/faster_multitypetransitionsSEIR_tleap.R")
if(!dir.exists("./output")){dir.create("./output")}

#NOTES: 
#Input: faster_multitypetransitionsSEIR_tleap.R" sourced above
#Output: simulation/scenario output

# SCENARIOS #
  #Assumption: 100% of flock is protected at start
#Scenario 0: Baseline (no vaccination)
#Scenario 1: Flock size & introduction time
  #introduction time when vx coverage is 100%, 90%, 75%, etc 
#Scenario 2.1: Waning immunity with same strain & intro time 
#Scenario 2.2: Waning immunity with different strain & intro time


############################
##  1. Run scenarios    ####
############################

#Layers####
#Baseline parameters for layer flock from faster_multitypetransitionsSEIR_tleap.R ####
#Define number of types
itypes = 2; #Type 1  = not protected by vaccination; Type 2 = protected by vaccination

#List of all parameters
param.list.baseline.layer <- list(
  #Simulation settings####
  scenario = "Default params", 
  runs = 100, #number of runs
  deltat = 0.005, #time step; nb: 0.01 is faster but not accurate, 0.001 is very slow and only mild improvement in accuracy
  #Flock settings####
  length.round = 19*30,#length of the round. Simulation starts at time 0
  age_at_d0 = 120, #unit: days; model starts at t = 0 when chickens are 120 days/17 wks old (important for egg production curve and vx buildup)
  N0 = 46000, #population size
  #Vaccination settings ####
  itypes = itypes, #types
  p.protect = 1, #proportion of the flock that is vaccinated
  #Introduction settings ####
  no.introduction= 1, #initially infected at introduction
  intro.time = 0, #time during the round at which the infection is introduced
  #Infection settings ####
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
  mortRate = 0.0002, #per capita death rate for layers (Hobbelen et al. 2020; DOI: https://doi.org/10.1038/s41598-020-68623-w)
  #Egg settings ####
  eh = 0.57, # daily egg-laying rate for healthy chickens (S, E, R)
  ei = 0.57 * 0.5,# daily egg-laying rate for infected chickens (I) (50% reduction in egg laying)
  disfigured = 0.1, # daily removal rate for infected eggs (10% do not pass inspection)
  pickup_time = 3 #how often eggs are picked up 
)

### Scenario 0: Baseline Size and introduction time ###
scenario.list.size.baseline <- list()
for(i in c(1:3)){
  for(j in c(1:6)){
    param.list <- param.list.baseline.layer
    param.list$deltat <- 0.005
    param.list$runs <- 100
    param.list$p.protect <- 0
    param.list$N0 <-  c(15000,46000,61000)[i] # from farm size data (q1,q2,q3)
    param.list$intro.time<- c(0, 259, 283.5, 311, 371.5, 504)[j]
    param.list$scenario <- paste0("Baseline_",param.list$N0,"_IntroTime_",param.list$intro.time)
    scenario.list.size.baseline[[6*(i-1)+j]]<- param.list
  }
}

for(i in  c(1:length(scenario.list.size.baseline))){
  print(scenario.list.size.baseline[[i]]$scenario);
  
  sim.out<- sim.multitypeSEIR_tleap(scenario.list.size.baseline[[i]], seed = 1440)
  #record output and parameters, save as RDS
  op <- list(out = sim.out, pars = scenario.list.size.baseline[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y_%m_%d_"),scenario.list.size.baseline[[i]]$scenario,".RDS"))
}

### Scenario 1: Size and introduction time ###
scenario.list.size.introtime <- list()
for(i in c(1:3)){
  for(j in c(1:6)){
    param.list <- param.list.baseline.layer
    param.list$deltat <- 0.005
    param.list$runs <- 100
    param.list$N0 <-  c(15000,46000,61000)[i] # from farm size data (q1,q2,q3)
    param.list$intro.time<- c(0, 259, 283.5, 311, 371.5, 504)[j]
    param.list$scenario <- paste0("S1_LayerSize_",param.list$N0,"_IntroTime_",param.list$intro.time)
    scenario.list.size.introtime[[6*(i-1)+j]]<- param.list
  }
}

for(i in  c(1:length(scenario.list.size.introtime))){
  print(scenario.list.size.introtime[[i]]$scenario);

  sim.out<- sim.multitypeSEIR_tleap(scenario.list.size.introtime[[i]], seed = 1440)
  #record output and parameters, save as RDS
  op <- list(out = sim.out, pars = scenario.list.size.introtime[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y_%m_%d_"),scenario.list.size.introtime[[i]]$scenario,".RDS"))
}

### Scenario 2.1:  Waning same strain ###
scenario.list.homo.waning.introtime <- list()
for(i in c(1:6)){
  param.list <- param.list.baseline.layer
  param.list$intro.time<- c(0, 259, 283.5, 311, 371.5, 504)[i]
  param.list$scenario <- paste0("S2.1_HomoWaning_IntroTime_",param.list$intro.time,"_LayerSize_46000")
  scenario.list.homo.waning.introtime[[i]]<- param.list
}

for(i in  c(1:length(scenario.list.homo.waning.introtime))){
  print(scenario.list.homo.waning.introtime[[i]]$scenario);

sim.out<- sim.multitypeSEIR_tleap(scenario.list.homo.waning.introtime[[i]], seed = 1440)
  #record output and parameters, save as RDS
  op <- list(out = sim.out, pars = scenario.list.homo.waning.introtime[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y_%m_%d_"),scenario.list.homo.waning.introtime[[i]]$scenario,".RDS"))
}

### Scenario 2.2:  Waning different strain ###
scenario.list.waning.hetero.introtime <- list()
for(i in c(1:6)){
  param.list <- param.list.baseline.layer;
  param.list$trans.mean.wane = 280; #only considers transition from type 2 to 1
  param.list$trans.var.wane = 140^2; #only considers transition from type 2 to 1
  param.list$intro.time<- c(0, 259, 283.5, 311, 371.5, 504)[i]
  param.list$scenario <- paste0("S2.2_HeteroWaning_IntroTime_",param.list$intro.time,"_LayerSize_46000");
  scenario.list.waning.hetero.introtime[[i]]<- param.list
}


for(i in  c(1:length(scenario.list.waning.hetero.introtime))){
  print(scenario.list.waning.hetero.introtime[[i]]$scenario);
  
sim.out<- sim.multitypeSEIR_tleap(scenario.list.waning.hetero.introtime[[i]], seed = 1440)
  
  #record output and parameters, save as RDS
  op <- list(out = sim.out, pars = scenario.list.waning.hetero.introtime[[i]])
  saveRDS(op, file = paste0("output/",format(Sys.Date(),"%Y_%m_%d_"), scenario.list.waning.hetero.introtime[[i]]$scenario,".RDS"))
}

#Optional:
############################
##  2. View a scenario   ###
############################

#change file name as necessary
output.dir <- "./output/"
S.results <- file.path(output.dir, "2025_03_20_S1_LayerSize_61000_IntroTime_0.RDS")
S.sim <- readRDS(S.results)
S.out <- S.sim$out
View(as.data.frame(S.out))


############################
##  3. Plot scenarios   ####
############################

#make folder for saving figures
if(!dir.exists("./figures")){dir.create("./figures")}


#load scenario data as a combined data frame
load_scenario_data <- function(pattern, output_dir = "./output/") {
  rds_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  scenario_data <- lapply(rds_files, function(file) {
    data <- readRDS(file)
    #convert simulation output to a data frame
    df <- as.data.frame(data$out)
    #add scenario name and run identifier
    df$scenario <- data$pars$scenario
    df$run <- data$out$run
    #aggregate compartments
    df <- df %>%
      mutate(
        S = S.1 + S.2,
        L = L.1 + L.2,
        I = I.1 + I.2,
        R = R.1 + R.2
      ) %>%
      select(time, run, scenario, S, L, I, R)
    return(df)
  })
  combined_data <- do.call(rbind, scenario_data)
  return(combined_data)
}

#output directory
output_dir <- "./output/"

#load data for each scenario as data frames
scenario1_df   <- load_scenario_data(pattern = "S1_LayerSize_61000_IntroTime_505.*\\.RDS$", output_dir = output_dir)
scenario2_1_df <- load_scenario_data(pattern = "S2\\.1_HomoWaning_.*\\.RDS$", output_dir = output_dir)
scenario2_2_df <- load_scenario_data(pattern = "S2\\.2_HeteroWaning_.*\\.RDS$", output_dir = output_dir)

### Plotting Scenario 1 ####
scenario1_long <- scenario1_df %>%
  pivot_longer(
    cols = c(S, L, I, R),
    names_to = "variable",
    values_to = "value"
  )

scenario1_plot <- ggplot(scenario1_long, aes(x = time, y = value, colour = variable, group = interaction(run, variable))) +
  geom_path() +
  facet_wrap(~scenario, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(
    title = "SEIR Chickens Across Scenarios",
    subtitle = "Scenario 1: Layer Size and Introduction Time",
    x = "Time (days)",
    y = "Number of Chickens",
    colour = "Compartment"
  )

ggsave(filename = "./figures/scenario1_plot.png", plot = scenario1_plot, width = 10, height = 8, dpi = 300)


### Plotting Scenario 2.1 ####
scenario2_1_long <- scenario2_1_df %>%
  pivot_longer(
    cols = c(S, L, I, R),
    names_to = "variable",
    values_to = "value"
  )

scenario2_1_plot <- ggplot(scenario2_1_long, aes(x = time, y = value, colour = variable, group = interaction(run, variable))) +
  geom_path() +
  facet_wrap(~scenario, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(
    title = "SEIR Chickens Across Scenarios",
    subtitle = "Scenario 2.1: Homo Waning and Introduction Time",
    x = "Time (days)",
    y = "Number of Chickens",
    colour = "Compartment"
  )

ggsave(filename = "./figures/scenario2.1_plot.png", plot = scenario2_1_plot, width = 10, height = 8, dpi = 300)


### Plotting Scenario 2.2 ####
scenario2_2_long <- scenario2_2_df %>%
  pivot_longer(
    cols = c(S, L, I, R),
    names_to = "variable",
    values_to = "value"
  )

scenario2_2_plot <- ggplot(scenario2_2_long, aes(x = time, y = value, colour = variable, group = interaction(run, variable))) +
  geom_path() +
  facet_wrap(~scenario, scales = "free", ncol = 3) +
  theme_minimal() +
  labs(
    title = "SEIR Chickens Across Scenarios",
    subtitle = "Scenario 2.2: Hetero Waning and Introduction Time",
    x = "Time (days)",
    y = "Number of Chickens",
    colour = "Compartment"
  )

ggsave(filename = "./figures/scenario2.2_plot.png", plot = scenario2_2_plot, width = 10, height = 8, dpi = 300)
  



############################
##  4. Plot final size   ###
############################

#load and process .RDS files for all runs ####
load_and_process_final_size <- function(pattern, output_dir = "./output/") {
  rds_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  
  scenario_data <- lapply(rds_files, function(file) {
    data <- readRDS(file)
    df <- as.data.frame(data$out)         
    df$scenario <- data$pars$scenario       
    df$run <- data$out$run                  
    
    #calculate final size 
    df <- df %>%
      mutate(
        final_size = (R.1 + R.2) + (DI.1 + DI.2) + (DR.1 + DR.2)
      ) %>%
      group_by(run, scenario) %>%
      summarise(fs = max(final_size, na.rm = TRUE), .groups = "drop")
    
    return(df)
  })
  
  combined_data <- do.call(rbind, scenario_data)
  return(combined_data)
}

#extract parameters from one .RDS file for a scenario ####
extract_scenario_parameters <- function(pattern, output_dir = "./output/") {
  rds_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(rds_files) == 0) stop("No files found for the pattern!")
  
  #extract parameters
  param_list <- readRDS(rds_files[1])$pars
  return(param_list)
}

#extract parameters for each scenario ####
scenario1_params <- extract_scenario_parameters(pattern = "S1_LayerSize_.*\\.RDS$")
scenario2_1_params <- extract_scenario_parameters(pattern = "S2.1_HomoWaning_.*\\.RDS$")
scenario2_2_params <- extract_scenario_parameters(pattern = "S2.2_HeteroWaning_.*\\.RDS$")

#load and process data for each scenario ####
scenario1_fs <- load_and_process_final_size(pattern = "S1_LayerSize_.*\\.RDS$")
scenario2_1_fs <- load_and_process_final_size(pattern = "S2.1_HomoWaning_.*\\.RDS$")
scenario2_2_fs <- load_and_process_final_size(pattern = "S2.2_HeteroWaning_.*\\.RDS$")

#plotting function for histograms ####
plot_tile_histograms <- function(data, param_list, title, binwidth = 2000) {
  #calculate small outbreaks (<10% of population) per scenario
  small_outbreaks_df <- data %>%
    group_by(scenario) %>%
    summarise(
      small_outbreaks = sum(fs < 0.1 * param_list$N0),
      total_runs = n(),
      .groups = "drop"
    )
  
  #make subtitle for each scenario; 
  scenario_details <- small_outbreaks_df %>%
    mutate(
      subtitle = paste0("# small outbreaks: ", small_outbreaks, 
                        " out of ", total_runs, " runs")
    )
  custom_labeller <- as_labeller(setNames(scenario_details$subtitle, scenario_details$scenario))
  
  #make histogram plot
  p <- ggplot(data, aes(x = fs)) +
    geom_histogram(binwidth = binwidth, fill = "salmon", color = "black", boundary = 0) +
    facet_wrap(~scenario, scales = "free_y", ncol = 3, labeller = custom_labeller) +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "Outbreak Size",
      y = "Frequency"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  return(p)
}

# Plot Scenario 1 ####
scenario1_histo <- plot_tile_histograms(
  scenario1_fs,
  param_list = list(N0 = scenario1_params$N0, runs = scenario1_params$runs),
  title = "Distribution of Outbreak Sizes: Scenario 1"
)
ggsave(filename = "./figures/scenario1_histo.png", plot = scenario1_histo, width = 10, height = 8, dpi = 300)

# Plot Scenario 2.1 ####
scenario2.1_histo <- plot_tile_histograms(
  scenario2_1_fs,
  param_list = list(N0 = scenario2_1_params$N0, runs = scenario2_1_params$runs),
  title = "Distribution of Outbreak Sizes: Scenario 2.1"
)
ggsave(filename = "./figures/scenario2.1_histo.png", plot = scenario2.1_histo, width = 10, height = 8, dpi = 300)

# Plot Scenario 2.2 ####
scenario2.2_histo <- plot_tile_histograms(
  scenario2_2_fs,
  param_list = list(N0 = scenario2_2_params$N0, runs = scenario2_2_params$runs),
  title = "Distribution of Outbreak Sizes: Scenario 2.2"
)
ggsave(filename = "./figures/scenario2.2_histo.png", plot = scenario2.2_histo, width = 10, height = 8, dpi = 300)


#printing fs
cat("Final sizes for Scenario 1:\n")
print(scenario1_fs)

cat("\nFinal sizes for Scenario 2.1:\n")
print(scenario2_1_fs)

cat("\nFinal sizes for Scenario 2.2:\n")
print(scenario2_2_fs)



#########################
##  5. Plot deaths   ####
#########################

#load .RDS files and combine them into a data frame
load_scenario_data <- function(pattern, output_dir = "./output/") {
  rds_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  scenario_data <- lapply(rds_files, function(file) {
    data <- readRDS(file)
    df <- as.data.frame(data$out)           
    df$scenario <- data$pars$scenario         
    df$run <- data$out$run                    
    #select only the columns for fs calculation
    df <- df %>%
      select(time, run, scenario, DI.1, DI.2, DR.1, DR.2)
    return(df)
  })
  combined_data <- do.call(rbind, scenario_data)
  return(combined_data)
}

#load scenario data
scenario1_data <- load_scenario_data(pattern = "S1_LayerSize_61000_IntroTime_505.*\\.RDS$")
scenario2_1_data <- load_scenario_data(pattern = "S2.1_HomoWaning_.*\\.RDS$")
scenario2_2_data <- load_scenario_data(pattern = "S2.2_HeteroWaning_.*\\.RDS$")

#create the plot from a data frame
plot_scenario <- function(data, title, subtitle) {
  data_long <- data %>%
    pivot_longer(
      cols = c(DI.1, DI.2, DR.1, DR.2),
      names_to = "variable",
      values_to = "value"
    )
  
  p <- ggplot(data_long, aes(x = time, y = value, colour = variable, group = interaction(run, variable))) +
    geom_path() +
    facet_wrap(~scenario, scales = "free", ncol = 3) +
    theme_minimal() +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Time (days)",
      y = "Number of Chickens",
      colour = "Compartment"
    )
  
  return(p)
}

# Plot Scenario 1 
p1 <- plot_scenario(
  scenario1_data,
  title = "DI/DR Chickens Across Scenarios",
  subtitle = "Scenario 1: Layer Size and Introduction Time"
)
print(p1)

# Plot Scenario 2.1 
p2 <- plot_scenario(
  scenario2_1_data,
  title = "DI/DR Chickens Across Scenarios",
  subtitle = "Scenario 2.1: Homo Waning and Introduction Time"
)
print(p2)

# Plot Scenario 2.2 
p3 <- plot_scenario(
  scenario2_2_data,
  title = "DI/DR Chickens Across Scenarios",
  subtitle = "Scenario 2.2: Hetero Waning and Introduction Time"
)
print(p3)






