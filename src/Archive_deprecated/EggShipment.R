#Egg shipment 
source("./src/loadLibraries.R") 
library(dplyr)

output.dir <- "./output/"
#change the pattern below to load different scenarios.
sim_files <- list.files(output.dir, pattern = "S1_LayerSize.*\\.RDS$", full.names = TRUE)
detection_files <- list.files(output.dir, pattern = "Detection_results_complete_S1.*\\.RDS$", full.names = TRUE)

# User-defined parameters (should match your SEIR simulation settings)
pickup_time <- 3    # Egg pickup occurs every 3 days
deltat <- 0.005     # Time step used in simulation

# Function to process one scenario by processing each run separately
processScenario <- function(sim_file, det_file) {
  
  # Load simulation object and detection results
  sim_obj <- readRDS(sim_file)
  det_df <- readRDS(det_file)
  
  # Convert simulation output to a data frame.
  # We assume the simulation output is stored in an element named 'out'
  # with columns: time, run, Egg_H, Egg_I, etc.
  sim_df <- as.data.frame(sim_obj$out)
  
  # Pre-split simulation data by run to reduce repeated filtering on a large data frame
  sim_list <- split(sim_df, sim_df$run)
  
  # Initialize a list to store results for each run
  results_list <- vector("list", nrow(det_df))
  
  for(i in seq_len(nrow(det_df))) {
    run_id <- det_df$run[i]
    min_det_time <- det_df$min.det.time[i]
    
    # Extract simulation data for the current run from the pre-split list
    sim_run <- sim_list[[as.character(run_id)]]
    
    # Subset simulation data up to the detection time
    sim_run_subset <- sim_run %>% filter(time <= min_det_time)
    
    # ----- Calculate Holding Eggs at Detection Day -----
    # Select rows where the day (floor of time) equals min_det_time.
    # If multiple records exist, take the maximum values.
    holding_data <- sim_run_subset %>% filter(floor(time) == min_det_time)
    if(nrow(holding_data) > 0) {
      holding_Egg_H <- max(holding_data$Egg_H, na.rm = TRUE)
      holding_Egg_I <- max(holding_data$Egg_I, na.rm = TRUE)
    } else {
      holding_Egg_H <- NA
      holding_Egg_I <- NA
    }
    
    # ----- Reconstruct Shipped Eggs from Completed Pickup Cycles -----
    # Mark pickup events when (time %% pickup_time < deltat) and time > deltat.
    cycles_df <- sim_run_subset %>% 
      mutate(
        pickup_event = (time %% pickup_time < deltat & time > deltat),
        cycle = floor(time / pickup_time)
      )
    
    # Identify cycles in which at least one pickup event occurred
    complete_cycles <- cycles_df %>% 
      group_by(cycle) %>% 
      summarise(complete = any(pickup_event), .groups = "drop") %>% 
      filter(complete) %>% 
      pull(cycle)
    
    if (length(complete_cycles) > 0) {
      shipped_summary <- cycles_df %>% 
        filter(cycle %in% complete_cycles) %>%
        group_by(cycle) %>%
        summarise(shipped_Egg_H = max(Egg_H, na.rm = TRUE),
                  shipped_Egg_I = max(Egg_I, na.rm = TRUE),
                  .groups = "drop")
      
      shipped_totals <- shipped_summary %>% 
        summarise(shipped_Egg_H = sum(shipped_Egg_H),
                  shipped_Egg_I = sum(shipped_Egg_I))
      
      shipped_Egg_H <- shipped_totals$shipped_Egg_H
      shipped_Egg_I <- shipped_totals$shipped_Egg_I
    } else {
      shipped_Egg_H <- 0
      shipped_Egg_I <- 0
    }
    
    # Store the results for the current run
    results_list[[i]] <- data.frame(run = run_id,
                                    min.det.time = min_det_time,
                                    holding_Egg_H = holding_Egg_H,
                                    holding_Egg_I = holding_Egg_I,
                                    shipped_Egg_H = shipped_Egg_H,
                                    shipped_Egg_I = shipped_Egg_I,
                                    stringsAsFactors = FALSE)
    
    # Optionally, run garbage collection to free memory:
    gc()
  }
  
  results_final <- bind_rows(results_list) %>%
    mutate(simulation_file = basename(sim_file))
  
  return(results_final)
}

# Process all scenarios by pairing simulation and detection files.
# This assumes that sim_files and detection_files are in the same order (adjust if necessary).
all_results <- map2_dfr(sim_files, detection_files, processScenario)

# View the final results: for each run, you'll have min.det.time,
# holding eggs (healthy and unhealthy) on that day, and cumulative shipped eggs up to that time.
print(all_results)

