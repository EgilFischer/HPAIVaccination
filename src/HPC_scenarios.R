##########################################################################
# Create R-file for scenario's to run on HPC                             #
##########################################################################
# converts list to a string with name and value pairs
list_to_string <- function(my_list) {
  # Initialize a character vector to hold the string representations
  list_strings <- character(length(my_list))
  names(list_strings) <- names(my_list)
  
  # Convert each element to a string
  for (i in seq_along(my_list)) {
    element <- my_list[[i]]
    if (is.character(element)) {
      list_strings[i] <- paste0('"', element, '"')
    } else if (is.matrix(element)) {
      # Convert matrix to a string representation
      list_strings[i] <- paste0("matrix(c(", paste(as.vector(element), collapse = ", "), "), ncol = ", ncol(element), ")")
    } else if (length(element)>1) {
      # Recursively convert lists
      list_strings[i] <-  paste0("c(",paste0( element,  collapse = ", "),")")
    } else if (is.numeric(element) || is.logical(element)) {
      list_strings[i] <- as.character(element)
    } else if (is.null(element)) {
      stop("Unsupported NULL data type in the list.")
    } else {
      stop("Unsupported data type in the list.")}
  }
  
  # Combine into a single string
  result <- paste(
    sapply(names(list_strings), function(name) {
      paste0(name, "=", list_strings[name] )
    }),
    collapse = ",\n"
  )
  
  return (result)
}

#function producing R-scripts
generate_r_script <- function(param_list, scenario_seed = 12345) {
  # Define the content of the R script using the provided param_list
  script_content <- paste0(
    "# Scenario:", param_list$scenario,"\n",
    " param_list <- list(", list_to_string(param_list), ") \n",
    "# Source the simulator\n",
    "source(\"/home/uu_vet_te/efischer/HPAI_vaccination/src/faster_multitypetransitionsSEIR_tleap.R\")\n",
    "# Run the simulator with the parameter set\n",
    "sim.out <- sim.multitypeSEIR_tleap(param_list, seed = " ,scenario_seed,")\n",
    "op <- list(out = sim.out, pars = param_list)\n",
    'saveRDS(op, file = paste0("/home/uu_vet_te/efischer/HPAI_vaccination/output/", format(Sys.Date(), "%Y_%m_%d_"), param_list$scenario, ".RDS"))'
  )
  
  # Create a file name based on the scenario and current date
  file_name <- paste0("./HPC_scenarios/script_", format(Sys.Date(), "%Y_%m_%d_"), param_list$scenario, ".R")
  
  # Write the content to the file
  writeLines(script_content, con = file_name)
  
  # Return the name of the generated file
  return(file_name)
}


# baseline parameters
# Scenario: baseline layer farm ####
param_list_baseline_layer <- list(
  #Simulation settings####
  scenario="baseline_layer_farm", 
  runs = 100, #number of runs
  deltat = 0.005, #time step; nb: 0.01 is faster but not accurate, 0.001 is very slow and only mild improvement in accuracy
  #Flock settings####
  length.round = 19*30,#length of the round. Simulation starts at time 0
  age_at_d0 = 120, #unit: days; model starts at t = 0 when chickens are 120 days/17 wks old (important for egg production curve and vx buildup)
  N0 = 46000, #population size
  #Vaccination settings ####
  itypes = 2, #types
  p.protect = 0, #proportion of the flock that is vaccinated
  #Introduction settings ####
  no.introduction= c(1,0), #initially infected at introduction forcing initial infection to be of low-titre 
  intro.time = 0, #time during the round at which the infection is introduced
  #Infection settings ####
  beta = matrix(c(1.99, 1.99,0.23,0.23), ncol = 2), #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
  latency.period = c(1,1), #latency period of 1 day
  k.latency =2, #k parameter Erlang distribution
  infectious.period = c(4.28,3.59), #Duration infectious period 
  k.infectious = 2, #k parameter Erlang distribution, k=1 gives sim results comparable to pmaj predictions for exponential 
  trans.mean.wane = Inf, #only considers transition from type 2 to 1; for waning off set to 20000
  trans.var.wane = Inf, #only considers transition from type 2 to 1
  trans.mean.buildup = Inf, #mean time to build-up immunity
  trans.var.buildup = Inf, #variance of time to build-up immunity
  pdie = c(1.0,0.01),#probability of dying at end of infectious period
  mortRate = 2e-4, #per capita death rate for layers (Hobbelen et al. 2020; DOI: https://doi.org/10.1038/s41598-020-68623-w)
  #Egg settings ####
  eh = 0.57, # daily egg-laying rate for healthy chickens (S, E, R)
  ei = 0.285,# daily egg-laying rate for infected chickens (I) (50% reduction in egg laying)
  disfigured = 0.1, # daily removal rate for infected eggs (10% do not pass inspection)
  pickup_time = 3 #how often eggs are picked up 
)

# Generate scripts for simulations ####

generate_r_script(param_list_baseline_layer)

#Scenario baseline preproduction farm ####
param_list_baseline_preproduction <-param_list_baseline_layer
param_list_baseline_preproduction$scenario <-"baseline_preproduction"
param_list_baseline_preproduction$length.round <- 120 #length of the round. Simulation starts at time 0
param_list_baseline_preproduction$age_at_d0 <-0 #age at day 0 = 0
param_list_baseline_preproduction$trans.mean.buildup<-Inf
param_list_baseline_preproduction$trans.var.buildup <- Inf

generate_r_script(param_list_baseline_preproduction)

#fitted build up parameters ####
# Group      Mean  Variance
# 1      CEVA  7.266798  210.4728
# 2        BI 42.174651 1045.1162
# 3 BIBOOSTER 41.403477  844.1632

# Group Treatment          mean_protected
# < chr > < chr >                       < dbl >
#   1 1     CEVA                        0.961
# 2 2     BI                          0.923
# 3 3     BI + booster                  0.977
# 4 4     Neg. controle BI            0
# 5 4     Neg. controle CEVA          0.075
# > 

# Scenario: vaccination CEVA
param_list_CEVA_layer <- param_list_baseline_layer
param_list_CEVA_layer$scenario <- "CEVA_layer_farm"
param_list_CEVA_layer$trans.mean.buildup <- 7.266798
param_list_CEVA_layer$trans.var.buildup <- 210.4728
param_list_CEVA_layer$p.protect <- 0.961
param_list_CEVA_layer$beta = matrix(c(1.94, 1.94,0.25,0.25), ncol = 2) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_CEVA_layer$infectious.period = c(4.18,3.86) #Duration infectious period 
generate_r_script(param_list_CEVA_layer)

# Scenario: vaccination BI
param_list_BI_layer <- param_list_baseline_layer
param_list_BI_layer$scenario <- "BI_layer_farm"
param_list_BI_layer$trans.mean.buildup <- 42.174651
param_list_BI_layer$trans.var.buildup <- 1045.1162
param_list_BI_layer$p.protect <- 0.923
param_list_BI_layer$beta = matrix(c(2.13, 2.13,0.21,0.21), ncol = 2) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_BI_layer$infectious.period = c(4.42,3.40) #Duration infectious period 

generate_r_script(param_list_BI_layer)
# Scenario: vaccination BIBOOSTER
param_list_BIBOOSTER_layer <- param_list_baseline_layer
param_list_BIBOOSTER_layer$scenario <- "BIBOOSTER_layer_farm"
param_list_BIBOOSTER_layer$trans.mean.buildup <- 41.403477
param_list_BIBOOSTER_layer$trans.var.buildup <- 844.1632
param_list_BIBOOSTER_layer$p.protect <- 0.977
param_list_BIBOOSTER_layer$beta = matrix(c(2.13, 2.13,0.21,0.21), ncol = 2) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_BIBOOSTER_layer$infectious.period = c(4.42,3.40) #Duration infectious period 

generate_r_script(param_list_BIBOOSTER_layer)

# Scenario: vaccination CEVA preproduction
param_list_CEVA_preproduction <- param_list_CEVA_layer
param_list_CEVA_preproduction$scenario <- "CEVA_preproduction"
param_list_CEVA_preproduction$length.round <- 120 #length of the round. Simulation starts at time 0
param_list_CEVA_preproduction$age_at_d0 <-0 #age at day 0 = 0
generate_r_script(param_list_CEVA_preproduction)

# Scenario: vaccination BI preproduction
param_list_BI_preproduction <- param_list_BI_layer
param_list_BI_preproduction$scenario <- "BI_preproduction"
param_list_BI_preproduction$length.round <- 120 #length of the round. Simulation starts at time 0
param_list_BI_preproduction$age_at_d0 <-0 #age at day 0 = 0
generate_r_script(param_list_BI_preproduction)

# Scenario: vaccination BIBOOSTER preproduction
param_list_BIBOOSTER_preproduction <- param_list_BIBOOSTER_layer
param_list_BIBOOSTER_preproduction$scenario <- "BIBOOSTER_preproduction"
param_list_BIBOOSTER_preproduction$length.round <- 120 #length of the round. Simulation starts at time 0
param_list_BIBOOSTER_preproduction$age_at_d0 <-0 #age at day 0 = 0
generate_r_script(param_list_BIBOOSTER_preproduction)

# introduction half way preporduction period ####
#Scenario baseline preproduction farm ####
param_list_baseline_preproduction_half <-param_list_baseline_preproduction
param_list_baseline_preproduction_half$scenario <-"baseline_preproduction_half"
param_list_baseline_preproduction_half$intro.time <- 60 #introduction at half way the preproduction period
generate_r_script(param_list_baseline_preproduction_half)

# Scenario: vaccination CEVA preproduction
param_list_CEVA_preproduction_half <- param_list_CEVA_preproduction
param_list_CEVA_preproduction_half$scenario <- "CEVA_preproduction_half"
param_list_CEVA_preproduction_half$intro.time <- 60 #introduction at half way the preproduction period
generate_r_script(param_list_CEVA_preproduction_half)

# Scenario: vaccination BI preproduction
param_list_BI_preproduction_half <- param_list_BI_preproduction
param_list_BI_preproduction_half$scenario <- "BI_preproduction_half"
param_list_BI_preproduction_half$intro.time <- 60 
generate_r_script(param_list_BI_preproduction_half)

# Scenario: vaccination BIBOOSTER preproduction
param_list_BIBOOSTER_preproduction_half <- param_list_BIBOOSTER_preproduction
param_list_BIBOOSTER_preproduction_half$scenario <- "BIBOOSTER_preproduction_half"
param_list_BIBOOSTER_preproduction_half$intro.time <- 60
generate_r_script(param_list_BIBOOSTER_preproduction_half)


#distinct different virus
#assume a new virus with less protection against transmission
# consider a shift of 2-log titre based on ... 

#only for BIBOOSTER and CEVA
# Scenario: vaccination CEVA distinct virus
param_list_CEVA_distinct_layer <- param_list_CEVA_layer
param_list_CEVA_distinct_layer$scenario <- "CEVA_layer_farm_distinct_virus"
param_list_CEVA_distinct_layer$p.protect <- 0.669
generate_r_script(param_list_CEVA_distinct_layer)

# Scenario: vaccination BIBOOSTER distinct virus
param_list_BIBOOSTER_distinct_layer <- param_list_BIBOOSTER_layer
param_list_BIBOOSTER_distinct_layer$scenario <- "BIBOOSTER_layer_farm_distinct_virus"
param_list_BIBOOSTER_distinct_layer$p.protect <- 0.797
generate_r_script(param_list_BIBOOSTER_distinct_layer)

# Sensitivity analyses different cut-off values ####
# cut-off <= 6
# p protect
# # A tibble: 5 × 4
# Group Treatment          mean_protected_co6 mean_protected_co7
# <chr> <chr>                           <dbl>              <dbl>
#   1 1     CEVA                            0.862              0.669
# 2 2     BI                              0.798              0.599
# 3 3     BI+booster                      0.924              0.797

#build-up parameters

# Group     Mean  Variance
# 1      CEVA 15.67698  716.7565
# 2        BI 66.51786 2354.5685
# 3 BIBOOSTER 61.14218 1185.4970


param_list_CEVA06_layer_farm <-param_list_baseline_layer
param_list_CEVA06_layer_farm$scenario <-"CEVA06_layer_farm"
param_list_CEVA06_layer_farm$p.protect <- 0.862
param_list_CEVA06_layer_farm$trans.mean.buildup<-15.67698
param_list_CEVA06_layer_farm$trans.var.buildup <- 716.7565
param_list_CEVA06_layer_farm$beta <- matrix(c(1.92, 1.92,0.15,0.15), ncol = 2) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_CEVA06_layer_farm$infectious.period = c(4.33,3.70) #Duration infectious period 

generate_r_script(param_list_CEVA06_layer_farm)



param_list_BI06_layer_farm <-param_list_baseline_layer
param_list_BI06_layer_farm$scenario <-"BI06_layer_farm"
param_list_BI06_layer_farm$p.protect <- 0.798
param_list_BI06_layer_farm$trans.mean.buildup<-66.51786
param_list_BI06_layer_farm$trans.var.buildup <- 2354.5685
param_list_BI06_layer_farm$beta <- matrix(c(2.07, 2.07,0.17,0.17), ncol = 2) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_BI06_layer_farm$infectious.period = c(4.29,3.35) #Duration infectious period 

generate_r_script(param_list_BI06_layer_farm)



param_list_BIBOOSTER06_layer_farm <-param_list_baseline_layer
param_list_BIBOOSTER06_layer_farm$scenario <-"BIBOOSTER06_layer_farm"
param_list_BIBOOSTER06_layer_farm$p.protect <- 0.924
param_list_BIBOOSTER06_layer_farm$trans.mean.buildup<-15.67698
param_list_BIBOOSTER06_layer_farm$trans.var.buildup <- 716.7565
param_list_BIBOOSTER06_layer_farm$beta <- matrix(c(2.07, 2.07,0.17,0.17), ncol = 2) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_BIBOOSTER06_layer_farm$infectious.period = c(4.29,3.35) #Duration infectious period 

generate_r_script(param_list_BIBOOSTER06_layer_farm)

# cut-off <= 7
# p protect
# # A tibble: 5 × 4
# Group Treatment          mean_protected_co6 mean_protected_co7
# <chr> <chr>                           <dbl>              <dbl>
#   1 1     CEVA                          0.862              0.669
# 2 2     BI                              0.798              0.599
# 3 3     BI+booster                      0.924              0.797

#build-up parameters

# Group      Mean Variance
# 1      CEVA  50.05946 3404.345
# 2        BI 112.46372 8161.928
# 3 BIBOOSTER  79.65847 1451.801

param_list_CEVA07_layer_farm <-param_list_baseline_layer
param_list_CEVA07_layer_farm$scenario <-"CEVA07_layer_farm"
param_list_CEVA07_layer_farm$p.protect <- 0.669
param_list_CEVA07_layer_farm$trans.mean.buildup<-50.05946
param_list_CEVA07_layer_farm$trans.var.buildup <- 3404.345
param_list_CEVA07_layer_farm$beta <- matrix(c(1.22, 1.22,0.14,0.14), ncol = 2) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_CEVA07_layer_farm$infectious.period = c(4.09,3.66) #Duration infectious period 

generate_r_script(param_list_CEVA07_layer_farm)

param_list_BI07_layer_farm <-param_list_baseline_layer
param_list_BI07_layer_farm$scenario <-"BI07_layer_farm"
param_list_BI07_layer_farm$p.protect <- 0.599
param_list_BI07_layer_farm$trans.mean.buildup<-112.46372
param_list_BI07_layer_farm$trans.var.buildup <- 8161.928
param_list_BI07_layer_farm$beta <- matrix(c(2.13, 2.13,0.11,0.11), ncol = 2) #transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_BI07_layer_farm$infectious.period = c(3.68,3.40) #Duration infectious period 

generate_r_script(param_list_BI07_layer_farm)

param_list_BIBOOSTER07_layer_farm <-param_list_baseline_layer
param_list_BIBOOSTER07_layer_farm$scenario <-"BIBOOSTER07_layer_farm"
param_list_BIBOOSTER07_layer_farm$p.protect <- 0.797
param_list_BIBOOSTER07_layer_farm$trans.mean.buildup<-79.65847
param_list_BIBOOSTER07_layer_farm$trans.var.buildup <- 1451.801
param_list_BIBOOSTER07_layer_farm$beta <- matrix(c(2.13, 2.13,0.11,0.11), ncol = 2)#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976  
param_list_BIBOOSTER07_layer_farm$infectious.period = c(3.68,3.40) #Duration infectious period 

generate_r_script(param_list_BIBOOSTER07_layer_farm)
