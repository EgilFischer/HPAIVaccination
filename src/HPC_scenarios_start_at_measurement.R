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
generate_r_script <- function(param_list, sub_dir = NULL, scenario_seed = 12345) {
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
  file_name <- paste0("./HPC_scenarios/",sub_dir,"script_", format(Sys.Date(), "%Y_%m_%d_"), param_list$scenario, ".R")
  
  # Write the content to the file
  writeLines(script_content, con = file_name)
  
  # Return the name of the generated file
  return(file_name)
}

#create bash file to can be run as array job
bash_script <- function(dir)
{
  script_content <- paste0(
    "#!/bin/bash","\n",
    "#SBATCH --job-name=RunRScriptArray","\n",
    "#SBATCH --output=slurm_output_%A_%a.out","\n",
    "#SBATCH --error=slurm_error_%A_%a.err","\n",
    "#SBATCH --nodes=1","\n",
    "#SBATCH --ntasks-per-node=1","\n",
    "#SBATCH --time=02:30:00","\n",
    "#SBATCH --mem=5G","\n",
    "","\n",
    "# Define the directory where R scripts are located","\n",
    "R_SCRIPTS_DIR=\"/home/uu_vet_te/efischer/HPAI_vaccination/scenario_files/\" ","\n",
    "","\n",
    "# List of scenario files (adjust as needed)","\n",
    "SCENARIO_FILES=(\"",paste(list.files(path = dir, pattern ="R"),collapse = "\" \"")," \")",
    "","\n",
    "# Use the job array index to select the correct file","\n",
    "SCENARIO_INDEX=$SLURM_ARRAY_TASK_ID","\n",
    "R_SCRIPT=\"${R_SCRIPTS_DIR}${SCENARIO_FILES[$SCENARIO_INDEX]}\"" ,"\n",
    "","\n",
    "# Check if the R script exists","\n",
    "if [ ! -f \"$R_SCRIPT\" ]; then","\n",
    "echo \"Error: File $R_SCRIPT not found.\"","\n",
    "exit 1","\n",
    "fi","\n",
    "","\n",
    "# Run the R script" ,"\n",
    "Rscript \"$R_SCRIPT\" ")
  
  # Create a file name with path to scripts
  file_name <- paste0(dir,"/filled_run_script.sh")
  
  # Write the content to the file
  writeLines(script_content, con = file_name)
  
  # Return the name of the generated file
  return(list(file_name,paste("Number of files: ", length(list.files(path = dir, pattern ="R"))-1)))
}


# load default parameters
if(!exists("default_parameter_list")){source("src/default_parameter_list.R")}

#load table with proportion high titre at different measurements
if(!exists("data_field_study_cuttoffs")) source("src/ProcessFieldExperiment.R")

#create a named list with parameters per treatment ####
#for cut_off >= 6
list_vaccin_specs_co6 <- list(
  BIBOOSTER = list(beta = matrix(c(2.07, 2.07,0.17,0.17), ncol = 2),
                 infectious.period =  c(4.29,3.35)),
  BI =list(beta = matrix(c(2.07, 2.07,0.17,0.17), ncol = 2),
           infectious.period =  c(4.29,3.35)),
  CEVA =list(beta = matrix(c(1.92, 1.92,0.15,0.15), ncol = 2),
             infectious.period =  c(4.33,3.70))
  
)

#for cut_off >= 5
list_vaccin_specs_co5 <- list(
  BIBOOSTER = list(beta = matrix(c(2.13, 2.13,0.21,0.21), ncol = 2),
                   infectious.period =  c(4.42,3.40)),
  BI =list(beta = matrix(c(2.13, 2.13,0.21,0.21), ncol = 2),
           infectious.period =  c(4.42,3.40)),
  CEVA =list(beta = matrix(c(1.94, 1.94,0.25,0.25), ncol = 2),
             infectious.period =  c(4.18,3.86))
  
)
#for cut_off >= 7
list_vaccin_specs_co7 <- list(
  BIBOOSTER = list(beta = matrix(c(2.06, 2.06,0.11,0.11), ncol = 2),
                   infectious.period =  c(3.68,3.40)),
  BI =list(beta = matrix(c(2.06, 2.06,0.11,0.11), ncol = 2),
           infectious.period =  c(3.68,3.40)),
  CEVA =list(beta = matrix(c(1.22, 1.22,0.14,0.14), ncol = 2),
             infectious.period =  c(4.09,3.66))
  
)

#for no vaccination - here we use the average low-titre of cut-off >= 5 because this has least dilution of protected animals
list_vaccin_specs_novac <- list(
  NOVAC = list(beta = matrix(c(2.035, 2.035,2.035,2.035), ncol = 2),
                   infectious.period =  c(3.885,3.885))
)

#iterate over the table with measurements per treatment ####
param_list <- default_parameter_list;
#constant immunity
param_list$trans.mean.wane <- Inf;
param_list$trans.var.wane <- Inf;
param_list$trans.mean.buildup <- Inf;
param_list$trans.var.buildup <- Inf;
param_list$pdie <- c(0.2,0.01); #only 20% vaccinated die
param_list$N0 <- 46000; #medium sized farm
param_list$intro.time <- 0; #immediate introduction

#for cut_off 6
for(TREAT in c("BI","BIBOOSTER","CEVA")){
  #specific values for this vaccine
  param_list$beta <- list_vaccin_specs_co6[[TREAT]]$beta
  param_list$infectious.period <- list_vaccin_specs_co6[[TREAT]]$infectious.period
  #set moment specific values by iterating over the measurements
  measurements_high_titre <- data_field_study_cut_offs%>%
    filter(Group!=4 & Treatment == TREAT)
  for(i in c(1:nrow(measurements_high_titre))){
    param_list$p.protect <- as.numeric(measurements_high_titre[i,"O6"])
    param_list$age_at_d0<- as.numeric(measurements_high_titre[i,"time_point"])
    param_list$scenario <-paste(TREAT,"CO6","Time", measurements_high_titre[i,"time_point"],"Farm",measurements_high_titre[i,"Farm"],sep = "_")
    #create r script
    generate_r_script(param_list, sub_dir = "CO6/")
  }
}
bash_script("C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_scenarios/CO6")
#for cut_off 5
for(TREAT in c("BI","BIBOOSTER","CEVA")){
  #specific values for this vaccine
  param_list$beta <- list_vaccin_specs_co5[[TREAT]]$beta
  param_list$infectious.period <- list_vaccin_specs_co5[[TREAT]]$infectious.period
  #set moment specific values by iterating over the measurements
  measurements_high_titre <- data_field_study_cut_offs%>%
    filter(Group!=4 & Treatment == TREAT)
  for(i in c(1:nrow(measurements_high_titre))){
    param_list$p.protect <- as.numeric(measurements_high_titre[i,"O5"])
    param_list$age_at_d0<- as.numeric(measurements_high_titre[i,"time_point"])
    param_list$scenario <-paste(TREAT,"CO5","Time", measurements_high_titre[i,"time_point"],"Farm",measurements_high_titre[i,"Farm"],sep = "_")
    #create r script
    generate_r_script(param_list, sub_dir = "CO5/")
  }
}

bash_script("C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_scenarios/CO5")


#for cut_off 7
for(TREAT in c("BI","BIBOOSTER","CEVA")){
  #specific values for this vaccine
  param_list$beta <- list_vaccin_specs_co7[[TREAT]]$beta
  param_list$infectious.period <- list_vaccin_specs_co7[[TREAT]]$infectious.period
  #set moment specific values by iterating over the measurements
  measurements_high_titre <- data_field_study_cut_offs%>%
    filter(Group!=4 & Treatment == TREAT)
  for(i in c(1:nrow(measurements_high_titre))){
    param_list$p.protect <- as.numeric(measurements_high_titre[i,"O7"])
    param_list$age_at_d0<- as.numeric(measurements_high_titre[i,"time_point"])
    param_list$scenario <-paste(TREAT,"CO7","Time", measurements_high_titre[i,"time_point"],"Farm",measurements_high_titre[i,"Farm"],sep = "_")
    #create r script
    generate_r_script(param_list, sub_dir = "CO7/")
  }
}

bash_script("C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_scenarios/CO7")


#for no vaccination
param_list$pdie <- c(1.0,1.0); #unvaccinated will always die
for(TREAT in c("NOVAC")){
  #specific values for this vaccine
  param_list$beta <- list_vaccin_specs_novac[[TREAT]]$beta
  param_list$infectious.period <- list_vaccin_specs_novac[[TREAT]]$infectious.period
  #set moment specific values by iterating over the measurements
  
  measurements_high_titre <- data_field_study_cut_offs%>%
    filter( Group == 4)
  #only need is the time points
  time_points <- unique(measurements_high_titre[,c("Farm","time_point")])
  
  
  for(i in c(1:nrow(time_points))){
    param_list$p.protect <- 0
    param_list$age_at_d0<- as.numeric(time_points[i,"time_point"])
    param_list$scenario <-paste("NOVAC","Time", time_points[i,"time_point"],"Farm",time_points[i,"Farm"],sep = "_")
    #create r script
    generate_r_script(param_list, sub_dir = "NOVAC/")
  }
}

bash_script("C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_scenarios/NOVAC")
