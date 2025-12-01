#########################################################
#                                                        
#                 Data post-processing                           
#                   
#########################################################
#Input: DailyData and Detection Results 

source("./src/loadLibraries.R") 
library(scales)

#output directory
output_dir <- "./output/"
#make folder for saving figures
if(!dir.exists("./figures")){dir.create("./figures")}

############################################
##  General plotting setup for DailyData  ##
############################################
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
    return(df)
  })
  combined_data <- do.call(rbind, scenario_data)
  return(combined_data)
}

#output directory
output_dir <- "./HPC_output/"

##################################################
##  1. Plot infection dynamics per scenario   ####
##################################################
#output: Line plots of SLIR compartments over time for each scenario (within scenarios: one tile per rds file)

#load data for each scenario as data frames
scenario_preproduction_df <- load_scenario_data(pattern = "*preproduction*\\.RDS$", output_dir = output_dir)
scenario_layer_farm_df <- load_scenario_data(pattern = "*layer_farm*\\.RDS$", output_dir = output_dir)

### Plotting Scenario's preproduction ###
scenario_preproduction_long <- scenario_preproduction_df %>%select(
    time, run, scenario, S, L, I, R
  ) %>%
  pivot_longer(
    cols = c(S, L, I, R),
    names_to = "variable",
    values_to = "value"
  )

scenario_preproduction_plot <- ggplot(scenario_preproduction_long, aes(x = time, y = value, 
                                                                       colour = variable, 
                                                                       group = interaction(run, variable))) +
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

ggsave(filename = "./figures/scenario_preproduction_plot.png", plot = scenario_preproduction_plot, width = 10, height = 8, dpi = 300)


### Plotting Scenario 2.1 ###
scenario2_1_long <- scenario2_1_df %>%
  pivot_longer(
    cols = c(S, L, I, R),
    names_to = "variable",
    values_to = "value"
  )

scenario2_1_plot <- ggplot(scenario2_1_long, aes(x = day, y = value, colour = variable, group = interaction(run, variable))) +
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


### Plotting Scenario 2.2 ###
scenario2_2_long <- scenario2_2_df %>%
  pivot_longer(
    cols = c(S, L, I, R),
    names_to = "variable",
    values_to = "value"
  )

scenario2_2_plot <- ggplot(scenario2_2_long, aes(x = day, y = value, colour = variable, group = interaction(run, variable))) +
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

###########################
##    2. Plot eggs     ###
###########################
#output: on-farm egg dynamics for I and H eggs over time per simulation (messy if you do my scenario)
plot_egg_dynamics <- function(data, scenario_title) {
  egg_long <- data %>%
    select(day, run, Egg_I_farm, Egg_H_farm) %>%
    pivot_longer(
      cols = c(Egg_I_farm, Egg_H_farm),
      names_to = "Egg_Type",
      values_to = "Egg_Count"
    )
  
  # Create the egg dynamics plot
  egg_plot <- ggplot(egg_long, aes(x = day, y = Egg_Count, colour = Egg_Type,
                                   group = interaction(run, Egg_Type))) +
    geom_line(alpha = 0.5) +
    labs(
      title = scenario_title,
      x = "Day",
      y = "Egg Count",
      colour = "Egg Type"
    ) +
    theme_minimal()
  
  return(egg_plot)
}

# Load data for each scenario !better to load individual not entire scenarios 
egg_df <- load_scenario_data(pattern = "DailyData.2025_03_30_S2.1_HomoWaning_IntroTime_504_LayerSize_46000.*\\.RDS$", output_dir = output_dir)

# Plot egg dynamics 
egg_plots <- plot_egg_dynamics(egg_df, "Egg Dynamics: Scenario 2.1 (HomoWaning) IntroTime X")
ggsave(filename = "./figures/egg_dynamics_scenario2.1.png", plot = egg_plots,
       width = 10, height = 8, dpi = 300)

############################
##  3. Plot final size   ###
############################
#load and process .RDS files for all runs
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

#extract parameters from one .RDS file for a scenario
extract_scenario_parameters <- function(pattern, output_dir = "./output/") {
  rds_files <- list.files(output_dir, pattern = pattern, full.names = TRUE)
  if (length(rds_files) == 0) stop("No files found for the pattern!")
  
  #extract parameters
  param_list <- readRDS(rds_files[1])$pars
  return(param_list)
}

#extract parameters for each scenario
scenario1_params <- extract_scenario_parameters(pattern = "Daily\\Data.*S1_LayerSize.*\\.RDS$")
scenario2_1_params <- extract_scenario_parameters(pattern = "Daily\\Data.*HomoWaning.*\\.RDS$")
scenario2_2_params <- extract_scenario_parameters(pattern = "Daily\\Data.*HeteroWaning.*\\.RDS$")

#load and process data for each scenario
scenario1_fs <- load_and_process_final_size(pattern = "Daily\\Data.*S1_LayerSize.*\\.RDS$")
scenario2_1_fs <- load_and_process_final_size(pattern = "Daily\\Data.*HomoWaning.*\\.RDS$")
scenario2_2_fs <- load_and_process_final_size(pattern = "Daily\\Data.*HeteroWaning.*\\.RDS$")

#plotting function for histograms
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

# Plot Scenario 1 
scenario1_histo <- plot_tile_histograms(
  scenario1_fs,
  param_list = list(N0 = scenario1_params$N0, runs = scenario1_params$runs),
  title = "Distribution of Outbreak Sizes: Scenario 1"
)
ggsave(filename = "./figures/scenario1_histo.png", plot = scenario1_histo, width = 10, height = 8, dpi = 300)

# Plot Scenario 2.1 
scenario2.1_histo <- plot_tile_histograms(
  scenario2_1_fs,
  param_list = list(N0 = scenario2_1_params$N0, runs = scenario2_1_params$runs),
  title = "Distribution of Outbreak Sizes: Scenario 2.1"
)
ggsave(filename = "./figures/scenario2.1_histo.png", plot = scenario2.1_histo, width = 10, height = 8, dpi = 300)

# Plot Scenario 2.2 
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

############################################
## 4. Detection Method by outbreak size  ##
###########################################
#output: proportion of outbreaks (minor and major) detected by active vs passive vs not detected 

# Load files S2.1 and S2.2 detection module output
detection_file_S2.1 <- list.files(output_dir, pattern = "Detection_results_complete_S2.1_Homo\\.RDS$", full.names = TRUE)
detection_file_S2.2 <- list.files(output_dir, pattern = "Detection_results_complete_S2.2_Hetero\\.RDS$", full.names = TRUE)

detection_S2.1_df <- readRDS(detection_file_S2.1)
detection_S2.2_df <- readRDS(detection_file_S2.2)

plot_detection_proportions <- function(df, scenario_name) {
  df %>%
    #filter(outbreak.size == "MAJOR") %>%  #if only want to look at major
    mutate(detection.method = case_when(
      is.infinite(min.det.time) ~ "Undetected",
      det.since.intro.ac < det.since.intro.pas ~ "Active",
      det.since.intro.ac > det.since.intro.pas ~ "Passive",
      det.since.intro.ac == det.since.intro.pas ~ "No difference",
      TRUE ~ "Other"
    )) %>%
    group_by(outbreak.size, detection.method) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(outbreak.size) %>%
    mutate(proportion = count / sum(count)) %>%
    ggplot(aes(x = outbreak.size, y = proportion, fill = detection.method)) +
    geom_col(position = "fill") +
    labs(
      title = paste("Detection Method Proportions –", scenario_name),
      x = "Outbreak Size",
      y = "Proportion",
      fill = "Detection Method"
    ) +
    theme_minimal()
}
plot_S2.1 <- plot_detection_proportions(detection_S2.1_df, "S2.1 Homo")
ggsave(filename = "./figures/detection_proportions_S2.1_Homo.png", plot = plot_S2.1, width = 8, height = 6, dpi = 300)

plot_S2.2 <- plot_detection_proportions(detection_S2.2_df, "S2.2 Hetero")
ggsave(filename = "./figures/detection_proportions_S2.2_Hetero.png",plot = plot_S2.2, width = 8, height = 6, dpi = 300)




######################################
##  5. Boxplots: active detection  ##
#####################################
detection_file_S2.1 <- list.files(output_dir, pattern = "Detection_results_complete_S2.1_Homo\\.RDS$", full.names = TRUE)
detection_file_S2.2 <- list.files(output_dir, pattern = "Detection_results_complete_S2.2_Hetero\\.RDS$", full.names = TRUE)

detection_S2.1_df <- readRDS(detection_file_S2.1)
detection_S2.2_df <- readRDS(detection_file_S2.2)

plot_detection_boxplot <- function(df, scenario_name, filename) {
  df_major <- df %>%
    filter(outbreak.size == "MAJOR") %>%
    mutate(
      det.time.ac.capped = ifelse(is.infinite(det.since.intro.ac), 100, det.since.intro.ac),
      time.interval.ac = factor(time.interval.ac, levels = c(7, 14, 30))
    )
  
  intro_counts <- df_major %>%
    count(intro.time, name = "n_outbreaks")
  
  df_major_labeled <- df_major %>%
    left_join(intro_counts, by = "intro.time") %>%
    mutate(
      intro.time.labeled = paste0(intro.time, " (n=", n_outbreaks, ")"),
      intro.time.labeled = factor(intro.time.labeled, levels = unique(intro.time.labeled))
    )
  
  p <- ggplot(df_major_labeled, aes(x = intro.time.labeled, y = det.time.ac.capped, fill = time.interval.ac)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    scale_fill_manual(values = c("7" = "#1f78b4", "14" = "#a6cee3", "30" = "#b2df8a")) +
    labs(
      title = paste("Detection Time Distribution by Introduction Time and Sampling Interval"),
      subtitle = paste("Major outbreaks only (", scenario_name, "); Inf = 100 days", sep = ""),
      x = "Introduction Time",
      y = "Detection Time (days)",
      fill = "Sampling Interval"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save the plot
  ggsave(filename = filename, plot = p, width = 10, height = 6, dpi = 300)
}

plot_detection_boxplot(detection_S2.1_df, "S2.1 Homo", "./figures/boxplot_detection_S2.1_Homo.png")
plot_detection_boxplot(detection_S2.2_df, "S2.2 Hetero", "./figures/boxplot_detection_S2.2_Hetero.png")


#####################################################
##  6. Major outbreak detection by sampling freq  ###
#####################################################
detection_file_S2.1 <- list.files(output_dir, pattern = "Detection_results_complete_S2.1_Homo\\.RDS$", full.names = TRUE)
detection_file_S2.2 <- list.files(output_dir, pattern = "Detection_results_complete_S2.2_Hetero\\.RDS$", full.names = TRUE)

detection_S2.1_df <- readRDS(detection_file_S2.1)
detection_S2.2_df <- readRDS(detection_file_S2.2)

plot_detection_proportions <- function(df, scenario_name) {
  df %>%
    filter(outbreak.size == "MAJOR") %>% #hash out to have separate bar for minor
    mutate(detection.method = case_when(
      is.infinite(min.det.time) ~ "Undetected",
      det.since.intro.ac < det.since.intro.pas & time.interval.ac == 7  ~ "Active 7",
      det.since.intro.ac < det.since.intro.pas & time.interval.ac == 14 ~ "Active 14",
      det.since.intro.ac < det.since.intro.pas & time.interval.ac == 30 ~ "Active 30",
      det.since.intro.ac > det.since.intro.pas ~ "Passive",
      det.since.intro.ac == det.since.intro.pas ~ "No difference",
      TRUE ~ "Other"
    )) %>%
    group_by(outbreak.size, detection.method) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(outbreak.size) %>%
    mutate(proportion = count / sum(count)) %>%
    ggplot(aes(x = outbreak.size, y = proportion, fill = detection.method)) +
    geom_col(position = "fill") +
    scale_fill_manual(
      values = c(
        "Active 7" = "#1b9e77",
        "Active 14" = "#66c2a5",
        "Active 30" = "#a6dba0",
        "Passive" = "#a6cee3",
        "Undetected" = "#999999",
        "No difference" = "#e7298a"
      )
    ) +
    labs(
      title = paste("Detection Method Proportions –", scenario_name),
      x = "Outbreak Size",
      y = "Proportion",
      fill = "Detection Method"
    ) +
    theme_minimal()
}

plot_S2.1 <- plot_detection_proportions(detection_S2.1_df, "S2.1 Homo")
ggsave(filename = "./figures/detection_proportions_S2.1_Homo.png", plot = plot_S2.1, width = 8, height = 6, dpi = 300)

plot_S2.2 <- plot_detection_proportions(detection_S2.2_df, "S2.2 Hetero")
ggsave(filename = "./figures/detection_proportions_S2.2_Hetero.png", plot = plot_S2.2, width = 8, height = 6, dpi = 300)


# #####################################################
# ##  7. Major outbreak detection by passive method ###
# #####################################################
##output: threshold vs ratio "winner" 
# detection_file_S2.1 <- list.files(output_dir, pattern = "Detection_results_complete_S2.1_Homo\\.RDS$", full.names = TRUE)
# detection_file_S2.2 <- list.files(output_dir, pattern = "Detection_results_complete_S2.2_Hetero\\.RDS$", full.names = TRUE)
# 
# detection_S2.1_df <- readRDS(detection_file_S2.1)
# detection_S2.2_df <- readRDS(detection_file_S2.2)
# #S2.1
# df_major <- detection_S2.1_df %>%
#   filter(outbreak.size == "MAJOR")
# 
# df_major <- df_major %>%
#   mutate(
#     pas.winner = case_when(
#       is.infinite(pas.det.time) ~ "Undetected",
#       pas.det.time.threshold < pas.det.time.ratio ~ "Threshold",
#       pas.det.time.threshold > pas.det.time.ratio ~ "Ratio",
#       pas.det.time.threshold == pas.det.time.ratio ~ "Tie",
#       TRUE ~ "Other"
#     )
#   )
# df_passive_summary <- df_major %>%
#   group_by(intro.time, pas.winner) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(intro.time) %>%
#   mutate(proportion = n / sum(n))
# 
# plot_passive_S2.1 <- ggplot(df_passive_summary, aes(x = factor(intro.time), y = proportion, fill = pas.winner)) +
#   geom_col(position = "stack") +
#   labs(
#     title = "Proportion of Major Outbreaks Detected by Passive Method",
#     subtitle = "(S2.1 Homo)",
#     x = "Introduction Time",
#     y = "Proportion",
#     fill = "Detection Winner"
#   ) +
#   theme_classic() +
#   theme(axis.line = element_line(color = "black"))
# 
# ggsave(filename = "./figures/passive_detection_winners_S2.1.png", plot = plot_passive_S2.1, width = 8, height = 6, dpi = 300)
# 
# #S2.2
# df_major_S2.2 <- detection_S2.2_df %>%
#   filter(outbreak.size == "MAJOR")
# 
# # Determine passive winner
# df_major_S2.2 <- df_major_S2.2 %>%
#   mutate(
#     pas.winner = case_when(
#       is.infinite(pas.det.time) ~ "Undetected",
#       pas.det.time.threshold < pas.det.time.ratio ~ "Threshold",
#       pas.det.time.threshold > pas.det.time.ratio ~ "Ratio",
#       pas.det.time.threshold == pas.det.time.ratio ~ "Tie",
#       TRUE ~ "Other"
#     )
#   )
# 
# # Summarize proportions
# df_passive_summary_S2.2 <- df_major_S2.2 %>%
#   group_by(intro.time, pas.winner) %>%
#   summarise(n = n(), .groups = "drop") %>%
#   group_by(intro.time) %>%
#   mutate(proportion = n / sum(n))
# 
# # Make the plot
# plot_passive_S2.2 <- ggplot(df_passive_summary_S2.2, aes(x = factor(intro.time), y = proportion, fill = pas.winner)) +
#   geom_col(position = "stack") +
#   labs(
#     title = "Proportion of Major Outbreaks Detected by Passive Method",
#     subtitle = "(S2.2 Hetero)",
#     x = "Introduction Time",
#     y = "Proportion",
#     fill = "Detection Winner"
#   ) +
#   theme_classic() +
#   theme(axis.line = element_line(color = "black"))
# 
# # Save the plot
# ggsave(filename = "./figures/passive_detection_winners_S2.2.png", plot = plot_passive_S2.2, width = 8,height = 6, dpi = 300)

###############################################################
##  8. Different Major outbreak detection by passive method  ###
################################################################
#output: box and whisker of threshold vs ratio detection time

detection_file_S2.1 <- list.files(output_dir, pattern = "Detection_results_complete_S2.1_Homo\\.RDS$", full.names = TRUE)
detection_file_S2.2 <- list.files(output_dir, pattern = "Detection_results_complete_S2.2_Hetero\\.RDS$", full.names = TRUE)

detection_S2.1_df <- readRDS(detection_file_S2.1)
detection_S2.2_df <- readRDS(detection_file_S2.2)

plot_passive_detection_boxplot <- function(df, scenario_name, filename) {
  df_major <- df %>%
    filter(outbreak.size == "MAJOR") %>%
    select(intro.time, pas.det.time.threshold, pas.det.time.ratio) %>%
    mutate(
      threshold = ifelse(is.infinite(pas.det.time.threshold), 100, pas.det.time.threshold),
      ratio = ifelse(is.infinite(pas.det.time.ratio), 100, pas.det.time.ratio)
    ) %>%
    pivot_longer(cols = c(threshold, ratio), names_to = "method", values_to = "detection.time")
  
  intro_counts <- df_major %>%
    count(intro.time, name = "n_outbreaks")
  
  df_labeled <- df_major %>%
    left_join(intro_counts, by = "intro.time") %>%
    mutate(
      intro.time.labeled = paste0(intro.time, " (n=", n_outbreaks, ")"),
      intro.time.labeled = factor(intro.time.labeled, levels = unique(intro.time.labeled)),
      method = factor(method, levels = c("threshold", "ratio"))
    )
  
  p <- ggplot(df_labeled, aes(x = intro.time.labeled, y = detection.time, fill = method)) +
    geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
    scale_fill_manual(values = c("threshold" = "#1f78b4", "ratio" = "#b2df8a"),
                      labels = c("Threshold", "Ratio")) +
    labs(
      title = "Passive Detection Time by Introduction Time and Method",
      subtitle = paste("Major outbreaks only (", scenario_name, "); Inf = 100 days", sep = ""),
      x = "Introduction Time",
      y = "Detection Time (days)",
      fill = "Passive Detection Method"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(filename = filename, plot = p, width = 10, height = 6, dpi = 300)
}

plot_passive_detection_boxplot(detection_S2.1_df, "S2.1 Homo", "./figures/boxplot_passive_S2.1_Homo.png")
plot_passive_detection_boxplot(detection_S2.2_df, "S2.2 Hetero", "./figures/boxplot_passive_S2.2_Hetero.png")



################################################
##  9. Plot eggs on farm/shipped by min det  ###
################################################
detection_file_S2.1 <- list.files(output_dir, pattern = "Detection_results_complete_S2.1_Homo\\.RDS$", full.names = TRUE)
detection_file_S2.2 <- list.files(output_dir, pattern = "Detection_results_complete_S2.2_Hetero\\.RDS$", full.names = TRUE)

detection_S2.1_df <- readRDS(detection_file_S2.1)
detection_S2.2_df <- readRDS(detection_file_S2.2)

#S2.1
s2.1 <- ggplot(detection_S2.1_df, aes(x = as.factor(intro.time), y = detect.Egg_I_shipped, fill = outbreak.size)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "Shipped eggs laid by infectious birds at time of earliest detection",
    subtitle = "Scenario 2.1 Homologous waning",
    x = "Introduction Time",
    y = "Infected Eggs Shipped at Detection Time",
    fill = "Outbreak Size"
  ) +
  theme_minimal()
ggsave(filename = "./figures/eggs_shipped_by_introtime_S2.1.png", plot = s2.1,
       width = 10, height = 6, dpi = 300)
#may get warning of "Removed x rows containing non-finite values" because detection time fell outside sim time ie. no eggs

#S2.2
s2.2 <- ggplot(detection_S2.2_df, aes(x = as.factor(intro.time), y = detect.Egg_I_shipped, fill = outbreak.size)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    title = "Shipped eggs laid by infectious birds at time of earliest detection",
    subtitle = "Scenario 2.2 Heterologous waning",
    x = "Introduction Time",
    y = "Infected Eggs Shipped at Detection Time",
    fill = "Outbreak Size"
  ) +
  theme_minimal()
ggsave(filename = "./figures/eggs_shipped_by_introtime_S2.2.png", plot = s2.2,
       width = 10, height = 6, dpi = 300)
