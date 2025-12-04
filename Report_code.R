### "Surveillance in a vaccinated poultry population" ####
# Author e.a.j.fischer@uu.nl 

# set working directories and directories with simulation outputs ####
setwd("C:/Users/fisch106/SurfDrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination_private_version")
hpc_directory <- "./HPC_output"
output.dir.detection <- "./HPC_output/Detection/";
output.dir.daily_data <- "./HPC_output/DailyData/";
  
# Source required scripts that load and process the data ####
  
## load / set libraries, functions and labels####
#run scripts to load required libraries and calculation of probability of a major outbreak
suppressPackageStartupMessages(source("./src/loadLibraries.R"))
suppressMessages(source("./src/pmaj_calculation.R"))

#function to calculate clopper-pearson confidence interval
clopper_pearson <- function(x, n, conf.level = 0.95) {
  alpha <- 1 - conf.level
  lower <- qbeta(alpha / 2, x, n - x + 1)
  upper <- qbeta(1 - alpha / 2, x + 1, n - x)
  return(c(lower, upper))
}

# Okabe-Ito palette = coluor blind friendly
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")

#registered symbol
reg_symbol = "\U00AE"

#labels
ceva_labels <- c('CEVA' = paste0("VECTORMUNE",reg_symbol),
                 'NOVAC' = "No vaccination",
                 'Neg_controle CEVA' = "Control group")
bi_labels <- c('BI' = "VAXXITEK HVT + IBD + H5",
               'BIBOOSTER' = paste0("VAXXITEK HVT + IBD + H5 \n + Volvac",reg_symbol," B.E.S.T. AI + ND"),
               'NOVAC' = "No vaccination",
               'Neg_controle BI' = "Control group")




#####################################
## run analyses or load analysis ####
#####################################

#set those sections to run. Setting these to false will result in using previous runs of these sections or scripts. This speeds up when processing the data for editing text and figures etc. 
run_detection_module = FALSE
run_process_field_experiment = TRUE
run_combine_daily_data = FALSE

#process field experiment ####
if(run_process_field_experiment){
  suppressMessages(source("./src/ProcessFieldExperiment.R"))
}

#load data of field experiment ####
data_field_study_cut_offs <- readRDS("./input/FieldExperiment/Processed_Data/data_field_study_cut_offs.RDS")


# detection module ####
if(run_detection_module){ #change to FALSE to prevent running this code
  #run analysis
  suppressMessages(source("./src/pre_processing_detModule.R"))
  suppressMessages(source("./src/detectionModule.R"))
}



## load the simulations ####
if(!run_combine_daily_data){
  # Load the data
  file.list <- list.files(output.dir.daily_data, pattern = "DailyData.*\\.RDS$", full.names = TRUE);
  all_scenarios_sim <- NULL
  for(i in 1:length(file.list)){
    #create one data.frame
    next_scenario <- readRDS(file.list[i]);
    all_scenarios_sim <- rbind(all_scenarios_sim, cbind(next_scenario$out, scenario = next_scenario$pars$scenario,
                                                        intro_time = next_scenario$pars$age_at_d0,
                                                        time_point = next_scenario$pars$age_at_d0,
                                                        p_protect = next_scenario$pars$p.protect));
  }
  all_scenarios_sim <- all_scenarios_sim %>%
    mutate(Farm = str_extract(str_extract(scenario,"Farm_A|Farm_B"),c("A|B")),
           Treatment = str_extract(scenario,c("BI_|BIBOOSTER|CEVA|NOVAC")),
           cut_off = str_extract(scenario, "O5|O6|O7"))
  all_scenarios_sim$cut_off <-replace_na(all_scenarios_sim$cut_off, "NOVAC") 
  all_scenarios_sim$Treatment <- str_replace(all_scenarios_sim$Treatment, "BI_","BI")
  
  all_scenarios_sim <- all_scenarios_sim%>%mutate(production_phase = ifelse(time_point<120,"preproduction","production"))
  saveRDS(all_scenarios_sim,file = paste0(output.dir.daily_data, "/Combined/Combined_DailyData.RDS"))
}else{
  # Load the data
  all_scenarios_sim <-readRDS(file = paste0(output.dir.daily_data, "/Combined/Combined_DailyData.RDS"))
}

# get the values of transmission parameters of all the treatments ####
all_files <- data.frame(files = unlist(list.files(hpc_directory, include.dirs = FALSE)))%>%filter(str_detect(files,pattern = "Farm"))
all_files_selection <- data.frame(files = all_files$files,
                                  Treatment = str_extract(all_files$files, "NOVAC|CEVA|BI|BIBOOSTER"),
                                  cut_off = str_extract(all_files$files, "O5|O6|O7"))%>%
  reframe(
    .by = c(Treatment, cut_off),
    file = sample(files,1))

#get parameters from file 
get_parameters_from_file <- function(file_name, path = NULL, parameter = NULL, index = NULL){
  tmp_pars <- readRDS(paste0(path,"/",file_name))$pars
  #get parameter indicated
  if(!is.null(parameter)){
    tmp_pars <- tmp_pars[[parameter]]
    #get value at index 
    if(!is.null(index)){
      tmp_pars <- tmp_pars[index]
    }
  }
  return(tmp_pars)
}

#get values 
all_files_selection <- all_files_selection%>%mutate(beta_low_titre = sapply(all_files_selection$file,
                                                                            FUN = get_parameters_from_file,
                                                                            path = hpc_directory, 
                                                                            parameter = "beta", 
                                                                            index = 1),
                                                    Tinf_low_titre = sapply(all_files_selection$file,
                                                                            FUN = get_parameters_from_file,
                                                                            path = hpc_directory, 
                                                                            parameter = "infectious.period", 
                                                                            index = 1),
                                                    k_Tinf_low_titre = sapply(all_files_selection$file,
                                                                              FUN = get_parameters_from_file,
                                                                              path = hpc_directory, 
                                                                              parameter = "k.infectious", 
                                                                              index = 1),
                                                    beta_high_titre = sapply(all_files_selection$file,
                                                                             FUN = get_parameters_from_file,
                                                                             path = hpc_directory, 
                                                                             parameter = "beta", 
                                                                             index = 3),
                                                    Tinf_high_titre = sapply(all_files_selection$file,
                                                                             FUN = get_parameters_from_file,
                                                                             path = hpc_directory, 
                                                                             parameter = "infectious.period", 
                                                                             index = 2),
                                                    k_Tinf_high_titre = sapply(all_files_selection$file,
                                                                               FUN = get_parameters_from_file,
                                                                               path = hpc_directory, 
                                                                               parameter = "k.infectious", 
                                                                               index =1))

# Create a data frame for the transmission coefficients from parameter list
# create input to the table - check with input!
input_table <- with(all_files_selection,{data.frame(
  Treatment = Treatment,
  cut_off = replace_na(cut_off,"-"),
  cutoff_name = replace_na(str_replace(cut_off,pattern = "O", replacement = " &gt;= "),"&#45;"),
  beta_low_titre = beta_low_titre,
  Tinf_low_titre = Tinf_low_titre,
  Tinf_025p_low_titre = round(qgamma(0.025, shape =2, scale =Tinf_low_titre/k_Tinf_low_titre),digits =2), 
  Tinf_975p_low_titre =  round(qgamma(0.975, shape =2, scale =Tinf_low_titre/k_Tinf_low_titre),digits =2), 
  beta_high_titre = beta_high_titre,
  Tinf_high_titre = Tinf_high_titre,
  Tinf_025p_high_titre =  round(qgamma(0.025, shape =2, scale =Tinf_high_titre/k_Tinf_high_titre),digits =2), 
  Tinf_975p_high_titre =  round(qgamma(0.975, shape =2, scale =Tinf_high_titre/k_Tinf_high_titre/2),digits =2)
)})

row.names(input_table)<-NULL;

#calculate the prob. of major outbreak for different inputs for use in the report ####
# get all input parameters
input_pars <- input_table
#add row(s) for BIBOOSTER 
input_pars <- rbind(input_pars, 
                    input_pars%>%filter(Treatment == "BI")%>%mutate(Treatment = "BIBOOSTER"))

pmajor_scenarios <- NULL;

#for each vaccine and cut-off value determine R and p_major for each sample point
for(i in 1:nrow(input_pars)){
  tmp_param_list <- with(input_pars[i,],{
    #set the essential parameter list
    list(
      N0 = 46000,
      beta = matrix(c(beta_low_titre,beta_low_titre,beta_high_titre,beta_high_titre),ncol =2),
      infectious.period  = c(Tinf_low_titre, Tinf_high_titre),
      k.infectious = c(2,2),
      mortRate = 2e-04
    )
  })
  #take care of the unvaccinated parameter input which is NA
  if(is.na(tmp_param_list$beta[,2])[1]){tmp_param_list$beta[,2]<- tmp_param_list$beta[,1]}
  if(is.na(tmp_param_list$infectious.period[2])){tmp_param_list$infectious.period[2]<- tmp_param_list$infectious.period[1]}
  
  #now determine probability of major outbreak and R for each data point in the field study
  #first select the right values
  treatment = input_pars$Treatment[i]
  tmp_field_data <- data_field_study_cut_offs%>%filter(Treatment ==input_pars$Treatment[i])
  if(input_pars$Treatment[i]== "NOVAC") {
    time_points_farm <- data_field_study_cut_offs%>%filter(Group == 4)
    pmajor_scenarios <-rbind(pmajor_scenarios,cbind(tibble(time_point = time_points_farm$time_point,Farm = time_points_farm$Farm,Group =0,Treatment = "NOVAC",cut_off = "NOVAC"),
                                                    pmajor(tmp_param_list,p = 0, n = 1 )))}else{
                                                      for(j in 1:nrow(tmp_field_data))
                                                      {
                                                        pmajor_scenarios <-rbind(pmajor_scenarios,cbind(tmp_field_data[j,c("time_point","Farm","Group","Treatment")],
                                                                                                        tibble(cut_off = input_pars$cut_off[i]),
                                                                                                        pmajor(tmp_param_list,p = as.numeric(tmp_field_data[j,input_pars$cut_off[i]]), n = c(0,1) , average_by_p = TRUE)))
                                                      }
                                                    }
}

#arrange simulation outcomes for use in the report ####
#combine with the pmajor calculations
all_scenarios_sim<-(left_join(all_scenarios_sim,pmajor_scenarios%>%mutate(Treatment = str_replace(Treatment,"No vaccination","NOVAC"))%>%select(time_point,Farm,cut_off, Treatment,p), by = c("time_point","Farm","cut_off", "Treatment"), relationship = "many-to-many"))

#split data set
scenarios_sim <- all_scenarios_sim%>%filter(cut_off %in% c("O6","NOVAC"))
scenarios_sim_sensitivity <- all_scenarios_sim%>%filter(cut_off %in% c("O6","NOVAC"))

#get the detection results####
S1.results_all <- readRDS(list.files(output.dir.detection, pattern = "*_Detection_results_complete_EU*\\.RDS$", full.names = TRUE));
#split scenario in farm type and treatment
S1.results_all <- S1.results_all %>%
  mutate(farm = str_extract(scenario,c("Farm_A|Farm_B")),
         Treatment = str_extract(scenario,c("BI_|BI0|BIBOOSTER|CEVA|baseline")),
         cut_off = str_extract(scenario, "O5|O6|O7"),
         time_point = as.numeric(str_extract(scenario, paste0(unique(pmajor_scenarios$time_point),collapse = "|"))))

#remove underscore
S1.results_all$Treatment <- str_replace(S1.results_all$Treatment, "BI_","BI")
S1.results_all$Treatment <- replace_na(S1.results_all$Treatment, "NOVAC")
#replace NA by NOVAC
S1.results_all$cut_off <- replace_na(S1.results_all$cut_off,"NOVAC")

#add production phase
S1.results_all <- S1.results_all%>%mutate(production_phase = ifelse(time_point<120,"preproduction","production"))
S1.results_all$outbreak.size.ofconcern <- ifelse(S1.results_all$final.size > 0.001*46000, "LARGE","SMALL")

#join to have the proportion high-titre per scenario
S1.results_all<-(left_join(S1.results_all,pmajor_scenarios%>%mutate(farm = paste0("Farm_",Farm))%>%select(time_point,farm,cut_off, Treatment, p), by = c("time_point","farm","cut_off", "Treatment"), relationship = "many-to-many"))

#correct the half egg per day
S1.results_all$detect.Egg_I_shipped <-S1.results_all$detect.Egg_I_shipped /0.57  

# get the fraction of small outbreaks
small_outbreaks <- S1.results_all%>%reframe(.by = c(Treatment, production_phase),
                                            fraction_small = sum(final.size<= 0.001*46000)/n(),
                                            n_small = sum(final.size<= 0.001*46000)/n(),
                                            n_total  = n())


S1.results_all <- S1.results_all %>%
  mutate(det.since.intro.survey = ifelse(confirm.survey.result %in% c("NA", "NEGATIVE"), Inf, first.ac.detection.survey.time))
#define totally missed outbreaks
S1.results_all <- S1.results_all %>% mutate(
  missed = (is.infinite(det.since.intro.pas) & is.infinite(det.since.intro.ac) & is.infinite(det.since.intro.survey)) | (confirm.result == "NEGATIVE" & confirm.survey.result == "NEGATIVE" & live.confirmation.result == "NEGATIVE" & live.confirmation.survey.result == "NEGATIVE"),
  missed_nosurvey = (is.infinite(det.since.intro.pas) & is.infinite(det.since.intro.ac) ) | (confirm.result == "NEGATIVE" & confirm.survey.result == "NEGATIVE" & live.confirmation.result == "NEGATIVE"),
)
#add a variable determining the detection method
S1.results_all <- S1.results_all %>% 
  mutate(det.method = case_when(
    is.infinite(det.since.intro.pas) & is.infinite(det.since.intro.ac) & is.infinite(det.since.intro.survey)  ~ "none",
    is.infinite(det.since.intro.pas) & is.infinite(det.since.intro.ac)   ~ "none w/o survey",
    det.since.intro.ac <= det.since.intro.pas & det.since.intro.ac <= det.since.intro.survey & (confirm.result == "POSITIVE" | live.confirmation.result == "POSITIVE")~ "active",
    det.since.intro.pas <= det.since.intro.ac & det.since.intro.pas <=det.since.intro.survey ~ "passive",
    det.since.intro.survey <= det.since.intro.pas & det.since.intro.survey <= det.since.intro.ac & confirm.survey.result == "POSITIVE" ~ "survey",
    TRUE ~ "none"
  ))

S1.results_all<-S1.results_all%>%mutate(
  det.since.intro.ac.confirmed = ifelse(confirm.result == "POSITIVE"|live.confirmation.result == "POSITIVE", det.since.intro.ac, Inf),
  det.since.intro.survey.confirmed = ifelse(confirm.survey.result == "POSITIVE"|live.confirmation.survey.result == "POSITIVE", det.since.intro.survey, Inf))



#split the results in S1.results and S1.results.sensitivity
S1.results.sensitivity <- S1.results_all%>%filter(cut_off %in% c("O5","O7","NOVAC"))
S1.results <- S1.results_all%>%filter(cut_off %in% c("O6","NOVAC"))


#Input-tables as in report ####
input_table_beta_T_k_R <- with(all_files_selection,{data.frame(
  Treatment = Treatment,
  cut_off = replace_na(cut_off,"-"),
  cutoff_name = replace_na(str_replace(cut_off,pattern = "O", replacement = " &gt;= "),"&#45;"),
  beta_low_titre = beta_low_titre,
  Tinf_low_titre =Tinf_low_titre,
  k_Tinf_low_titre =k_Tinf_low_titre,
  R_low_titre = beta_low_titre * Tinf_low_titre,
  beta_low_titre = beta_high_titre,
  Tinf_low_titre =Tinf_high_titre,
  k_Tinf_low_titre =k_Tinf_high_titre,
  R_high_titre = beta_high_titre *  Tinf_high_titre
  
)})%>%mutate(cut_off = factor(cut_off,levels = c("-", "O6","O5", "O7")))%>%arrange(cut_off)
row.names(input_table_beta_T_k_R)<-NULL;


#################################################################################################
#  create tables and figure with the inputs in, the outputs of the simulations and calculations #
#################################################################################################


#CEVA
#Transmission parameters from transmission experiments used in this study.  All scenarios are run with the values for cut-off >=6.  
suppressMessages(knitr::kable(x= input_table_beta_T_k_R %>% filter(Treatment != "BI")%>%
                                filter(cut_off %in% c("-","O6"))%>%
                                select(!cutoff_name),
  escape= FALSE,col.names = c("Vaccination","cut-off", "$\\beta_{low titre}$","Mean inf. per","k","R_low","$\\beta_{high titre}$","Mean inf. per","k","R_high")))

#BI
#Transmission parameters from transmission experiments used in this study.  All scenarios are run with the values for cut-off >=6.  
suppressMessages(knitr::kable(x= input_table_beta_T_k_R %>% filter(Treatment != "CEVA")%>%
                                filter(cut_off %in% c("-","O6"))%>%
                                select(!cutoff_name),
                              escape= FALSE,col.names = c("Vaccination","cut-off", "$\\beta_{low titre}$","Mean inf. per","k","R_low","$\\beta_{high titre}$","Mean inf. per","k","R_high")))



#### Proportion high titre ####

#times <- seq( 0,588,1) # remove redundatn?

plot_field_data <- data_field_study_cut_offs%>%
  filter(!str_detect(Treatment, "Neg. controle"))%>%
  pivot_longer(cols = starts_with("O"), names_to = "cut_off")%>%
  mutate(cut_off = factor(cut_off, levels = c("O6","O5","O7")))%>%arrange(cut_off)

##### CEVA

# "Proportion birds with high titre (log10(HI-titre)>= 5, 6 or 7) after vaccination (in days). Lines are fitted curves, dots are data points. ####
ggplot(data = data.frame(plot_field_data%>%filter(time_point>=120)%>%filter(cut_off == "O6" & !str_detect(Treatment ,"BI"))))+
  geom_point(aes(x = time_point, y = value, colour = Treatment))+
  scale_colour_discrete(labels = ceva_labels)+
  scale_y_continuous( labels = function(x) paste0(round(x * 100), "%"))+
  labs(x = "Time since prime vaccination (days)", y = "Percentage chickens with high titre ")+
  theme_minimal()+
  theme(legend.title = element_blank())

suppressMessages(ggsave("./figures/buildup_ceva_6.png"))


#### BI ####
#"Proportion birds with high titre (log10(HI-titre)>= 5, 6 or 7) after vaccination (in days). Lines are fitted curves, dots are data points."}

ggplot(data = data.frame(plot_field_data%>%filter(time_point>=120)%>%filter(cut_off == "O6" & !str_detect(Treatment ,"CEVA"))))+
  geom_point(aes(x = time_point, y = value, colour = Treatment))+
  scale_colour_discrete(labels = bi_labels)+
  scale_y_continuous( labels = function(x) paste0(round(x * 100), "%"))+
  labs(x = "Time since prime vaccination (days)", y = "Percentage chickens with high titre ")+
  theme_minimal()+
  theme(legend.title = element_blank())


suppressMessages(ggsave("./figures/buildup_bi_6.png"))


### Parameters for consequences
#"Egg production curve with in black normal curve and red line non-disfigured alive eggs produced by infected birds."}
#Egg : Flock production function governing percentage of birds laying eggs per day at time t
flock_production <- function(t, age_at_d0, eh) {
      t_shift <- (t + age_at_d0)/7 #function assumes week 0 = birth, model uses day 0 = first day of maturity
      a <- 103
      b <- 0.0016
      c <- 1.16
      d <- 20.75
      eh*((a * exp(-b * t_shift)) / (1 + exp(-c * (t_shift - d))))/100
    }

ggplot(aes(t), data = data.frame(t = seq(0, 19*30,0.1)))+
  geom_function(fun = flock_production, args = list(age_at_d0 =120, eh =0.57 ))+
  geom_function(fun = flock_production, args = list(age_at_d0 =120, eh =0.285*(1-0.1)), colour = "red" )+
  xlab("Time since arrival at layer farm i.e. age + 120 (days)")+
  ylab("Eggs per day per layer")#+
  #ggtitle("Egg production curve.  ")

suppressMessages(ggsave("./figures/eggproductioncurve.png"))



# Results ####
## Reproduction number and critical high-titre proportion
pmajor_scenarios_crit <- pmajor_scenarios%>%
  mutate(production_phase = ifelse(time_point<120,"preproduction","production"),
  above_pc = as.numeric(p>pc))

average_pmajor_scenarios_crit<- pmajor_scenarios_crit%>%filter(cut_off == "O6")%>%reframe(.by = c(Treatment, production_phase),
                           prop_above = sum(above_pc)/n(),
                           prop_above_ll = clopper_pearson( sum(above_pc),n())[1],
                           prop_above_ul  = clopper_pearson( sum(above_pc),n())[2],
                           n_above = sum(above_pc),
                           n_total = n()  )%>%
  pivot_wider(names_from = production_phase, 
              values_from = c(prop_above,prop_above_ll,prop_above_ul, n_above, n_total))

reduction_R <- pmajor_scenarios_crit%>%filter(cut_off == "O6")%>%reframe(.by = c(Treatment, production_phase, above_pc),
                                                                         meanR0 = mean(R0),
                                                                         meanRv = mean(Rv),
                                                                         mean_reduction = 100*mean(1 -(Rv/R0)),
                                                                         min_reduction = 100*min(1- (Rv/R0)),
                                                                         max_reduction = 100*max(1 - (Rv/R0)),
                                                                         sample_times = n())


### CEVA

#"The reproduction number at sampling moments of the field study for a cut-off value of log 6 per farm. Dashed line indicates the line R = 1, which is the threshold for herd-immunity. "}
ggplot(pmajor_scenarios%>%filter(time_point>=120)%>%filter(!str_detect(Treatment,"BI"))%>%filter(cut_off %in% c("O6","NOVAC")))+
  geom_point(aes(x = time_point, y = Rv, colour = Treatment, shape = Farm))+
  scale_y_continuous(transform = "log10",breaks = c(0.5,1,1.5,2,3,5,10), limits = c(0.5,10))+
  scale_colour_discrete(labels = ceva_labels)+
  scale_shape_discrete(labels = list('A' = "Farm A", 'B' = "Farm B"))+
  labs(x = "Sample time (days)", y = expression(italic(R[v])))+
  geom_hline(aes(yintercept =1), linetype = "dashed")+
  theme_minimal()+
  theme(legend.title = element_blank())
  

suppressMessages(ggsave("./figures/R_ceva.png"))


average_pmajor_scenarios_crit%>%filter(!str_detect(Treatment,"BI"))%>%
  select(1,3,5,7) %>%
  knitr::kable(digits =2, caption = "Fraction and 95%-confidence interval of time points with a proportion high titre above the critical fraction  during production phase (age above 120 days).", col.names = c("Vaccination","production", "2.5%", "97.5%"))

reduction_R%>%filter(!str_detect(Treatment,"BI") & production_phase == "production")%>%
  #select(1,3,5,7) %>%
  knitr::kable(digits =2, 
               caption = "Reduction of R.")
# Proportion of runs resulting in negligible outbreaks
small_outbreaks%>%filter(!str_detect(Treatment,"BI"))%>%
  knitr::kable(digits =2, 
               caption = "Proportion of introduction with negligible spread (i.e. <= 0.1% of farm).")


### BI

# "The reproduction number at sampling moments of the field study for a cut-off value of log 6 per farm. Dashed line indicates the line R = 1, which is the threshold for herd-immunity. "}
ggplot(pmajor_scenarios%>%filter(time_point>=120)%>%filter(!str_detect(Treatment,"CEVA"))%>%filter(cut_off %in% c("O6","NOVAC")))+
  geom_point(aes(x = time_point, y = Rv, colour = Treatment, shape = Farm))+
  scale_y_continuous(transform = "log10",breaks = c(0.5,1,1.5,2,3,5,10), limits = c(0.5,10))+
  scale_colour_discrete(labels = bi_labels)+
  scale_shape_discrete(labels = list('A' = "Farm A", 'B' = "Farm B"))+
  labs(x = "Sample time (days)", y = expression(italic(R[v])))+
  geom_hline(aes(yintercept =1), linetype = "dashed")+
  theme_minimal()+
  theme(legend.title = element_blank())


suppressMessages(ggsave("./figures/R_BI.png"))

average_pmajor_scenarios_crit%>%filter(!str_detect(Treatment,"CEVA"))%>%
  select(1,3,5,7) %>%
  knitr::kable(digits =2, caption = "Fraction and 95%-confidence interval of time points with a proportion high titre above the critical fraction  during production phase (age above 120 days).", col.names = c("Vaccination","production", "2.5%", "97.5%"))

reduction_R%>%filter(!str_detect(Treatment,"CEVA") & production_phase == "production")%>%
  #select(1,3,5,7) %>%
  knitr::kable(digits =2, 
               caption = "Reduction of R.")

# Proportion of runs resulting in negligible outbreaks
small_outbreaks%>%filter(!str_detect(Treatment,"CEVA"))%>%
  knitr::kable(digits =2, 
               caption = "Proportion of introduction with negligible spread (i.e. <= 0.1% of farm).")





## Detection results ####
#choose interval   
perc_low = 0.25;
perc_high = 0.75;
### Fraction large and small outbreaks detected

#create variable to determine detection  by survey 

table_data <- S1.results%>%reframe(.by = c(production_phase,Treatment,outbreak.size.ofconcern,time.interval.ac),
  passive = sum(as.numeric(!is.infinite(det.since.intro.pas)))/n(),
  active  = sum(as.numeric(confirm.result == "POSITIVE" | live.confirmation.result == "POSITIVE" )*as.numeric(!is.infinite(det.since.intro.ac)))/n(),
  survey = sum(as.numeric(confirm.survey.result == "POSITIVE" | live.confirmation.survey.result == "POSITIVE")*as.numeric(!is.infinite(det.since.intro.survey)))/n(),
  missed = 1 - sum(missed_nosurvey)/n(),
 # missed_nosurvey = missed_nosurvey/n(),
  shortest_passive = sum((det.since.intro.pas<det.since.intro.ac.confirmed) & (det.since.intro.pas<det.since.intro.survey.confirmed), na.rm = TRUE)/n(),
  shortest_active = sum((det.since.intro.pas>det.since.intro.ac.confirmed) & (det.since.intro.ac<det.since.intro.survey.confirmed), na.rm = TRUE)/n(),
  shortest_survey = sum((det.since.intro.pas>det.since.intro.survey.confirmed) & (det.since.intro.ac>det.since.intro.survey.confirmed), na.rm = TRUE)/n()
  )%>%pivot_wider(names_from = c(outbreak.size.ofconcern),
                  values_from = c(active, 
                                  passive, 
                                  survey, 
                                  missed,
                                  shortest_passive,
                                  shortest_active,
                                  shortest_survey))

detection_times_vs_interval <- S1.results%>%filter(outbreak.size.ofconcern == "LARGE" & is.finite(min.det.time))%>%reframe(.by = c(Treatment, production_phase, time.interval.ac),                                                    
                                                    p1_detection_time = quantile(min.det.time,perc_low),
                                                    median_detection_time = median(min.det.time),
                                                    p2_detection_time = quantile(min.det.time,perc_high)
                                                      )

correct_table_data <- table_data
for(j in unique(table_data$production_phase))
  {
    for(k in unique(table_data$Treatment))
    {
              correct_table_data[is.finite(correct_table_data$time.interval.ac) & correct_table_data$production_phase == j & correct_table_data$Treatment == k, "survey_SMALL"] =correct_table_data[is.infinite(correct_table_data$time.interval.ac) & correct_table_data$production_phase == j & correct_table_data$Treatment == k, c("survey_SMALL")] ;
        correct_table_data[is.finite(correct_table_data$time.interval.ac) & correct_table_data$production_phase == j & correct_table_data$Treatment == k, "survey_LARGE"] =correct_table_data[is.infinite(correct_table_data$time.interval.ac) & correct_table_data$production_phase == j & correct_table_data$Treatment == k, c("survey_LARGE")] 
     
    }
  }


#size detected v.s  undetected small outbreaks
size_det_vs_undet <-S1.results%>%reframe(.by = c(production_phase,Treatment,outbreak.size.ofconcern,time.interval.ac,missed),
                                         min_outbreak_size = min(detect.N.dead.detectables + detect.N.live.detectables), #number of dead and live infected at detection.
                                         median_outbreak_size = median(detect.N.dead.detectables + detect.N.live.detectables),
                                         max_outbreak_size = quantile(detect.N.dead.detectables+ detect.N.live.detectables,perc_high))%>%
  arrange(outbreak.size.ofconcern,Treatment, production_phase, time.interval.ac)


#create range of detection times without vaccination
det_interval_novac <- detection_times_vs_interval%>%filter(Treatment == "NOVAC" & is.infinite(time.interval.ac ))


#consequences
consequences <- S1.results%>%filter(production_phase == "production")%>%reframe(
  .by = c(Treatment,  time.interval.ac),
  fraction_runs_with_transport = sum(detect.Egg_I_shipped>=1,na.rm = TRUE)/n(),
  median_I_eggs = median(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),na.rm = TRUE), 
  perc_low_I_eggs = quantile(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),perc_low,na.rm = TRUE),
  perc_high_I_eggs = quantile(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),perc_high, na.rm = TRUE),
  median_I_dead = median(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),na.rm = TRUE), 
  perc_low_I_dead = quantile(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),perc_low,na.rm = TRUE),
  perc_high_I_dead = quantile(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),perc_high,na.rm = TRUE),
  median_I_alive = median(round(detect.N.live.detectables[detect.N.live.detectables>=1]),na.rm = TRUE), 
  perc_low_I_alive = quantile(round(detect.N.live.detectables[detect.N.live.detectables>=1]),perc_low,na.rm = TRUE),
  perc_high_I_alive = quantile(round(detect.N.live.detectables[detect.N.live.detectables>=1]),perc_high,na.rm = TRUE)
  )

#combine outcomes to useable table.

output_within_table <- full_join(full_join(correct_table_data%>%filter(production_phase == "production")%>%select(Treatment, time.interval.ac, missed_LARGE), 
          detection_times_vs_interval%>%filter(production_phase == "production")%>%select(Treatment, time.interval.ac,p1_detection_time, median_detection_time, p2_detection_time)),
          consequences)




### CEVA
ceva_tab_output <- output_within_table%>%
  filter(!str_detect(Treatment, "BI"))

ceva_tab_output%>%knitr::kable(digits = 2, col.namess = c("Treatment","interval","P-det","l_det_time","med_det_time","h_det_time",
                                                                                 "frac.with","med I", "l I","h I","med dead","l_dead","h-dead","med live", "l live","h-live"))

sink("output/consequences_ceva.html")
cat(ceva_tab_output%>%knitr::kable(digits = 2, format = "html", col.namess = c("Treatment","interval","P-det","l_det_time","med_det_time","h_det_time",
                                                                               "frac.with","med I", "l I","h I","med dead","l_dead","h-dead","med live", "l live","h-live"))
)
sink()
### BI
bi_tab_output <- output_within_table%>%
  filter(!str_detect(Treatment, "CEVA"))

bi_tab_output%>%knitr::kable(digits = 2, col.namess = c("Treatment","interval","P-det","l_det_time","med_det_time","h_det_time",
                                                                                 "frac.with","med I", "l I","h I","med dead","l_dead","h-dead","med live", "l live","h-live"))
sink("output/consequences_bi.html")
cat(bi_tab_output%>%knitr::kable(digits = 2, format = "html", col.namess = c("Treatment","interval","P-det","l_det_time","med_det_time","h_det_time",
                                                                               "frac.with","med I", "l I","h I","med dead","l_dead","h-dead","med live", "l live","h-live"))
)
sink()


#input into introduction and between-farm model
#probability of missing a vaccinated farm
p_miss_farm <- table_data%>%reframe(.by = c(Treatment, time.interval.ac),
                                    p_miss = mean(missed_LARGE))

#infectious period detected farms with LARGE outbreak based on average of CEVA and BI BOOSTER 
fit_data <- S1.results%>%
  filter(outbreak.size.ofconcern == "LARGE" & is.finite(min.det.time) & Treatment != "BI")%>%
  reframe(.by = c(time.interval.ac,Treatment),
          vaccinated = c("Unvaccinated","Vaccinated")[(as.numeric(Treatment != "NOVAC")+1)],
          min.det.time)
                                                                                            
gamma_distribution <- fit_data%>%reframe(.by = c(time.interval.ac, vaccinated),
                                         alpha = mean(min.det.time)^2 / var(min.det.time),
                                         beta = mean(min.det.time) / var(min.det.time))
names(gamma_distribution)<- c("time.interval.ac", "Treatment",       "alpha",            "beta")
gamma_distribution <- rbind(gamma_distribution,
                            fit_data%>%reframe(.by = c(time.interval.ac, Treatment),
                                         alpha = mean(min.det.time)^2 / var(min.det.time),
                                         beta = mean(min.det.time) / var(min.det.time)))


#for ease of plotting create values of fitted functions
gamma_distribution_to_plot <- gamma_distribution %>%
  mutate(
    min.det.time = list(0:max(fit_data$min.det.time)),
    y = map2(alpha, beta, ~ dgamma(0:max(fit_data$min.det.time), .x, .y))
  ) %>%
  unnest(cols = c(min.det.time, y))%>%select(Treatment, time.interval.ac, min.det.time, y)


ggplot() +
  #Histogram
  geom_histogram(
    data = fit_data,
    aes(x = min.det.time, y = after_stat(density), fill = as.factor(time.interval.ac)),
    bins = 30,
    alpha = 0.6 ) +
  #Gamma density
  geom_path(
    data =gamma_distribution_to_plot,
    aes(x = min.det.time, y = y, colour = as.factor(time.interval.ac)),
    linewidth = 1 ) +
  # Facets
  facet_grid(rows = vars(Treatment)) +
  # Labels and theme
  labs(
    title = "Detection Time Distribution by Vaccination Status",
    x = "Detection Time",
    y = "Density"
  ) 


# "Parameters for the gammadistribution of detection times"}
knitr::kable(gamma_distribution%>%filter(Treatment == "Vaccinated")%>%mutate(mean = alpha/beta,
                            var = alpha/beta^2), digits = 2)

#"Parameters for the gammadistribution of detection times"}
knitr::kable(gamma_distribution%>%filter(str_detect(Treatment, "CEVA"))%>%mutate(mean = alpha/beta,
                            var = alpha/beta^2), digits = 2)

#"Parameters for the gammadistribution of detection times"}
knitr::kable(gamma_distribution%>%filter(str_detect(Treatment, "BI"))%>%mutate(mean = alpha/beta,
                            var = alpha/beta^2), digits = 2)

#################################################################################################
# SENSITIVITY ANALYSES 
# create tables and figure with the inputs in, the outputs of the simulations and calculations #
#################################################################################################


#CEVA for sensitivity
suppressMessages(knitr::kable(x= input_table_beta_T_k_R %>% filter(Treatment != "BI")%>%
                                #  filter(cut_off %in% c("-","O6"))%>%
                                select(!cutoff_name),
                              escape= FALSE,col.names = c("Vaccination","cut-off", "$\\beta_{low titre}$","Mean inf. per","k","R_low","$\\beta_{high titre}$","Mean inf. per","k","R_high")))


#BI for sensitivity
suppressMessages(knitr::kable(x= input_table_beta_T_k_R %>% filter(Treatment != "CEVA")%>%
                                #  filter(cut_off %in% c("-","O6"))%>%
                                select(!cutoff_name),
                              escape= FALSE,col.names = c("Vaccination","cut-off", "$\\beta_{low titre}$","Mean inf. per","k","R_low","$\\beta_{high titre}$","Mean inf. per","k","R_high")))


#### Proportion high titre ####
plot_field_data <- data_field_study_cut_offs%>%
  filter(!str_detect(Treatment, "Neg. controle"))%>%
  pivot_longer(cols = starts_with("O"), names_to = "cut_off")%>%
  mutate(cut_off = factor(cut_off, levels = c("O6","O5","O7")))%>%arrange(cut_off)

##### CEVA

# "Proportion birds with high titre (log10(HI-titre)>= 5, 6 or 7) after vaccination (in days). Lines are fitted curves, dots are data points. ####
ggplot(data = data.frame(plot_field_data%>%filter(time_point>=120)%>%filter(str_detect(Treatment ,"BI"))))+
  geom_point(aes(x = time_point, y = value, colour = Treatment))+
  scale_colour_discrete(labels = ceva_labels)+
  scale_y_continuous( labels = function(x) paste0(round(x * 100), "%"))+
  labs(x = "Time since prime vaccination (days)", y = "Percentage chickens with high titre ")+
  theme_minimal()+
  theme(legend.title = element_blank())

suppressMessages(ggsave("./figures/buildup_ceva_sens.png"))


#### BI ####
#"Proportion birds with high titre (log10(HI-titre)>= 5, 6 or 7) after vaccination (in days). Lines are fitted curves, dots are data points."}

ggplot(data = data.frame(plot_field_data%>%filter(time_point>=120)%>%filter(!str_detect(Treatment ,"CEVA"))))+
  geom_point(aes(x = time_point, y = value, colour = Treatment))+
  scale_colour_discrete(labels = bi_labels)+
  scale_y_continuous( labels = function(x) paste0(round(x * 100), "%"))+
  labs(x = "Time since prime vaccination (days)", y = "Percentage chickens with high titre ")+
  theme_minimal()+
  theme(legend.title = element_blank())


suppressMessages(ggsave("./figures/buildup_bi_sens.png"))


# Results ####

## Reproduction number and critical high-titre proportion
average_pmajor_scenarios_crit<- pmajor_scenarios_crit%>%reframe(.by = c(Treatment, production_phase,cut_off),
                                                                                          prop_above = sum(above_pc)/n(),
                                                                                          prop_above_ll = clopper_pearson( sum(above_pc),n())[1],
                                                                                          prop_above_ul  = clopper_pearson( sum(above_pc),n())[2],
                                                                                          n_above = sum(above_pc),
                                                                                          n_total = n()  )%>%
  pivot_wider(names_from = production_phase, values_from = c(prop_above,prop_above_ll,prop_above_ul, n_above, n_total))

reduction_R <- pmajor_scenarios_crit%>%reframe(.by = c(Treatment, cut_off,production_phase, above_pc),
                                                                         meanR0 = mean(R0),
                                                                         meanRv = mean(Rv),
                                                                         mean_reduction = 100*mean(1 -(Rv/R0)),
                                                                         min_reduction = 100*min(1- (Rv/R0)),
                                                                         max_reduction = 100*max(1 - (Rv/R0)),
                                                                         sample_times = n())


### CEVA

#"The reproduction number at sampling moments of the field study for a cut-off value of log 6 per farm. Dashed line indicates the line R = 1, which is the threshold for herd-immunity. "}
ggplot(pmajor_scenarios%>%filter(time_point>=120)%>%
         filter(!str_detect(Treatment,"BI")))+
  geom_point(aes(x = time_point, y = Rv, colour = Treatment, shape = Farm))+
  scale_y_continuous(transform = "log10",breaks = c(0.5,1,1.5,2,3,5,10), limits = c(0.5,10))+
  scale_colour_discrete(labels = ceva_labels)+
  scale_shape_discrete(labels = list('A' = "Farm A", 'B' = "Farm B"))+
  labs(x = "Sample time (days)", y = expression(italic(R[v])))+
  geom_hline(aes(yintercept =1), linetype = "dashed")+
  theme_minimal()+
  theme(legend.title = element_blank())+facet_grid(cut_off~.)


suppressMessages(ggsave("./figures/R_ceva_sens.png"))


average_pmajor_scenarios_crit%>%filter(!str_detect(Treatment,"BI"))%>%
  select(1,2,3,5,7) %>%
  knitr::kable(digits =2, caption = "Fraction and 95%-confidence interval of time points with a proportion high titre above the critical fraction  during production phase (age above 120 days).", col.names = c("Vaccination","cutoff","production", "2.5%", "97.5%"))

reduction_R%>%filter(!str_detect(Treatment,"BI") & production_phase == "production")%>%
  #select(1,3,5,7) %>%
  knitr::kable(digits =2, 
               caption = "Reduction of R.")
# # Proportion of runs resulting in negligible outbreaks
# small_outbreaks%>%filter(!str_detect(Treatment,"BI"))%>%
#   knitr::kable(digits =2, 
#                caption = "Proportion of introduction with negligible spread (i.e. <= 0.1% of farm).")


### BI

# "The reproduction number at sampling moments of the field study for a cut-off value of log 6 per farm. Dashed line indicates the line R = 1, which is the threshold for herd-immunity. "}
ggplot(pmajor_scenarios%>%filter(time_point>=120)%>%filter(!str_detect(Treatment,"CEVA")))+
  geom_point(aes(x = time_point, y = Rv, colour = Treatment, shape = Farm))+
  scale_y_continuous(transform = "log10",breaks = c(0.5,1,1.5,2,3,5,10), limits = c(0.5,10))+
  scale_colour_discrete(labels = bi_labels)+
  scale_shape_discrete(labels = list('A' = "Farm A", 'B' = "Farm B"))+
  labs(x = "Sample time (days)", y = expression(italic(R[v])))+
  geom_hline(aes(yintercept =1), linetype = "dashed")+
  theme_minimal()+
  theme(legend.title = element_blank())+facet_grid(cut_off~.)


suppressMessages(ggsave("./figures/R_BI_sens.png"))

average_pmajor_scenarios_crit%>%filter(!str_detect(Treatment,"CEVA"))%>%
  select(1,2,3,5,7) %>%
  knitr::kable(digits =2, caption = "Fraction and 95%-confidence interval of time points with a proportion high titre above the critical fraction  during production phase (age above 120 days).", col.names = c("Vaccination","cut-off","production", "2.5%", "97.5%"))

reduction_R%>%filter(!str_detect(Treatment,"CEVA") & production_phase == "production")%>%
  #select(1,3,5,7) %>%
  knitr::kable(digits =2, 
               caption = "Reduction of R.")

# # Proportion of runs resulting in negligible outbreaks
# small_outbreaks%>%filter(!str_detect(Treatment,"CEVA"))%>%
#   knitr::kable(digits =2, 
#                caption = "Proportion of introduction with negligible spread (i.e. <= 0.1% of farm).")





## Detection results ####
#choose interval   
perc_low = 0.25;
perc_high = 0.75;
### Fraction large and small outbreaks detected

#create variable to determine detection  by survey 

table_data <- S1.results.sensitivity%>%reframe(.by = c(production_phase,cut_off, Treatment,outbreak.size.ofconcern,time.interval.ac),
                                   passive = sum(as.numeric(!is.infinite(det.since.intro.pas)))/n(),
                                   active  = sum(as.numeric(confirm.result == "POSITIVE" | live.confirmation.result == "POSITIVE" )*as.numeric(!is.infinite(det.since.intro.ac)))/n(),
                                   survey = sum(as.numeric(confirm.survey.result == "POSITIVE" | live.confirmation.survey.result == "POSITIVE")*as.numeric(!is.infinite(det.since.intro.survey)))/n(),
                                   missed = 1 - sum(missed_nosurvey)/n(),
                                   # missed_nosurvey = missed_nosurvey/n(),
                                   shortest_passive = sum((det.since.intro.pas<det.since.intro.ac.confirmed) & (det.since.intro.pas<det.since.intro.survey.confirmed), na.rm = TRUE)/n(),
                                   shortest_active = sum((det.since.intro.pas>det.since.intro.ac.confirmed) & (det.since.intro.ac<det.since.intro.survey.confirmed), na.rm = TRUE)/n(),
                                   shortest_survey = sum((det.since.intro.pas>det.since.intro.survey.confirmed) & (det.since.intro.ac>det.since.intro.survey.confirmed), na.rm = TRUE)/n()
)%>%pivot_wider(names_from = c(outbreak.size.ofconcern),
                values_from = c(active, 
                                passive, 
                                survey, 
                                missed,
                                shortest_passive,
                                shortest_active,
                                shortest_survey))

detection_times_vs_interval <- S1.results.sensitivity%>%
  filter(outbreak.size.ofconcern == "LARGE" & is.finite(min.det.time))%>%
  reframe(.by = c(Treatment,cut_off, production_phase, time.interval.ac),                                                    
                                                                                                                           p1_detection_time = quantile(min.det.time,perc_low),
                                                                                                                           median_detection_time = median(min.det.time),
                                                                                                                           p2_detection_time = quantile(min.det.time,perc_high)
)

correct_table_data <- table_data
for(j in unique(table_data$production_phase))
{
  for(k in unique(table_data$Treatment))
  {
    for(l in unique(table_data$cut_off))
    {
    correct_table_data[is.finite(correct_table_data$time.interval.ac) & correct_table_data$production_phase == j & correct_table_data$Treatment == k & correct_table_data$cut_off ==l, "survey_SMALL"] =correct_table_data[is.infinite(correct_table_data$time.interval.ac) & correct_table_data$production_phase == j & correct_table_data$Treatment == k& correct_table_data$cut_off ==l, c("survey_SMALL")] ;
    correct_table_data[is.finite(correct_table_data$time.interval.ac) & correct_table_data$production_phase == j & correct_table_data$Treatment == k& correct_table_data$cut_off ==l, "survey_LARGE"] =correct_table_data[is.infinite(correct_table_data$time.interval.ac) & correct_table_data$production_phase == j & correct_table_data$Treatment == k& correct_table_data$cut_off ==l, c("survey_LARGE")] 
    }
  }
}


#size detected v.s  undetected small outbreaks
size_det_vs_undet <-S1.results.sensitivity%>%reframe(.by = c(production_phase,cut_off,Treatment,outbreak.size.ofconcern,time.interval.ac,missed),
                                         min_outbreak_size = min(detect.N.dead.detectables + detect.N.live.detectables), #number of dead and live infected at detection.
                                         median_outbreak_size = median(detect.N.dead.detectables + detect.N.live.detectables),
                                         max_outbreak_size = quantile(detect.N.dead.detectables+ detect.N.live.detectables,perc_high))%>%
  arrange(outbreak.size.ofconcern,Treatment, production_phase, time.interval.ac)


#create range of detection times without vaccination
det_interval_novac <- detection_times_vs_interval%>%filter(Treatment == "NOVAC" & is.infinite(time.interval.ac ))


#consequences
consequences <- S1.results.sensitivity%>%filter(production_phase == "production")%>%reframe(
  .by = c(Treatment, cut_off, time.interval.ac),
  fraction_runs_with_transport = sum(detect.Egg_I_shipped>=1,na.rm = TRUE)/n(),
  median_I_eggs = median(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),na.rm = TRUE), 
  perc_low_I_eggs = quantile(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),perc_low,na.rm = TRUE),
  perc_high_I_eggs = quantile(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),perc_high, na.rm = TRUE),
  median_I_dead = median(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),na.rm = TRUE), 
  perc_low_I_dead = quantile(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),perc_low,na.rm = TRUE),
  perc_high_I_dead = quantile(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),perc_high,na.rm = TRUE),
  median_I_alive = median(round(detect.N.live.detectables[detect.N.live.detectables>=1]),na.rm = TRUE), 
  perc_low_I_alive = quantile(round(detect.N.live.detectables[detect.N.live.detectables>=1]),perc_low,na.rm = TRUE),
  perc_high_I_alive = quantile(round(detect.N.live.detectables[detect.N.live.detectables>=1]),perc_high,na.rm = TRUE)
)

#combine outcomes to useable table.

output_within_table <- full_join(full_join(correct_table_data%>%filter(production_phase == "production")%>%select(Treatment,cut_off, time.interval.ac, missed_LARGE), 
                                           detection_times_vs_interval%>%filter(production_phase == "production")%>%select(Treatment,cut_off, time.interval.ac,p1_detection_time, median_detection_time, p2_detection_time)),
                                 consequences)




### CEVA
ceva_tab_output <- output_within_table%>%
  filter(!str_detect(Treatment, "BI"))

ceva_tab_output%>%knitr::kable(digits = 2, col.namess = c("Treatment","cut_off","interval","P-det","l_det_time","med_det_time","h_det_time",
                                                          "frac.with","med I", "l I","h I","med dead","l_dead","h-dead","med live", "l live","h-live"))

sink("output/consequences_ceva_sensitivity.html")
cat(ceva_tab_output%>%knitr::kable(digits = 2, format = "html", col.namess = c("Treatment","cut_off","interval","P-det","l_det_time","med_det_time","h_det_time",
                                                                               "frac.with","med I", "l I","h I","med dead","l_dead","h-dead","med live", "l live","h-live"))
)
sink()
### BI
bi_tab_output <- output_within_table%>%
  filter(!str_detect(Treatment, "CEVA"))

bi_tab_output%>%knitr::kable(digits = 2, col.namess = c("Treatment","cut_off","interval","P-det","l_det_time","med_det_time","h_det_time",
                                                        "frac.with","med I", "l I","h I","med dead","l_dead","h-dead","med live", "l live","h-live"))
sink("output/consequences_bi_sensitivity.html")
cat(bi_tab_output%>%knitr::kable(digits = 2, format = "html", col.namess = c("Treatment","cut_off","interval","P-det","l_det_time","med_det_time","h_det_time",
                                                                             "frac.with","med I", "l I","h I","med dead","l_dead","h-dead","med live", "l live","h-live"))
)
sink()


# #input into introduction and between-farm model
# #probability of missing a vaccinated farm
# p_miss_farm <- table_data%>%reframe(.by = c(Treatment, cut_off,time.interval.ac),
#                                     p_miss = mean(missed_LARGE))
# 
# #infectious period detected farms with LARGE outbreak based on average of CEVA and BI BOOSTER 
# fit_data <- S1.results.sensitivity%>%
#   filter(outbreak.size.ofconcern == "LARGE" & is.finite(min.det.time) & Treatment != "BI")%>%
#   reframe(.by = c(cut_off,time.interval.ac,Treatment),
#           vaccinated = c("Unvaccinated","Vaccinated")[(as.numeric(Treatment != "NOVAC")+1)],
#           min.det.time)
# 
# gamma_distribution <- fit_data%>%reframe(.by = c(cut_off,time.interval.ac, vaccinated),
#                                          alpha = mean(min.det.time)^2 / var(min.det.time),
#                                          beta = mean(min.det.time) / var(min.det.time))
# names(gamma_distribution)<- c("cut_off","time.interval.ac", "Treatment",       "alpha",            "beta")
# gamma_distribution <- rbind(gamma_distribution,
#                             fit_data%>%reframe(.by = c(cut_off,time.interval.ac, Treatment),
#                                                alpha = mean(min.det.time)^2 / var(min.det.time),
#                                                beta = mean(min.det.time) / var(min.det.time)))
# 
# 
# #for ease of plotting create values of fitted functions
# gamma_distribution_to_plot <- gamma_distribution %>%
#   mutate(
#     min.det.time = list(0:max(fit_data$min.det.time)),
#     y = map2(alpha, beta, ~ dgamma(0:max(fit_data$min.det.time), .x, .y))
#   ) %>%
#   unnest(cols = c(min.det.time, y))%>%select(Treatment, cut_off,time.interval.ac, min.det.time, y)
# 
# 
# ggplot() +
#   #Histogram
#   geom_histogram(
#     data = fit_data,
#     aes(x = min.det.time, y = after_stat(density), fill = as.factor(time.interval.ac)),
#     bins = 30,
#     alpha = 0.6 ) +
#   #Gamma density
#   geom_path(
#     data =gamma_distribution_to_plot,
#     aes(x = min.det.time, y = y, colour = as.factor(time.interval.ac)),
#     linewidth = 1 ) +
#   # Facets
#   facet_grid(rows = vars(Treatment)) +
#   # Labels and theme
#   labs(
#     title = "Detection Time Distribution by Vaccination Status",
#     x = "Detection Time",
#     y = "Density"
#   ) 
# 
# 
# # "Parameters for the gammadistribution of detection times"}
# knitr::kable(gamma_distribution%>%filter(Treatment == "Vaccinated")%>%mutate(mean = alpha/beta,
#                                                                              var = alpha/beta^2), digits = 2)
# 
# #"Parameters for the gammadistribution of detection times"}
# knitr::kable(gamma_distribution%>%filter(str_detect(Treatment, "CEVA"))%>%mutate(mean = alpha/beta,
#                                                                                  var = alpha/beta^2), digits = 2)
# 
# #"Parameters for the gammadistribution of detection times"}
# knitr::kable(gamma_distribution%>%filter(str_detect(Treatment, "BI"))%>%mutate(mean = alpha/beta,
#                                                                                var = alpha/beta^2), digits = 2)


