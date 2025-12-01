### "Surveillance in a vaccinated poultry population" ####
# Author e.a.j.fischer@uu.nl 

setwd("C:/Users/fisch106/SurfDrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination - Copy (2)")
  
# Source required scripts that load and process the data ####
  
## libraries and functions ####
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

#set those sections to run. Setting these to false will result in using previous runs of these sections or scripts. This speeds up when processing the data for editing text and figures etc. 
run_detection_module = FALSE
run_process_field_experiment = TRUE
run_combine_daily_data = FALSE

#set directory where outputs are stored and run scripts to load required libraries and calculation of probability of a major outbreak
hpc_directory <- "C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_output"
suppressPackageStartupMessages(source("./src/loadLibraries.R"))
suppressMessages(source("./src/pmaj_calculation.R"))


#labels
ceva_labels <- c('CEVA' = paste0("VECTORMUNE",reg_symbol),
                 'NOVAC' = "No vaccination",
                 'Neg_controle CEVA' = "Control group")
bi_labels <- c('BI' = "VAXXITEK HVT + IBD + H5",
               'BIBOOSTER' = paste0("VAXXITEK HVT + IBD + H5 \n + Volvac",reg_symbol," B.E.S.T. AI + ND"),
               'NOVAC' = "No vaccination",
               'Neg_controle BI' = "Control group")

## run analyses or load analysis ####
if(run_process_field_experiment){
  suppressMessages(source("./src/ProcessFieldExperiment.R"))
}
if(run_detection_module){ #change to FALSE to prevent running this code
  #run analysis
  suppressMessages(source("./src/pre_processing_detModule.R"))
  suppressMessages(source("./src/detectionModule.R"))
}


## load the simulations ####
  
#choose to load previously combined data set or not
output.dir <- "./HPC_output/DailyData/";

if(!run_combine_daily_data){
  # Load the data
  file.list <- list.files(output.dir, pattern = "DailyData.*\\.RDS$", full.names = TRUE);
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
  saveRDS(all_scenarios_sim,file = paste0(output.dir, "/Combined/Combined_DailyData.RDS"))
}else
{
  # Load the data
  all_scenarios_sim <-readRDS(file = paste0(output.dir, "/Combined/Combined_DailyData.RDS"))
}

# get the parameters for use in the report 
# get the values of transmission parameters of all the treatments
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

#detection results####
output.dir <- "./HPC_output/Detection/";
S1.results_all <- readRDS(list.files(output.dir, pattern = "*_Detection_results_complete_EU*\\.RDS$", full.names = TRUE));
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
                                            fraction_small = sum(final.size<= 0.001*46000)/n())


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


#"Transmission parameters from transmission experiments used in this study. For BI and BIBOOSTER the same values are used. All scenarios are run with the values for cut-off >=6. ", echo = FALSE, include = TRUE} ####
suppressMessages(knitr::kable(x= input_table %>%filter(cut_off %in% c("-","O6"))%>%select(!cut_off)%>%select(!cutoff_name),
 # mutate(cut_off = factor(cut_off, levels = c("-", "O6","O5", "O7"))) %>%
  #arrange(cut_off)%>%select(!cut_off), 
  #format = "html", 
  escape= FALSE,col.names = c("Vaccination","$\\beta_{low titre}$","Mean inf. per","2.5%","97.5%","$\\beta_{high titre}$","Mean inf. per","2.5%","97.5%")))

# Reproduction numbers from transmission experiments used in this study. For BI and BIBOOSTER the same values are used."}
# Create a data frame for the transmission coefficients from parameter list
#create input to the table - check with input!
input_table_R <- with(all_files_selection,{data.frame(
  Treatment = Treatment,
  cut_off = replace_na(cut_off,"-"),
  cutoff_name = replace_na(str_replace(cut_off,pattern = "O", replacement = " &gt;= "),"&#45;"),
  beta_low_titre = beta_low_titre * Tinf_low_titre,
  beta_high_titre = beta_high_titre *  Tinf_high_titre
  
)})%>%mutate(cut_off = factor(cut_off,levels = c("-", "O6","O5", "O7")))%>%arrange(cut_off)
row.names(input_table_R)<-NULL;
#create a knitr table 
suppressMessages(knitr::kable(x= input_table_R%>%select(!cut_off ), 
                              #format = "html",
                              digits = 2, escape= FALSE,col.names = c("Vaccination","cut-off","$R_{low titre}$","$R_{high titre}$")))


#### Proportion high titre ####

#load data of field experiment
data_field_study_cut_offs <- readRDS("./input/FieldExperiment/Processed_Data/data_field_study_cut_offs.RDS")

times <- seq( 0,588,1)

plot_field_data <- data_field_study_cut_offs%>%filter(!str_detect(Treatment, "Neg. controle"))%>%pivot_longer(cols = starts_with("O"), names_to = "cut_off")%>%mutate(cut_off = factor(cut_off, levels = c("O6","O5","O7")))%>%arrange(cut_off)

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

## Detection modules

# Results ####

## Reproduction number and critical high-titre proportion


pmajor_scenarios_crit <- pmajor_scenarios%>%mutate(production_phase = ifelse(time_point<120,"preproduction","production"),
  above_pc = as.numeric(p>pc))

average_pmajor_scenarios_crit<- pmajor_scenarios_crit%>%filter(cut_off == "O6")%>%reframe(.by = c(Treatment, production_phase),
                           prop_above = sum(above_pc)/n(),
                           prop_above_ll = clopper_pearson( sum(above_pc),n())[1],
                           prop_above_ul  = clopper_pearson( sum(above_pc),n())[2])%>%pivot_wider(names_from = production_phase, values_from = c(prop_above,prop_above_ll,prop_above_ul))


### CEVA

#"The reproduction number at sampling moments of the field study for a cut-off value of log 6 per farm. Dashed line indicates the line R = 1, which is the threshold for herd-immunity. "}
ggplot(pmajor_scenarios%>%filter(time_point>=120)%>%filter(!str_detect(Treatment,"BI"))%>%filter(cut_off %in% c("O6","NOVAC")))+
  geom_point(aes(x = time_point, y = Rv, colour = Treatment, shape = Farm))+
  scale_y_continuous(transform = "log10",breaks = c(0.5,1,1.5,2,3,5,10), limits = c(0.5,10))+
  scale_colour_discrete(labels = ceva_labels)+
  scale_shape_discrete(labels = list('A' = "Farm A", 'B' = "Farm B"))+
  labs(x = "Sample time (days)", y = expression(italic(R[v])))+
  #geom_vline(aes(xintercept = 120))+
  geom_hline(aes(yintercept =1), linetype = "dashed")+
  theme_minimal()+
  theme(legend.title = element_blank())
  

suppressMessages(ggsave("./figures/R_ceva.png"))



average_pmajor_scenarios_crit%>%filter(!str_detect(Treatment,"BI"))%>%select(1,2,4,6,3,5,7) %>%knitr::kable(digits =2, 
                                                                      #format = "html",
                                                                      caption = "Fraction and 95%-confidence interval of time points with a proportion high titre above the critical fraction during the preproducton phase (age below 120 days) and during production phase (age above 120 days).", col.names = c("Vaccination","preproduction", "2.5%", "97.5%","production", "2.5%", "97.5%"))


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


average_pmajor_scenarios_crit%>%filter(!str_detect(Treatment,"CEVA"))%>%select(1,2,4,6,3,5,7) %>%knitr::kable(digits =2, 
                                                                      #format = "html",
                                                                      caption = "Fraction and 95%-confidence interval of time points with a proportion high titre above the critical fraction during the preproducton phase (age below 120 days) and during production phase (age above 120 days).", col.names = c("Vaccination","preproduction", "2.5%", "97.5%","production", "2.5%", "97.5%"))


## Detection probability, detection time and consequences when detected
## Outbreak sizes after introduction of one randomly selected infectious bird

### CEVA
# "Final size distributions per production phase (before/after age 120 days) and vaccination without detection."}
x_bins = 25;
ggplot(S1.results%>%filter(!str_detect(Treatment, "BI")))+
  geom_histogram(aes((final.size-1)/46000, after_stat(density/x_bins)),bins = x_bins)+
  labs(x = "Fraction population infected during outbreak","Simulation count")+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_grid(production_phase~Treatment, scales = "free")


#To do get these value
# Vaccination leads to a larger proportion of small outbreaks (less than
# 0.1% of the flock). For unvaccinated flocks, this was
# `r small_outbreaks[small_outbreaks$Treatment == "NOVAC",]$fraction_small[1]`.
# In the preproduction phase this changes to values between
# `r min(small_outbreaks[small_outbreaks$Treatment != "NOVAC" & small_outbreaks$production_phase == "preproduction","fraction_small"]%>%filter(!str_detect(Treatment, "BI")))`
# and
# `r max(small_outbreaks[small_outbreaks$Treatment != "NOVAC" & small_outbreaks$production_phase == "preproduction","fraction_small"]%>%filter(!str_detect(Treatment, "BI")))`.In
# the production phase this changes to values between
# `r min(small_outbreaks[small_outbreaks$Treatment != "NOVAC" & small_outbreaks$production_phase == "production","fraction_small"]%>%filter(!str_detect(Treatment, "BI")))`
# and
# `r max(small_outbreaks[small_outbreaks$Treatment != "NOVAC" & small_outbreaks$production_phase == "production","fraction_small"]%>%filter(!str_detect(Treatment, "BI")))`.




### BI
#The final size is the fraction of the total population infected during the entire course of the outbreak if it is not controlled.

# "Final size distributions per production phase (before/after age 120 days) and vaccination without detection."}
x_bins = 25;
ggplot(S1.results%>%filter(!str_detect(Treatment, "CEVA")))+
  geom_histogram(aes((final.size-1)/46000, after_stat(density/x_bins)),bins = x_bins)+
  labs(x = "Fraction population infected during outbreak","Simulation count")+
  ylim(0,1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  facet_grid(production_phase~Treatment, scales = "free")

 
# Vaccination leads to a larger proportion of small outbreaks (less than
# 0.1% of the flock). For unvaccinated flocks, this was
# `r small_outbreaks[small_outbreaks$Treatment == "NOVAC",]$fraction_small[1]`.
# In the preproduction phase this changes to values between
# `r min(small_outbreaks[small_outbreaks$Treatment != "NOVAC" & small_outbreaks$production_phase == "preproduction","fraction_small"]%>%filter(!str_detect(Treatment, "CEVA")))`
# and
# `r max(small_outbreaks[small_outbreaks$Treatment != "NOVAC" & !str_detect(small_outbreaks$Treatment, "CEVA") & small_outbreaks$production_phase == "preproduction","fraction_small"])`.In
# the production phase this changes to values between
# `r min(small_outbreaks[small_outbreaks$Treatment != "NOVAC" & small_outbreaks$production_phase == "production","fraction_small"]%>%filter(!str_detect(Treatment, "CEVA")))`
# and
# `r max(small_outbreaks[small_outbreaks$Treatment != "NOVAC" & small_outbreaks$production_phase == "production","fraction_small"]%>%filter(!str_detect(Treatment, "CEVA")))`.
# 


## Detection results ####

### Fraction large and small outbreaks detected

#create variable to determine detection  by survey 

table_data <- S1.results%>%reframe(.by = c(production_phase,Treatment,outbreak.size.ofconcern,time.interval.ac),
  passive = sum(as.numeric(!is.infinite(det.since.intro.pas)))/n(),
  active  = sum(as.numeric(confirm.result == "POSITIVE" | live.confirmation.result == "POSITIVE" )*as.numeric(!is.infinite(det.since.intro.ac)))/n(),
  survey = sum(as.numeric(confirm.survey.result == "POSITIVE" | live.confirmation.survey.result == "POSITIVE")*as.numeric(!is.infinite(det.since.intro.survey)))/n(),
  missed = sum(missed_nosurvey)/n(),
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
                                                    p1_detection_time = quantile(min.det.time,0.25),
                                                    median_detection_time = median(min.det.time),
                                                    p2_detection_time = quantile(min.det.time,0.75)
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
                                         max_outbreak_size = quantile(detect.N.dead.detectables+ detect.N.live.detectables,0.75))%>%arrange(outbreak.size.ofconcern,Treatment, production_phase, time.interval.ac)





### CEVA

# size of detected and undetected introduction with a small number of infected birds (<= 0.1% of farm)., echo = FAlSE}

ggplot(S1.results%>%filter(outbreak.size.ofconcern == "SMALL")%>% filter(str_detect(Treatment,"CEVA")  ), aes(x = as.factor(time.interval.ac), y =detect.N.dead.detectables + detect.N.live.detectables, fill = !missed ))+geom_boxplot()+scale_fill_manual(values = c("darkblue","darkgreen"),name = "Introduction detected" )+ labs( x = "Active surveillance interval", y = "Total number of infected birds (at fade-out or detection moment)")+facet_grid(Treatment~., scales = "free_y")

suppressWarnings( ggsave("./figures/false_positives_CEVA.png") )

#"Proportion detected large outbreaks (>0.1% flock infected) by passive, active (dead bird i.e. bucket) or live bird survey for different frequencies of active sampleing. Last point is no active sampling.  ", echo = FALSE}


ggplot(correct_table_data%>%filter(production_phase == "production")%>%filter(!str_detect(Treatment, "BI"))%>%filter(is.finite(time.interval.ac))) +
  geom_path(aes(x = time.interval.ac, y = passive_LARGE, color = "Passive")) +
# geom_point(aes(x = time.interval.ac, y = passive_LARGE, color = "Passive")) +
#  geom_path(aes(x = time.interval.ac, y = active_LARGE, color = "Active")) +
  geom_point(aes(x = time.interval.ac, y = active_LARGE, color = "Active"), size = 2) +
  geom_path(aes(x = time.interval.ac, y = survey_LARGE, color = "Survey")) +
#  geom_point(aes(x = time.interval.ac, y = survey_LARGE, color = "Survey")) +
  scale_x_continuous(breaks = unique(detection_times_vs_interval$time.interval.ac)) +
  labs(x = "Sampling interval", y = "Fraction detected") +
  scale_color_manual(
    values = c("Passive" = "black", "Active" = "darkblue", "Survey" = "darkorange"),
    name = "Detection Method"
  ) +
  facet_grid(production_phase ~ Treatment, labeller = as_labeller(ceva_labels))

suppressWarnings( ggsave("./figures/CSe_CEVA.png") )

#create a table
knitr::kable(correct_table_data%>%filter(production_phase == "production")%>%filter(!str_detect(Treatment, "BI")), digits = 2)


### BI


#size of detected and undetected introduction with a small number of infected birds (<= 0.1% of farm)., echo=FALSE}

ggplot(S1.results%>%filter(outbreak.size.ofconcern == "SMALL" & str_detect(Treatment,"BI") ), aes(x = as.factor(time.interval.ac), y =detect.N.dead.detectables + detect.N.live.detectables, fill = !missed ))+geom_boxplot()+scale_fill_manual(values = c("darkblue","darkgreen"),name = "Introduction detected" )+facet_grid(Treatment ~.)+ labs( x = "Active surveillance interval", y = "Total number of infected birds (at fade-out or detection moment)")

suppressWarnings( ggsave("./figures/false_positives_BI.png") )

#"Proportion detected large outbreaks (>0.1% flock infected) by passive, active (dead bird i.e. bucket) or live bird survey for different frequencies of active sampleing. Last point is no active sampling.  ", echo = FALSE}

ggplot(correct_table_data%>%filter(production_phase=="production")%>%filter(!str_detect(Treatment, "CEVA"))%>%filter(is.finite(time.interval.ac))) +
  geom_path(aes(x = time.interval.ac, y = passive_LARGE, color = "Passive")) +
# geom_point(aes(x = time.interval.ac, y = passive_LARGE, color = "Passive")) +
#  geom_path(aes(x = time.interval.ac, y = active_LARGE, color = "Active")) +
  geom_point(aes(x = time.interval.ac, y = active_LARGE, color = "Active"), size = 2) +
  geom_path(aes(x = time.interval.ac, y = survey_LARGE, color = "Survey")) +
#  geom_point(aes(x = time.interval.ac, y = survey_LARGE, color = "Survey")) +
  scale_x_continuous(breaks = unique(detection_times_vs_interval$time.interval.ac)) +
  labs(x = "Sampling interval", y = "Fraction detected") +
  scale_color_manual(
    values = c("Passive" = "black", "Active" = "darkblue", "Survey" = "darkorange"),
    name = "Detection Method"
  ) +
  facet_grid(.~ Treatment, labeller = as_labeller(bi_labels))

suppressWarnings( ggsave("./figures/CSe_BI.png") )

#create a table
correct_table_data%>%filter(production_phase == "production")%>%filter(!str_detect(Treatment, "CEVA"))


### Detection times 
#The detection times and the variability are higher for vaccinated flocks. Depending on the proportion high titres in the flock (see difference preproduction vs production and different vaccines), the median detection times are similar to unvaccinated flock when every 2 to 14 days carcasses are tested (bucket sampling). Detailed outcomes of every time point are given in the appendix.

#### CEVA

# "Detection times in vaccinated and unvaccinated flocks. Solid line and dots are median values, grey area depicts minimum and maximum range. Detection times without active surveillance for unvaccinated flocks are indicated with dashed lines for median and dotted lines for min-max range."}


#create range of detection times without vaccination
det_interval_novac <- detection_times_vs_interval%>%filter(Treatment == "NOVAC" & is.infinite(time.interval.ac ))

#plot
ggplot(detection_times_vs_interval%>%filter(production_phase=="production")%>%filter(!str_detect(Treatment,"BI")), aes(x = time.interval.ac))+
  geom_ribbon(aes(ymin = p1_detection_time, ymax = p2_detection_time), alpha = 0.1)+
  geom_path(aes(y = median_detection_time))+
  geom_point(aes(y = median_detection_time))+
  geom_hline(aes(yintercept = det_interval_novac$p1_detection_time[1]),linetype = "dotted" )+
  geom_hline(aes(yintercept = det_interval_novac$p2_detection_time[1]),linetype = "dotted") +
  geom_hline(aes(yintercept = det_interval_novac$median_detection_time[1]),linetype = "dashed")+
  scale_x_continuous(breaks = unique(detection_times_vs_interval$time.interval.ac) )+
  facet_grid(production_phase~Treatment, labeller = as_labeller( ceva_labels))+labs(y = "detection time",x = "sampling interval")

suppressWarnings( ggsave("./figures/detection_times_ceva.png") )


detection_times_vs_interval%>%filter(production_phase == "production")

# "Detection time distributions of outbreaks > 0.1% of total farm size for passive surveillance only."}
    ggplot(S1.results%>%filter(outbreak.size.ofconcern == "LARGE" &
                                 is.finite(pas.det.time))%>%filter(!str_detect(Treatment,"BI")),aes(x =factor(Treatment), y = pas.det.time, fill = Treatment,colour = Treatment))+
   geom_jitter(alpha = 0.1)+
   geom_boxplot(colour = "darkgrey", fill = NA,alpha =0.1,outlier.shape = NA)+ 
  scale_colour_manual(values = okabe_ito[1:2])+
  stat_summary(
    fun = median,
    geom = "point",
    shape = 18,
    size = 3,
    colour = "black",
    position = position_dodge(width = 0.9)
  )+  facet_grid(.  ~ production_phase)+ 
  theme(legend.position = "none")+
     labs(x = "",y = "Detection time")#+ 
    # ggtitle("Detection time distribution with passive detection for an outbreak larger than 0.001 of flock.") 

ggsave("./figures/passive_detection_times_ceva.png")


# "Fade out times for runs where the introduction is not detected."}
#for each run get an fade-out time
fade_out_time <- scenarios_sim%>%filter(L.1+L.2+I.1+I.2 > 0)%>%reframe(.by = c(scenario, run, Treatment, production_phase),
                                         fade_out_at = max(time))
# ggplot(fade_out_time)+
#   geom_histogram(aes(x = fade_out_at, y = after_stat(density), fill = Treatment))+facet_grid(production_phase~.)

#determine undetected runs 
undetected_runs <- unique(S1.results%>%filter(is.infinite(min.det.time)& outbreak.size.ofconcern=="LARGE")%>%reframe(scenario, run, Treatment, time.interval.ac))

undet_fade_out_time<- fade_out_time%>%semi_join(y = undetected_runs,by = c("scenario","run","Treatment"))

undet_fade_out_time <- left_join(undet_fade_out_time,undetected_runs)%>%mutate(right_censor = as.numeric(fade_out_at>=89))

mean_undet_fade_out_times <- undet_fade_out_time%>%reframe(.by = c("time.interval.ac"),
                                                           mean_fade_out_time = mean(fade_out_at),
                                                           var_fade_out_time = var(fade_out_at),
                                                           )%>%mutate(alpha = mean_fade_out_time^2 /  var_fade_out_time,
                                                                      beta = mean_fade_out_time /  var_fade_out_time)

plot_data <- mean_undet_fade_out_times %>%
  mutate(
    fade_out_at = list(0:90),
    y = map2(alpha, beta, ~ dgamma(0:90, .x, .y))
  ) %>%
  unnest(cols = c(fade_out_at, y))%>%select(time.interval.ac, fade_out_at, y)

ggplot(undet_fade_out_time%>%filter(fade_out_at <= 89)%>%filter(!str_detect(Treatment,"BI")))+
  geom_histogram(aes(x = fade_out_at, y = after_stat(density), fill = Treatment))+facet_grid(.~is.finite(time.interval.ac))

                                                                                       

#### BI

#"Detection times in vaccinated and unvaccinated flocks. Solid line and dots are median values, grey area depicts minimum and maximum range. Detection times without active surveillance for unvaccinated flocks are indicated with dashed lines for median and dotted lines for min-max range."}


#create range of detection times without vaccination
det_interval_novac <- detection_times_vs_interval%>%filter(Treatment == "NOVAC" & is.infinite(time.interval.ac ))

#plot
ggplot(detection_times_vs_interval%>%filter(production_phase=="production")%>%filter(!str_detect(Treatment,"CEVA")), aes(x = time.interval.ac))+
  geom_ribbon(aes(ymin = p1_detection_time, ymax = p2_detection_time), alpha = 0.1)+
  geom_path(aes(y = median_detection_time))+
  geom_point(aes(y = median_detection_time))+
  geom_hline(aes(yintercept = det_interval_novac$p1_detection_time[1]),linetype = "dotted" )+
  geom_hline(aes(yintercept = det_interval_novac$p2_detection_time[1]),linetype = "dotted") +
  geom_hline(aes(yintercept = det_interval_novac$median_detection_time[1]),linetype = "dashed")+
  scale_x_continuous(breaks = unique(detection_times_vs_interval$time.interval.ac) )+
  facet_grid(.~Treatment, labeller = as_labeller(bi_labels))+labs(y = "detection time",x = "sampling interval")
suppressWarnings( ggsave("./figures/detection_times_bi.png") )

#"Detection time distributions of outbreaks > 0.1% of total farm size for passive surveillance only."}
    ggplot(S1.results%>%filter(outbreak.size.ofconcern == "LARGE" &
                                 is.finite(pas.det.time))%>%filter(!str_detect(Treatment,"CEVA")),aes(x =factor(Treatment), y = pas.det.time, fill = Treatment,colour = Treatment))+
  geom_jitter(alpha = 0.1)+
  geom_boxplot(colour = "darkgrey",fill = NA, alpha =0.1,outlier.shape = NA)+ 
  scale_colour_manual(values = okabe_ito[c(3,4,2)])+
  scale_fill_manual(values = okabe_ito[c(3,4,2)])+
  stat_summary(
    fun = median,
    geom = "point",
    shape = 18,
    size = 3,
    colour = "black",
    position = position_dodge(width = 0.9)
  )+  facet_grid(.  ~ production_phase)+ 
  theme(legend.position = "none")+
     labs(x = "",y = "Detection time")#+ 
    # ggtitle("Detection time distribution with passive detection for an outbreak larger than 0.001 of flock.") 

ggsave("./figures/passive_detection_times_bi.png")

# "Fade out times for runs where the introduction is not detected."}
#for each run get an fade-out time
fade_out_time <- scenarios_sim%>%filter(L.1+L.2+I.1+I.2 > 0)%>%reframe(.by = c(scenario, run, Treatment, production_phase),
                                         fade_out_at = max(time))
# ggplot(fade_out_time)+
#   geom_histogram(aes(x = fade_out_at, y = after_stat(density), fill = Treatment))+facet_grid(production_phase~.)

#determine undetected runs 
undetected_runs <- unique(S1.results%>%filter(is.infinite(min.det.time)& outbreak.size.ofconcern=="LARGE")%>%reframe(scenario, run, Treatment, time.interval.ac))

undet_fade_out_time<- fade_out_time%>%semi_join(y = undetected_runs,by = c("scenario","run","Treatment"))

undet_fade_out_time <- left_join(undet_fade_out_time,undetected_runs)%>%mutate(right_censor = as.numeric(fade_out_at>=89))

mean_undet_fade_out_times <- undet_fade_out_time%>%reframe(.by = c("time.interval.ac"),
                                                           mean_fade_out_time = mean(fade_out_at),
                                                           var_fade_out_time = var(fade_out_at),
                                                           )%>%mutate(alpha = mean_fade_out_time^2 /  var_fade_out_time,
                                                                      beta = mean_fade_out_time /  var_fade_out_time)

plot_data <- mean_undet_fade_out_times %>%
  mutate(
    fade_out_at = list(0:90),
    y = map2(alpha, beta, ~ dgamma(0:90, .x, .y))
  ) %>%
  unnest(cols = c(fade_out_at, y))%>%select(time.interval.ac, fade_out_at, y)

ggplot(undet_fade_out_time%>%filter(fade_out_at <= 89)%>%filter(!str_detect(Treatment,"CEVA")))+
  geom_histogram(aes(x = fade_out_at, y = after_stat(density), fill = Treatment))+facet_grid(.~is.finite(time.interval.ac))


#### distributions for between-farm


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

# ### Epidemic curves
# 
# In all three vaccination scenarios both high and low titre birds are infected during a large outbreak.
# 
# #
# scenarios_sim <- scenarios_sim%>%mutate(cumulative_infected = N.total - (S.1 + S.2),
#                                         cumulative_infected_1 = N.total*(1-p) - S.1,
#                                         cumulative_infected_2 = N.total*(p) - S.2)
# 
# #loop over the detection results 
# S1.results$cumulative_infected_active = NA;
# S1.results$cumulative_infected_active_1 = NA;
# S1.results$cumulative_infected_active_2 = NA;
# for(i in c(1:length(S1.results$run ))){
#   mdt <- S1.results$min.det.time[i];
#   S1.results$cumulative_infected_active[i] = scenarios_sim%>%filter(run == S1.results$run[i] &
#                        Treatment == S1.results$Treatment[i] &
#                        intro_time == S1.results$time_point[i] & time <= mdt)%>%
#     reframe(cumulative_infected_active = last(cumulative_infected))%>%as.numeric();
#   S1.results$cumulative_infected_active_1[i] = scenarios_sim%>%filter(run == S1.results$run[i] &
#                        Treatment == S1.results$Treatment[i] &
#                        intro_time == S1.results$time_point[i] & time <= mdt)%>%
#     reframe(cumulative_infected_active_1 = last(cumulative_infected_1))%>%as.numeric();
#   S1.results$cumulative_infected_active_2[i] = scenarios_sim%>%filter(run == S1.results$run[i] &
#                        Treatment == S1.results$Treatment[i] &
#                        intro_time == S1.results$time_point[i] & time <= mdt)%>%
#     reframe(cumulative_infected_active_2 = last(cumulative_infected_2))%>%as.numeric();
# }
# 
# ggplot(S1.results)+geom_histogram(aes(x = cumulative_infected_active))+facet_grid(production_phase ~ Treatment, scales = "free_x")
# 
# ggplot(S1.results)+
#   geom_histogram(aes(x = cumulative_infected_active_1, colour = okabo_ito[7]))+
#   geom_histogram(aes(x = cumulative_infected_active_2, colour = okabo_ito[6]))+
#   facet_grid(production_phase ~ Treatment, scales = "free_x")
# 
# ```
# 
# 
# <!-- ```{r, echo = FALSE, fig.cap = "Epidemic curves per scenario. I.1 = number of low titre birds, I.2 = number of high titre birds."} -->
# 
# <!-- #check runs by visualizing daily data -->
# 
# <!-- ggplot(data =scenarios_sim%>%mutate(p_interval = cut(p,breaks = c(-1,0 ,0.1, 0.5, 0.9, 1), -->
# 
# <!--   labels = c("0.0","0.0-0.1", "0.1-0.5", "0.5-0.9","0.9-1.0"), -->
# 
# <!--   include.lowest = TRUE))%>%select(day, I.1,I.2,Treatment, Farm,production_phase, run,p_interval)%>%pivot_longer(cols = c("I.1","I.2"),names_to = "state", values_to = "I"),  -->
# 
# <!-- aes(x = day, y = I, colour = p_interval,group = as.factor(paste(run,p_interval)))) + -->
# 
# <!--   geom_path( alpha = 0.1) + -->
# 
# <!--   labs(title = paste("Infectious birds"), x = "Day since introduction", y = "Infectious") + -->
# 
# <!--   theme_minimal() + -->
# 
# <!--   scale_colour_discrete(name = "Proportion high titre")+ -->
# 
# <!--   facet_grid(Treatment~ Farm+state+production_phase, scales = "free") +xlim(0,90) -->
# 
# <!-- #ggsave(filename = "./figures/epicurveCO6.png") -->
# 
# <!-- ``` -->
# 
# For detection, incidence curves are more important. In particular, those
# animals that can be detected when sampling occurs. Therefore, a graph is
# shown below that presents the cumulative number of dead detectable birds
# in time.
# 
# ```{r echo = FALSE, include = FALSE}
# #check runs by visualizing daily data
# n = 2
# m = 7
# 
# get_indices <- function(t, n, m) {
#   indices <- c()
#   x <- 0
#   while (TRUE) {
#     # Current day and previous n days
#     current <- t - c(0, 1:(n-1))
#       # Shifted back by x*m days
#     shifted <- current - x * m 
#     # Combine and add to indices
#     indices <- unique(c(indices, current, shifted))
#     # Check if all shifted indices are greater than 0
#     if (all(shifted > 0)) {
#       x <- x + 1
#     } else {
#       break
#     }
#   }
#   # Return only positive indices
#   return(indices[indices > 0])
# }
# 
# 
# plot_data_detection_incidence <- scenarios_sim%>%
#   mutate(.by = c(run,scenario,Farm),
#     cumulative_dead_detectables_2 = sapply(
#       day,
#       function(current_time) {
#         # Generate a vector with all values that in interval c(current_time - n, current_time) for each value of m*x smaller than current_time
#         sum(daily.death.incidence.detectables[ get_indices(current_time, n,2)], na.rm = TRUE)
#         
#       }
#     ),
#     cumulative_dead_detectables_7 = sapply(
#       day,
#       function(current_time) {
#         # Generate a vector with all values that in interval c(current_time - n, current_time) for each value of m*x smaller than current_time
#         sum(daily.death.incidence.detectables[ get_indices(current_time, n,7)], na.rm = TRUE)
#         
#       }
#     ),
#     cumulative_dead_detectables_14 = sapply(
#       day,
#       function(current_time) {
#         # Generate a vector with all values that in interval c(current_time - n, current_time) for each value of m*x smaller than current_time
#         sum(daily.death.incidence.detectables[ get_indices(current_time, n,14)], na.rm = TRUE)
#         
#       }
#     ),
#     cumulative_dead_detectables_30 = sapply(
#       day,
#       function(current_time) {
#         # Generate a vector with all values that in interval c(current_time - n, current_time) for each value of m*x smaller than current_time
#         sum(daily.death.incidence.detectables[ get_indices(current_time, n,30)], na.rm = TRUE)
#         
#       }
#     )
#   )
# 
# 
# 
# ```
# 
# <!-- ```{r, echo = FALSE, fig.cap="Cumulative number of dead detectables. "} -->
# 
# <!-- ggplot(plot_data_detection_incidence%>%mutate(introduction_time= as.factor(10*round(intro_time/10))))+ -->
# 
# <!--   geom_path(aes(x = time,y = N.dead.detectables, group = paste(run, intro_time)))+ -->
# 
# <!--   facet_grid(production_phase  ~ Treatment, scales  = "free" ) -->
# 
# <!-- ``` -->
# #### CEVA
# 
# ```{r, echo = FALSE, fig.cap="Cumulative number of dead detectables in the last two days. "}
# ggplot(plot_data_detection_incidence%>%select(time, starts_with("cumulative_dead_detectables_"), Treatment, production_phase, run,intro_time)%>%pivot_longer(cols = starts_with("cumulative_dead_detectables_"))%>%mutate(interval = str_remove(name, "cumulative_dead_detectables_"))%>%filter(!str_detect(Treatment, "BI")))+
#   geom_path(aes(x = time,y = value, group =interaction(run, intro_time, interval),
#                  colour = interval))+
#                 scale_colour_discrete(type = RColorBrewer::brewer.pal(6, 'Blues')[2:6])+
#   facet_grid(production_phase ~ Treatment, scales = "free" )#+theme(legen = "Introduction")
# ```
# 
# #### BI
# ```{r, echo = FALSE, fig.cap="Cumulative number of dead detectables in the last two days. "}
# ggplot(plot_data_detection_incidence%>%select(time, starts_with("cumulative_dead_detectables_"), Treatment, production_phase, run,intro_time)%>%pivot_longer(cols = starts_with("cumulative_dead_detectables_"))%>%mutate(interval = str_remove(name, "cumulative_dead_detectables_"))%>%filter(!str_detect(Treatment, "CEVA")))+
#   geom_path(aes(x = time,y = value, group =interaction(run, intro_time, interval),
#                  colour = interval))+
#                 scale_colour_discrete(type = RColorBrewer::brewer.pal(6, 'Blues')[2:6])+
#   facet_grid(production_phase ~ Treatment, scales = "free" )#+theme(legen = "Introduction")
# ```
# ## Consequences at detection
# 
# ### Infected eggs shipped off premises
# 
# Infected eggs could pose a risk for human and animal health. Only when
# assuming a distinct virus (i.e. less protection) a substantial number of
# eggs can be shipped before detection.
# 
# ```{r tab_tranport_I_eggs, echo = FALSE,include = TRUE}
# # #visualize results
# # ggplot(S1.results%>%filter(farm_type == "layer_farm"))+
# #   geom_histogram(aes(log10(detect.Egg_I_shipped+1),fill = scenario), alpha = 0.5)+
# #   ggtitle("Eggs shipped before detection")+
# #   facet_grid(treatment ~ variant, scales = "free_x")

#Table
tmp_table <- knitr::kable(S1.results%>%reframe(
  .by = c(Treatment,  time.interval.ac),
  fraction_runs_with_transport = sum(detect.Egg_I_shipped>=1,na.rm = TRUE)/n(),
  median_I_eggs = median(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),na.rm = TRUE), 
  minimum_I_eggs = min(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),na.rm = TRUE),
  p25_I_eggs = quantile(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),0.25,na.rm = TRUE),
  p75_I_eggs = quantile(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),0.75, na.rm = TRUE),
  maximum_I_eggs = max(round(detect.Egg_I_shipped[detect.Egg_I_shipped>=1]),na.rm = TRUE)
  )
, digits = 2, 
#format = "html",
col.names = c("Vaccine","freq.","Runs with transport of infected eggs", "median","min.","25%","75%" ,"max."),
caption = "Shipped infected eggs. For the runs with at least 1 transported infected eggs the median, minimum and maximum number of shipped infected eggs.")

#as_image(tmp_table, file = "./figures/tab_shippedIeggs.png")

tmp_table


### Number of dead birds that were or are detectable at the moment of detection



# ggplot(S1.results)+
#   geom_histogram(aes(detect.N.dead.detectables,fill = scenario), alpha = 0.5)+
#   ggtitle("Number of dead detectable birds at detection")+
#   facet_grid(farm_type ~ treatment, scales = "free_x")

#Table
tmp_table <- knitr::kable(S1.results%>%reframe(
  .by = c(Treatment, time.interval.ac),
  median_I_eggs = median(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),na.rm = TRUE), 
  minimum_I_eggs = min(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),na.rm = TRUE),
  p25_I_eggs = quantile(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),0.25,na.rm = TRUE),
  p75_I_eggs = quantile(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),0.75,na.rm = TRUE),
  maximum_I_eggs = max(round(detect.N.dead.detectables[detect.N.dead.detectables>=1]),na.rm = TRUE)
  )
, digits = 0,  
#format = "html",
col.names = c("Vaccine","freq.", "median","min.","25%","75%" ,"max."),
caption = "Median, minimum and maximum number of dead detectable birds at moment of detection.")

#as_image(tmp_table, file = "./figures/tab_deadIbirds.png")

tmp_table


### Number of live birds that are detectable at the moment of detection
 
  #Table
tmp_table <- knitr::kable(S1.results%>%reframe(
  .by = c(Treatment, time.interval.ac),
  median_I_eggs = median(round(detect.N.live.detectables[detect.N.live.detectables>=1]),na.rm = TRUE), 
  minimum_I_eggs = min(round(detect.N.live.detectables[detect.N.live.detectables>=1]),na.rm = TRUE),
  p25_I_eggs = quantile(round(detect.N.live.detectables[detect.N.live.detectables>=1]),0.25,na.rm = TRUE),
  p75_I_eggs = quantile(round(detect.N.live.detectables[detect.N.live.detectables>=1]),0.75,na.rm = TRUE),
  maximum_I_eggs = max(round(detect.N.live.detectables[detect.N.live.detectables>=1]),na.rm = TRUE)
  ),
  #format = "html",
, digits = 0, 
col.names = c("Vaccine","Freq.", "median","min.","25%","75%" ,"max."),
caption = "Median, minimum and maximum number of live detectable birds at moment of detection.")

#as_image(tmp_table, file = "./figures/tab_liveIbirds.png")

tmp_table



# 
# 
# <!-- ## Sensitivity analyses -->
# # Appendix A Detailed results
# 
# ## Fraction of outbreaks detected
# 
# ```{r tab_det_small, echo = FALSE, tab.cap="Probability to detect a small outbreak (<= 0.1% of flock) by passive, active with confirmation. Freq. = sampling interval (Inf = no bucket sampling only live birds), Shortest = fraction shortest detection time is by active surveillance. "}
# 
#  knitr::kable(table_data[,c("production_phase","Treatment","time.interval.ac", "passive_SMALL","active_SMALL","survey_SMALL","shortest_active_SMALL","missed_SMALL")], digits = 2,
#              #format = "html",
#              align = 'c',col.names = c("Production phase","Vaccine","Freq.","Passive"," Active","Survey","Shortest","Missed"))
# 
# 
# 
# ```
# 
# ```{r tab_det_large, echo = FALSE, tab.cap="Probability to detect a large outbreak (> 0.1% of the flock) by passive, active with confirmation. Freq. = sampling interval (Inf = no bucket sampling only live birds), Shortest = fraction shortest detection time is by active surveillance. "}
# 
# 
#  knitr::kable( table_data[,c("production_phase","Treatment","time.interval.ac", "passive_LARGE","active_LARGE","survey_LARGE","shortest_active_LARGE","missed_LARGE")], digits = 2,
#              #format = "html",
#              align = 'c',col.names = c("Production phase","Vaccine","Freq.","Passive"," Active","Survey","Shortest","Missed"))
# 
# 
# 
# ```
# 
# ## Detection times per time point 
# 
# ```{r fig_det_dist_large_sample_moment,echo = FALSE, fig.cap = "Detection time distributions of large outbreaks for different active surveillance intervals at each time point. Black diamond indicates median value."}
# 
# # Custom labeller function
# custom_labeller <- function(x) {
#   ifelse(is.infinite(x), "No active surveillance", as.character(x))
# }
# 
# 
#     ggplot(S1.results%>%filter(outbreak.size.ofconcern == "LARGE" & is.finite(min.det.time)),
#            aes(x =as.factor(round(time_point/10)*10), y = min.det.time, fill = Treatment,colour = Treatment))+
#     geom_violin()+
#   stat_summary(
#     fun = median,
#     geom = "point",
#     shape = 18,
#     size = 1,
#     colour = "black",
#     position = position_dodge(width = 0.9)
#   )+ 
#     facet_grid(time.interval.ac  ~., scales = "free", labeller = as_labeller(custom_labeller))+ 
#     labs(x = "Moment of introduction",y = "Detection time")#+ 
#     #ggtitle("Detection time distribution with passive and active detection for a large outbreak.") 
# 
# #ggsave("./figures/detection_times.png")
# ```
# 
# ## outbreak size undetected 
# 
# ```{r fig_outbreak_size_undetected,echo = FALSE, fig.cap= "Outbreak size of undetected outbreaks by passive surveilllance, active bucket surveillance (different frequencies with Inf = no bucket surveilance) and live bird survey (30 day frequency) is in place."}
# ggplot(data = S1.results%>%filter(is.infinite(min.det.time)),aes(x = final.size, y = after_stat(density)))+
#   geom_histogram(aes(fill = as.factor(time.interval.ac)))+
#   scale_fill_discrete(name = "Freq. bucket")+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
#   facet_grid( outbreak.size.ofconcern~Treatment , scales = "free_y")
# 
# ```
# 
# 
# # Appendix B Sanity checks
# 
# In this section, a number of checks are reported to show whether
# simulations perform as expected.
# 
# ## Appendix B.1 probability of a major outbreak
# 
# The probability of a major outbreak was calculated according to
# [@Clancy2013TheInfection] [@Nandi2019StochasticSusceptibility] for
# gamma-distributed infectious periods and two types at infection given an
# introduction through an low-titre bird.
# 
# Without vaccination, the basic reproduction number is $7.9$ the
# probability of a major outbreak is over 90% (i.e.
# $$1 - \dfrac{1}{R_0} \approx 1 - \dfrac{1}{7.9} = 0.87$$ ).
# 
# ```{r echo=FALSE, fig.cap= "Probability of a major outbreak given vaccination. "}
# # Logit transformation function
# logit <- function(p) log(p / (1 - p))
# 
# ggplot(pmajor_scenarios)+
#   geom_point(aes(x = time_point, y = pmajor, colour = Treatment, shape = Farm))+
#   geom_vline(aes(xintercept = 120))+
#   ylim(0,1)+
#   #  scale_y_continuous(
#   #    trans = "logit",  # Use identity since we pre-transformed
#   #    breaks = c(0.0001,0.001,0.1,0.90,0.9999),
#   #    limits = c(0.0001, 0.9999)
#   #      )+
#    facet_grid(cut_off ~.)#, scale = "free")
# ```
# 
# Comparison between the fraction simulated and predicted (calculated)
# probability of a major outbreak $p_{major}$ shows a good
# correspondence.Difference occurs when the immunity is building up.
# Although a major outbreak is more likely to occur at that specific
# moment in time, the relatively fast build-up causes some lower numbers
# in the simulations (e.g. scenario CEVA_preproduction_half or
# BI_layer_farm).
# 
# ```{r, echo = FALSE, fig.cap = "Comparison between calculated and simulated probability of a major outbreak" }
# #calculate proportion of major outbreaks
# joined_S1_results <- left_join(S1.results %>%  reframe(.by =c(Treatment, cut_off, time_point),
#           p_major_sim = sum(final.size >=0.05*46000, na.rm = TRUE)/n()) , 
#                                pmajor_scenarios, by = c("cut_off","Treatment","time_point")) 
# 
# #table
# joined_S1_results%>%select(Treatment, cut_off,time_point,p_major_sim, pmajor,Rv)%>% knitr::kable(digits = 2, 
#                                                                                                  #format = "html", 
#                                                                                                  col.names = c("Treatment","Cut off","Time point","simulated $p_{major}$","predicted $p_{major}$", "Reproduction number") ,caption = "Proportion of major outbreaks per scenario.  ")#%>%as_image(file = "./figures/sim_vs_pred_pmajor.png")
# 
# #graph
# n_binom <- max(scenarios_sim$run)
# plot_data<- joined_S1_results%>%select(Treatment, cut_off,time_point,p_major_sim, pmajor)%>%mutate(.by = c(time_point,Treatment),
#   p_major_sim_ll =clopper_pearson(n_binom*p_major_sim,n_binom)[1],
#   p_major_sim_ul = clopper_pearson(n_binom*p_major_sim,n_binom)[2]
#   
# )#%>%pivot_longer(cols = c("p_major_sim","pmajor"),
#                                                                                           #          names_to = "sim_calc",
#                                                                                            #         values_to = "pm")
# 
# ggplot(plot_data)+
#    geom_point(aes(x = time_point, y = p_major_sim, colour = Treatment))+
#    geom_errorbar(aes(x = time_point, ymin = p_major_sim_ll,ymax = p_major_sim_ul, colour = Treatment))+
#   geom_point(aes(x = time_point, y = pmajor, colour = Treatment), shape =2)+
#   facet_grid(
#     .~Treatment)
#   
# # ggplot(plot_data)+
# #   geom_point(aes(x = time_point, y = pm, shape = sim_calc, colour = Treatment))+
# #   geom_line(aes(x = time_point, y = pm, group = paste(Treatment, time_point, cut_off), colour = Treatment))
# 
# #graph difference b4etween simulation and calculation
# # plot_data <- joined_S1_results%>%select(Treatment, cut_off,time_point,p_major_sim, pmajor)%>%mutate(dif_sim_calc = p_major_sim-pmajor)
# # ggplot(plot_data)+
# #   geom_point(aes(x = time_point, y = dif_sim_calc, colour = Treatment))+
# #   geom_segment(aes(x = time_point, xend = time_point,  y = dif_sim_calc,yend = 0,colour = Treatment))+
# #   facet_grid(cut_off~Treatment)+geom_vline(xintercept = 120)+ylim(c(-1,1)*1.05*max(abs(plot_data$dif_sim_calc)))
# 
# ```
# 
# # References
