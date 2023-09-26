#########################################################
#                                                        
#                 Spatial Scenarios                            
#                                          
#                                                        
#                  Author:  Egil Fischer
#                  Contact:                              
#                  Creation date                         
#########################################################
source("./src/loadLibraries.R")
source("./src/probMajorOutbreak.R")
source("./src/FittedWithinFarm.R")

#load location data ####
spatial.input<- read_excel("./input/20230404_AI_AnimalLocations_SizeType_v02.xlsx", sheet = 2)

#scale coordinates to km
spatial.input$Xkm <- spatial.input$X/1000
spatial.input$Ykm <- spatial.input$Y/1000
spatial.input$host_id <- c(1:nrow(spatial.input))

#screening and culling
culling.radius <- 1
culling.delay <- 1

# Define the transmission kernel and calculate the hazard matrix
# Parameters h0 used in Boender et al. 
h0 <- 0.002;
alpha <- 2.1;
r0 <- 1.9;
h_kernel <- function(r){h0/(1 + (r/r0)^alpha)} ; # transmission kernel as a function of r
beta<-1;

#set parameters
param.list.baseline.layer <- list(
  scenario = "baseline_Layer", #scenario
  runs = 10, #number of runs
  max.time = 17*30,#length of the run
  itypes = 2, #type
  N0 = 32000, #population size - agramatie website
  initial= 10 , #initially infected - choosen value
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = 2),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976 and Gemeraard et al 2023 #Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2/20, #Variance infectious period - Hobbelen et al uses shape parameter of 20 -> variance = m^2 / shape 
  #transRate = matrix(c(0,0.0,0.0,0), nrow = 2), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(1.00,0.001),#probability of dying at end of infectious period (for vaccinated -> 0 out of 20 -> 97.5%interval = 0.001)
  mortRate = 0.0005/7 #per capita death rate #Mortality events - based on performance reports and pers. com. mieke matthijs 0.5% per week +Gonzales & elbers figure 1 
)

#generic function for the introduction time which will affect the vaccination
intro_func<- function(type){
  if(type == "LAYER") return(runif(1,min = 0, max =17*30))
  if(type == "BROILER") return(runif(1,min = 0, max =42))
  else return(runif(1, min =0, max = 17*30))#assume duck and turkey same as LAYER
}


#introduction locations
set.seed(23011977)
spatial.input$density1KM <- pointdensity(spatial.input[,c("Xkm","Ykm")], eps = 1) 

hpda <- sample(c(spatial.input%>%
                   filter(spatial.input$density1KM>= quantile(spatial.input$density1KM,0.95))%>%select("host_id"))$host_id, 100, replace = TRUE)
#spatial.input%>%filter(nr %in% hpda)

lpda <- sample(c(spatial.input%>%
                   filter(spatial.input$density1KM<= quantile(spatial.input$density1KM,0.25))%>%select("host_id"))$host_id, 100, replace = TRUE)
#spatial.input%>%filter(nr %in% lpda)


#Simulations ####
######################################################
#  High Poultry Density Area                         #
######################################################
index.farm <- hpda
max.runs <- length(index.farm)
####################
#no vaccination ####
####################
#set functions ####
vac.func <- function(type, intro){0}
Tdist = Tdist.novac.pas
fadeout_func = fadeout_func.novac.pas
exposure.function <- exposure.function.novac.pas

#run the simulations ####
srm <- Sys.time()
set.seed(23011977)
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_novaccination_hdpa.RDS")
saveRDS(histories, file = "./output/histories_novaccination_hdpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_hdpa.RDS")
#join the histories
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joinedhistories_novaccination_hdpa.RDS")
saveRDS(joined.vacstatus, file = "./output/joinedvacstatus_novaccination_hdpa.RDS")
Sys.time()-srm


#read saved simulations####
infected_over_time_runs<- readRDS( file = "./output/infected_over_time_runs_novaccination_hdpa.RDS")
histories <- readRDS(file = "./output/histories_novaccination_hdpa.RDS")
V_stat_matrix <- readRDS(file = "./output/V_stat_matrix_hdpa.RDS")
joined.histories<- readRDS(file = "./output/joinedhistories_novaccination_hdpa.RDS")
joined.vacstatus<- readRDS( file = "./output/joinedvacstatus_novaccination_hdpa.RDS")


#length of infectious periods
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)
#add vaccination status, farm type and size 
inf.vac.novaccin<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))

#determine the exposure ###
inf.vac.novaccin$exposure <-exposure.function(inf.vac.novaccin)

inf.vac.novaccin$exposure[inf.vac.novaccin$exposure<0]<-500
saveRDS(inf.vac.novaccin, "./output/inf.vac.novaccin.RDS")
inf.vac.novaccin<-readRDS( "./output/inf.vac.novaccin.RDS")
exposure.novaccin.tot <- inf.vac.novaccin%>%reframe(.by = run,
                           tot.exposure = sum(exposure))

mean.exposure.novaccin <- mean(exposure.novaccin.tot$tot.exposure)


#############################################################################################
#50% high titre -> exposure per farm is slightly less, but detection times are increased ####
#############################################################################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",0.5,0)}
Tdist = Tdist.vac.pas
fadeout_func = fadeout_func.vac.pas
exposure.function <- exposure.function.vac.pas

#run simulations ####
srm <- Sys.time()
set.seed(23011977)
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_60vaccination_hdpa.RDS")
saveRDS(histories, file = "./output/histories_60vaccination_hdpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_60vaccination_hdpa.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,
                                                   spatial.input$TYPE,
                                                   spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joinedhistories_60vaccination_hdpa.RDS")
saveRDS(joined.vacstatus, file = "./output/joinedvacstatus_60vaccination_hdpa.RDS")

Sys.time()-srm

#read saved simulations####
infected_over_time_runs<- readRDS( file = "./output/infected_over_time_runs_60vaccination_hdpa.RDS")
histories <- readRDS(file = "./output/histories_60vaccination_hdpa.RDS")
V_stat_matrix <- readRDS(file = "./output/V_stat_matrix_60vaccination_hdpa.RDS")
joined.histories<- readRDS(file = "./output/joinedhistories_60vaccination_hdpa.RDS")
joined.vacstatus<- readRDS( file = "./output/joinedvacstatus_60vaccination_hdpa.RDS")

#calculate infectious period ####
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.60vac<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))


#determine exposure ####

inf.vac.60vac$exposure <- exposure.function(inf.vac.60vac)
saveRDS(inf.vac.60vac, "./output/inf.vac.60vac.RDS")
inf.vac.60vac<-readRDS( "./output/inf.vac.60vac.RDS")

exposure.60vaccin.tot <- inf.vac.60vac%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))

##################################################################################################################
#Clinical protection Vaccination has 50% high titre with reduced transmission rate and 50% death of low titre birds ####
#################################################################################################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",0.5,0)}
Tdist = Tdist.clinprot.pas
fadeout_func = fadeout_func.clinprot.pas
exposure.function <- exposure.function.clinprot.pas

#sims####
srm <- Sys.time()
# vac.func <- function(x){ifelse(x == "LAYER",10^-6 #does not effect pmajor but changes the infectious period
#                                ,0)} 
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_clinprot_hdpa.RDS")
saveRDS(histories, file = "./output/histories_clinprot_hdpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_clinprot_hdpa.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joined_histories_clinprot_hdpa.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_clinprot_hdpa.RDS")
Sys.time()-srm
             
#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_clinprot_hdpa.RDS")
histories<- readRDS(file = "./output/histories_clinprot_hdpa.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_clinprot_hdpa.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_clinprot_hdpa.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_clinprot_hdpa.RDS")


#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.clinprot<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))
inf.vac.clinprot$exposure<- exposure.function(inf.vac.clinprot)
saveRDS(inf.vac.clinprot, "./output/inf.vac.clinprot.RDS")
inf.vac.clinprot<-readRDS( "./output/inf.vac.clinprot.RDS")

exposure.clinprot.tot <- inf.vac.clinprot%>%reframe(.by = run,
                                                   tot.exposure = sum(exposure))


######################################
#Waning of immunity - homologuous ####
######################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",
                                         1-pgamma(intro,36,0.07),0)}
Tdist = Tdist.vac.wane.pas
fadeout_func = fadeout_func.wane.pas
exposure.function <- exposure.function.wane.pas

#sims####
srm <- Sys.time()
# vac.func <- function(type){ifelse(type == "LAYER",
#                                   1-pgamma(runif(1,min = 0,max = 540), shape = 36,rate =0.07),
#                                   0)}
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_wane_hdpa.RDS")
saveRDS(histories, file = "./output/histories_wane_hdpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_wane_hdpa.RDS")
saveRDS(Intro_matrix,file = "./output/Intro_matrix_wane_hdpa.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joined_histories_wane_hdpa.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_wane_hdpa.RDS")
Sys.time()-srm

#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_wane_hdpa.RDS")
histories<- readRDS(file = "./output/histories_wane_hdpa.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_wane_hdpa.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_wane_hdpa.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_wane_hdpa.RDS")
Intro_matrix<-readRDS(file = "./output/Intro_matrix_wane_hdpa.RDS")

#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.wane<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))
#add introduction times
intro.times.run.id <- NULL;
for(j in c(1:length(Intro_matrix[,1]))){
  intro.times.run.id<- rbind(intro.times.run.id,
                             data.frame(intro.time = Intro_matrix[j,],
                                        host_id = spatial.input$host_id,
                                        run = j))
  
}
inf.vac.wane<- left_join(inf.vac.wane,intro.times.run.id,by = c("host_id","run"))

#no fit 

inf.vac.wane$exposure<- exposure.function(inf.vac.wane) 
saveRDS(inf.vac.wane, "./output/inf.vac.wane.RDS")
inf.vac.wane<-readRDS( "./output/inf.vac.wane.RDS")

exposure.wane.tot <- inf.vac.wane%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))

######################################
#Waning of immunity - heterologuous ####
######################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",
                                         0.7*(1-pgamma(intro, shape = 4,rate =0.1/7)),0)}
Tdist = Tdist.vac.wane.ds.pas
fadeout_func = fadeout_func.wane.ds.pas
exposure.function <- exposure.function.wane.ds.pas

#sims####
srm <- Sys.time()

source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_wane_ds_hdpa.RDS")
saveRDS(histories, file = "./output/histories_wane_ds_hdpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_wane_ds_hdpa.RDS")
saveRDS(Intro_matrix, file = "./output/Intro_matrix_wane_ds_hdpa.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;joined.intro =NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,
                                                   spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
  joined.intro <- rbind(joined.intro,cbind(spatial.input$host_id,
                                                   spatial.input$TYPE,spatial.input$SIZE, Intro_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
names(joined.intro)<- c("host_id","TYPE","SIZE","intro.time","run")
saveRDS(joined.histories, file = "./output/joined_histories_wane_ds_hdpa.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_wane_ds_hdpa.RDS")
saveRDS(joined.intro,file = "./output/intro_wane_ds_hdpa.RDS")
Sys.time()-srm

#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_wane_ds_hdpa.RDS")
histories<- readRDS(file = "./output/histories_wane_ds_hdpa.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_wane_ds_hdpa.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_wane_ds_hdpa.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_wane_ds_hdpa.RDS")
joined.intro<- readRDS(file = "./output/intro_wane_ds_hdpa.RDS")
Intro_matrix<-readRDS(file = "./output/Intro_matrix_wane_ds_hdpa.RDS")

#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.wane.ds<- left_join(infectious.periods,joined.vacstatus[,  c("host_id","run","vac")], by = c("host_id","run"))
inf.vac.wane.ds<- left_join(inf.vac.wane.ds, joined.intro[,  c("host_id","run","intro.time")], by = c("host_id","run"))
inf.vac.wane.ds<- left_join(inf.vac.wane.ds,spatial.input[,c("host_id","TYPE","SIZE")], by = c("host_id"))
                  
#no fit 
inf.vac.wane.ds$exposure<- exposure.function(inf.vac.wane.ds) 
#replace negative numbers
inf.vac.wane.ds$exposure[inf.vac.wane.ds$exposure<100]<- 100
saveRDS(inf.vac.wane.ds, "./output/inf.vac.wane.ds.RDS")
inf.vac.wane.ds<-readRDS( "./output/inf.vac.wane.ds.RDS")


exposure.wane.ds.tot <- inf.vac.wane.ds%>%reframe(.by = run,
                                            tot.exposure = sum(exposure))


############################################################
# Combine SCENARIOS                                        #
############################################################

#summary plots for HDPA ###
infected_over_time_runs.scenarios<- rbind(cbind(readRDS(file = "./output/infected_over_time_runs_novaccination_hdpa.RDS"), scenario = "1. Baseline"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_clinprot_hdpa.RDS"), scenario = "2. Reduced transmission"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_60vaccination_hdpa.RDS"), scenario = "3. Clinical protection"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_wane_hdpa.RDS"), scenario = "4. Waning immunity \n homologous"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_wane_ds_hdpa.RDS"), scenario = "5. Waning immunity \n heterologous"))

plot.data.scenarios <- infected_over_time_runs.scenarios%>%reframe(.by = c(run,scenario),
                                            Detected = max(C)+max(PC))%>%reshape2::melt(id.vars=c("run","scenario"))
ggplot(plot.data.scenarios, aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
                   scale_fill_manual(values =c("Detected" = "red", "Screened" = "gray"), 
                                     labels =c("Detected" = "Infected", "Screened" = "Screened"),
                                     name = "Culling")+
                   labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(.~scenario)+
#  ggtitle("Infected farms")+
  theme(legend.position = "none")
ggsave("./output/figures/InfectedFarmsHPDA.png")

############################################################
# Combine human exposure with the no vaccination scenario   #
############################################################
hpda.exposure.total <- rbind(cbind(exposure.novaccin.tot, vac = 0, scenario = "1. Baseline"),
                             cbind(exposure.60vaccin.tot, vac = 0.6, scenario = "2. Reduced transmission"),
                             cbind(exposure.clinprot.tot, vac = 0.6, scenario = "3. Clinical protection"),
                             cbind(exposure.wane.tot, vac = 0.6, scenario = "4. Waning immunity \n homologous"),
                             cbind(exposure.wane.ds.tot, vac = 0.6, scenario = "5. Waning immunity \n heterologous"))
hpda.exposure.total$ratio <- hpda.exposure.total$tot.exposure/mean.exposure.novaccin

# Visualize Human exposure ####

ggplot(hpda.exposure.total, aes(ratio))+ 
  geom_histogram(aes(y =after_stat(count/sum(count))),position = "identity",binwidth = 0.5)+
  labs(x = "log10 ratio of exposure",y = "Proportion of runs")+
  geom_vline(xintercept = 1)+
  scale_x_log10()+
  facet_grid(.~scenario)
#  ggtitle("Exposure during an outbreak of multiple farms")
ggsave("./output/figures/ExposureSpatial.hpda.png", scale = 1.23)

hpda.exposure.total%>%reframe(.by = scenario,
                              min.ratio = min(ratio),
                              max.ratio = max(ratio))


###################################################
#         LPDA                                    #
###################################################
index.farm <- lpda
max.runs <- length(index.farm)
####################
#no vaccination ####
####################
#set functions ####
vac.func <- function(type, intro){0}
Tdist = Tdist.novac.pas
fadeout_func = fadeout_func.novac.pas
exposure.function <- exposure.function.novac.pas

#run the simulations ####
srm <- Sys.time()
set.seed(23011977)
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_novaccination_ldpa.RDS")
saveRDS(histories, file = "./output/histories_novaccination_ldpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_ldpa.RDS")
#join the histories
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joinedhistories_novaccination_ldpa.RDS")
saveRDS(joined.vacstatus, file = "./output/joinedvacstatus_novaccination_ldpa.RDS")
Sys.time()-srm

#read saved simulations####
infected_over_time_runs<- readRDS( file = "./output/infected_over_time_runs_novaccination_ldpa.RDS")
histories <- readRDS(file = "./output/histories_novaccination_ldpa.RDS")
V_stat_matrix <- readRDS(file = "./output/V_stat_matrix_ldpa.RDS")
joined.histories<- readRDS(file = "./output/joinedhistories_novaccination_ldpa.RDS")
joined.vacstatus<- readRDS( file = "./output/joinedvacstatus_novaccination_ldpa.RDS")




#length of infectious periods
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)
#add vaccination status, farm type and size 
inf.vac.novaccin<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))

#determine the exposure ###
inf.vac.novaccin$exposure <-exposure.function(inf.vac.novaccin)

inf.vac.novaccin$exposure[inf.vac.novaccin$exposure<0]<-500
exposure.novaccin.tot <- inf.vac.novaccin%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))

mean.exposure.novaccin <- mean(exposure.novaccin.tot$tot.exposure)
#############################################################################################
#50% high titre -> exposure per farm is slightly less, but detection times are increased ####
#############################################################################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",0.5,0)}
Tdist = Tdist.vac.pas
fadeout_func = fadeout_func.vac.pas
exposure.function <- exposure.function.vac.pas

#run simulations ####
srm <- Sys.time()
set.seed(23011977)
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_60vaccination_ldpa.RDS")
saveRDS(histories, file = "./output/histories_60vaccination_ldpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_60vaccination_ldpa.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,
                                                   spatial.input$TYPE,
                                                   spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joinedhistories_60vaccination_ldpa.RDS")
saveRDS(joined.vacstatus, file = "./output/joinedvacstatus_60vaccination_ldpa.RDS")

Sys.time()-srm

#read saved simulations####
infected_over_time_runs<- readRDS( file = "./output/infected_over_time_runs_60vaccination_ldpa.RDS")
histories <- readRDS(file = "./output/histories_60vaccination_ldpa.RDS")
V_stat_matrix <- readRDS(file = "./output/V_stat_matrix_60vaccination_ldpa.RDS")
joined.histories<- readRDS(file = "./output/joinedhistories_60vaccination_ldpa.RDS")
joined.vacstatus<- readRDS( file = "./output/joinedvacstatus_60vaccination_ldpa.RDS")

#calculate infectious period ####
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.60vac<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))

#determine exposure ####

inf.vac.60vac$exposure <- exposure.function(inf.vac.60vac)
exposure.60vaccin.tot <- inf.vac.60vac%>%reframe(.by = run,
                                                 tot.exposure = sum(exposure))

##################################################################################################################
#Clinical protection Vaccination has 50% high titre with reduced transmission rate and 50% death of low titre birds ####
#################################################################################################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",0.5,0)}
Tdist = Tdist.clinprot.pas
fadeout_func = fadeout_func.clinprot.pas
exposure.function <- exposure.function.clinprot.pas

#sims####
srm <- Sys.time()
# vac.func <- function(x){ifelse(x == "LAYER",10^-6 #does not effect pmajor but changes the infectious period
#                                ,0)} 
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_clinprot_ldpa.RDS")
saveRDS(histories, file = "./output/histories_clinprot_ldpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_clinprot_ldpa.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joined_histories_clinprot_ldpa.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_clinprot_ldpa.RDS")
Sys.time()-srm

#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_clinprot_ldpa.RDS")
histories<- readRDS(file = "./output/histories_clinprot_ldpa.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_clinprot_ldpa.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_clinprot_ldpa.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_clinprot_ldpa.RDS")


#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.clinprot<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))


inf.vac.clinprot$exposure<- exposure.function(inf.vac.clinprot)
exposure.clinprot.tot <- inf.vac.clinprot%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))


######################################
#Waning of immunity - homologuous ####
######################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",
                                         1-pgamma(intro,36,0.07),0)}
Tdist = Tdist.vac.wane.pas
fadeout_func = fadeout_func.wane.pas
exposure.function <- exposure.function.wane.pas

#sims####
srm <- Sys.time()
# vac.func <- function(type){ifelse(type == "LAYER",
#                                   1-pgamma(runif(1,min = 0,max = 540), shape = 36,rate =0.07),
#                                   0)}
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_wane_ldpa.RDS")
saveRDS(histories, file = "./output/histories_wane_ldpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_wane_ldpa.RDS")
saveRDS(Intro_matrix,file = "./output/V_stat_matrix_wane_ldpa.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joined_histories_wane_ldpa.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_wane_ldpa.RDS")
Sys.time()-srm

#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_wane_ldpa.RDS")
histories<- readRDS(file = "./output/histories_wane_ldpa.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_wane_ldpa.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_wane_ldpa.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_wane_ldpa.RDS")


#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.wane<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))
#add introduction times
intro.times.run.id <- NULL;
for(j in c(1:length(Intro_matrix[,1]))){
  intro.times.run.id<- rbind(intro.times.run.id,
                             data.frame(intro.time = Intro_matrix[j,],
                                        host_id = spatial.input$host_id,
                                        run = j))
  
}
inf.vac.wane<- left_join(inf.vac.wane,intro.times.run.id,by = c("host_id","run"))

#no fit 

inf.vac.wane$exposure<- exposure.function(inf.vac.wane) 
exposure.wane.tot <- inf.vac.wane%>%reframe(.by = run,
                                            tot.exposure = sum(exposure))

######################################
#Waning of immunity - heterologuous ####
######################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",
                                         0.7*(1-pgamma(intro, shape = 4,rate =0.1/7)),0)}
Tdist = Tdist.vac.wane.ds.pas
fadeout_func = fadeout_func.wane.ds.pas
exposure.function <- exposure.function.wane.ds.pas

#sims####
srm <- Sys.time()

source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_wane_ds_lpda.RDS")
saveRDS(histories, file = "./output/histories_wane_ds_lpda.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_wane_ds_lpda.RDS")
saveRDS(Intro_matrix, file = "./output/Intro_matrix_wane_ds_lpda.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;joined.intro =NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,
                                                   spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
  joined.intro <- rbind(joined.intro,cbind(spatial.input$host_id,
                                           spatial.input$TYPE,spatial.input$SIZE, Intro_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
names(joined.intro)<- c("host_id","TYPE","SIZE","intro.time","run")
saveRDS(joined.histories, file = "./output/joined_histories_wane_ds_lpda.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_wane_ds_lpda.RDS")
saveRDS(joined.intro,file = "./output/intro_wane_ds_lpda.RDS")
Sys.time()-srm
names(infected_over_time_runs)<- c("Time","I","C", "PC", "run")
#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_wane_ds_lpda.RDS")
histories<- readRDS(file = "./output/histories_wane_ds_lpda.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_wane_ds_lpda.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_wane_ds_lpda.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_wane_ds_lpda.RDS")
joined.intro<- readRDS(file = "./output/intro_wane_ds_lpda.RDS")


#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.wane.ds<- left_join(infectious.periods,joined.vacstatus[,  c("host_id","run","vac")], by = c("host_id","run"))
inf.vac.wane.ds<- left_join(inf.vac.wane.ds, joined.intro[,  c("host_id","run","intro.time")], by = c("host_id","run"))
inf.vac.wane.ds<- left_join(inf.vac.wane.ds,spatial.input[,c("host_id","TYPE","SIZE")], by = c("host_id"))

#no fit 
inf.vac.wane.ds$exposure<- exposure.function(inf.vac.wane.ds) 
#replace negative numbers
inf.vac.wane.ds$exposure[inf.vac.wane.ds$exposure<100]<- 100
exposure.wane.ds.tot <- inf.vac.wane.ds%>%reframe(.by = run,
                                                  tot.exposure = sum(exposure))


############################################################
# Combine SCENARIOS                                        #
############################################################

#summary plots for lpda ###
infected_over_time_runs.scenarios<- rbind(cbind(readRDS(file = "./output/infected_over_time_runs_novaccination_ldpa.RDS"), scenario = "1. Baseline"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_clinprot_lpda.RDS"), scenario = "3. Clinical protection"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_60vaccination_lpda.RDS"), scenario = "2.Reduce Transmission"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_wane_lpda.RDS"), scenario = "4. Waning \n homologous"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_wane_ds_lpda.RDS"), scenario = "5. Waning \n heterologous"))

plot.data.scenarios <- infected_over_time_runs.scenarios%>%reframe(.by = c(run,scenario),
                                                                   Detected = max(C)+max(PC))%>%reshape2::melt(id.vars=c("run","scenario"))
ggplot(plot.data.scenarios, aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
  scale_fill_manual(values =c("Detected" = "red", "Screened" = "gray"), 
                    labels =c("Detected" = "Infected", "Screened" = "Screened"),
                    name = "Culling")+
  labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(.~scenario)+
  #ggtitle("Infected farms (LPDA)")+
  xlim(0,60)+
  theme(legend.position = "none")
ggsave("./output/figures/InfectedFarmslpda.png", scale =  1.23)

############################################################
# Combine human exposure with the no vaccination scenario   #
############################################################
lpda.exposure.total <- rbind(cbind(exposure.novaccin.tot, vac = 0, scenario = "1. Baseline"),
                             cbind(exposure.60vaccin.tot, vac = 0.6, scenario = "2. Reduced transmission"),
                             cbind(exposure.clinprot.tot, vac = 0.6, scenario = "3. Clinical protection"),
                             cbind(exposure.wane.tot, vac = 0.6, scenario = "4. Waning immunity \n homologous"),
                             cbind(exposure.wane.ds.tot, vac = 0.6, scenario = "5. Waning immunity \n heterologous"))
lpda.exposure.total$ratio <- lpda.exposure.total$tot.exposure/mean.exposure.novaccin

# Visualize Human exposure ####

ggplot(lpda.exposure.total, aes(ratio+10^-5))+ 
  geom_histogram(aes(y =after_stat(count/sum(count))),position = "identity",binwidth = 0.5)+
  labs(x = "log10 ratio of exposure",y = "Proportion of runs")+
  geom_vline(xintercept = 1)+
  scale_x_log10(n.breaks = 5,limits = c(10^-6,25))+
  facet_grid(.~scenario)+
#  ggtitle("Exposure during an outbreak of multiple farms")+
  theme(legend.position = "none")

ggsave("./output/figures/ExposureSpatial.lpda.png", scale = 1.23)






#############################################################
#   Active surveillance                                     #
#############################################################
index.farm <- hpda
max.runs <- length(index.farm)
####################
#no vaccination ####
####################
#set functions ####
vac.func <- function(type, intro){0}
Tdist = Tdist.novac.min
fadeout_func = fadeout_func.generic
exposure.function <- exposure.function.novac.min

#run the simulations ####
srm <- Sys.time()
set.seed(23011977)
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_novaccination_hdpa_min.RDS")
saveRDS(histories, file = "./output/histories_novaccination_hdpa_min.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_hdpa_min.RDS")
#join the histories
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joinedhistories_novaccination_hdpa_min.RDS")
saveRDS(joined.vacstatus, file = "./output/joinedvacstatus_novaccination_hdpa_min.RDS")
Sys.time()-srm

#read saved simulations####
infected_over_time_runs<- readRDS( file = "./output/infected_over_time_runs_novaccination_hdpa_min.RDS")
histories <- readRDS(file = "./output/histories_novaccination_hdpa_min.RDS")
V_stat_matrix <- readRDS(file = "./output/V_stat_matrix_hdpa_min.RDS")
joined.histories<- readRDS(file = "./output/joinedhistories_novaccination_hdpa_min.RDS")
joined.vacstatus<- readRDS( file = "./output/joinedvacstatus_novaccination_hdpa_min.RDS")




#length of infectious periods
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)
#add vaccination status, farm type and size 
inf.vac.novaccin<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))

#determine the exposure ###
inf.vac.novaccin$exposure <-exposure.function(inf.vac.novaccin)

inf.vac.novaccin$exposure[inf.vac.novaccin$exposure<0]<-500
exposure.novaccin.tot <- inf.vac.novaccin%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))

mean.exposure.novaccin.min <- mean(exposure.novaccin.tot$tot.exposure)
#############################################################################################
#50% high titre -> exposure per farm is slightly less, but detection times are increased ####
#############################################################################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",0.5,0)}
Tdist = Tdist.vac.min
fadeout_func = fadeout_func.generic
exposure.function <- exposure.function.vac.min

#run simulations ####
srm <- Sys.time()
set.seed(23011977)
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_60vaccination_hpda_min.RDS")
saveRDS(histories, file = "./output/histories_60vaccination_hpda_min.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_60vaccination_hpda_min.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,
                                                   spatial.input$TYPE,
                                                   spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joinedhistories_60vaccination_hpda_min.RDS")
saveRDS(joined.vacstatus, file = "./output/joinedvacstatus_60vaccination_hpda_min.RDS")

Sys.time()-srm

#read saved simulations
infected_over_time_runs<- readRDS( file = "./output/infected_over_time_runs_60vaccination_hpda_min.RDS")
histories <- readRDS(file = "./output/histories_60vaccination_hpda_min.RDS")
V_stat_matrix <- readRDS(file = "./output/V_stat_matrix_60vaccination_hpda_min.RDS")
joined.histories<- readRDS(file = "./output/joinedhistories_60vaccination_hpda_min.RDS")
joined.vacstatus<- readRDS( file = "./output/joinedvacstatus_60vaccination_hpda_min.RDS")

#calculate infectious period ####
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.60vac<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))

#determine exposure ####
inf.vac.60vac$size <- inf.vac.60vac$SIZE
inf.vac.60vac$exposure <- exposure.function(inf.vac.60vac)
exposure.60vaccin.tot <- inf.vac.60vac%>%reframe(.by = run,
                                                 tot.exposure = sum(exposure))

##################################################################################################################
#Clinical protection Vaccination has 50% high titre with reduced transmission rate and 50% death of low titre birds ####
#################################################################################################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",0.5,0)}
Tdist = Tdist.clinprot.min
fadeout_func = fadeout_func.generic
exposure.function <- exposure.function.clinprot.min

#sims###
srm <- Sys.time()
# vac.func <- function(x){ifelse(x == "LAYER",10^-6 #does not effect pmajor but changes the infectious period
#                                ,0)} 
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_clinprot_hpda_min.RDS")
saveRDS(histories, file = "./output/histories_clinprot_hpda_min.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_clinprot_hpda_min.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joined_histories_clinprot_hpda_min.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_clinprot_hpda_min.RDS")
Sys.time()-srm

#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_clinprot_hpda_min.RDS")
histories<- readRDS(file = "./output/histories_clinprot_hpda_min.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_clinprot_hpda_min.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_clinprot_hpda_min.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_clinprot_hpda_min.RDS")


#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.clinprot<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))


inf.vac.clinprot$exposure<- exposure.function(inf.vac.clinprot)
exposure.clinprot.tot <- inf.vac.clinprot%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))


######################################
#Waning of immunity - homologuous ####
######################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",
                                         1-pgamma(intro,36,0.07),0)}
Tdist = Tdist.vac.wane.min
fadeout_func =fadeout_func.generic
exposure.function <- exposure.function.wane.min

#sims####
srm <- Sys.time()
# vac.func <- function(type){ifelse(type == "LAYER",
#                                   1-pgamma(runif(1,min = 0,max = 540), shape = 36,rate =0.07),
#                                   0)}
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_wane_hpda_min.RDS")
saveRDS(histories, file = "./output/histories_wane_hpda_min.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_wane_hpda_min.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joined_histories_wane_hpda_min.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_wane_hpda_min.RDS")
Sys.time()-srm

#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_wane_hpda_min.RDS")
histories<- readRDS(file = "./output/histories_wane_hpda_min.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_wane_hpda_min.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_wane_hpda_min.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_wane_hpda_min.RDS")


#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.wane<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))
#add introduction times
intro.times.run.id <- NULL;
for(j in c(1:length(Intro_matrix[,1]))){
  intro.times.run.id<- rbind(intro.times.run.id,
                             data.frame(intro.time = Intro_matrix[j,],
                                        host_id = spatial.input$host_id,
                                        run = j))
  
}
inf.vac.wane<- left_join(inf.vac.wane,intro.times.run.id,by = c("host_id","run"))

#no fit 

inf.vac.wane$exposure<- exposure.function(inf.vac.wane) 
exposure.wane.tot <- inf.vac.wane%>%reframe(.by = run,
                                            tot.exposure = sum(exposure))

######################################
#Waning of immunity - heterologuous ####
######################################
#set functions ####
vac.func <- function(type, intro){ifelse(type == "LAYER",
                                         0.7*(1-pgamma(intro, shape = 4,rate =0.1/7)),0)}
Tdist = Tdist.vac.wane.ds.min
fadeout_func = fadeout_func.generic
exposure.function <- exposure.function.wane.ds.min

#sims####
srm <- Sys.time()

source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_wane_ds_hpda_min.RDS")
saveRDS(histories, file = "./output/histories_wane_ds_hpda_min.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_wane_ds_hpda_min.RDS")
saveRDS(Intro_matrix, file = "./output/Intro_matrix_wane_ds_hpda_min.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;joined.intro =NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$host_id,
                                                   spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
  joined.intro <- rbind(joined.intro,cbind(spatial.input$host_id,
                                           spatial.input$TYPE,spatial.input$SIZE, Intro_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
names(joined.intro)<- c("host_id","TYPE","SIZE","intro.time","run")
saveRDS(joined.histories, file = "./output/joined_histories_wane_ds_hpda_min.RDS")
saveRDS(joined.vacstatus,file = "./output/vacstatus_wane_ds_hpda_min.RDS")
saveRDS(joined.intro,file = "./output/intro_wane_ds_hpda_min.RDS")
Sys.time()-srm

#read saved simulations ####
infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_wane_ds_hpda_min.RDS")
histories<- readRDS(file = "./output/histories_wane_ds_hpda_min.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_wane_ds_hpda_min.RDS")
joined.histories<- readRDS( file = "./output/joined_histories_wane_ds_hpda_min.RDS")
joined.vacstatus <- readRDS(file = "./output/vacstatus_wane_ds_hpda_min.RDS")
joined.intro<- readRDS(file = "./output/intro_wane_ds_hpda_min.RDS")


#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.wane.ds<- left_join(infectious.periods,joined.vacstatus[,  c("host_id","run","vac")], by = c("host_id","run"))
inf.vac.wane.ds<- left_join(inf.vac.wane.ds, joined.intro[,  c("host_id","run","intro.time")], by = c("host_id","run"))
inf.vac.wane.ds<- left_join(inf.vac.wane.ds,spatial.input[,c("host_id","TYPE","SIZE")], by = c("host_id"))

#no fit 
inf.vac.wane.ds$exposure<- exposure.function(inf.vac.wane.ds) 
#replace negative numbers
inf.vac.wane.ds$exposure[inf.vac.wane.ds$exposure<100]<- 100
exposure.wane.ds.tot <- inf.vac.wane.ds%>%reframe(.by = run,
                                                  tot.exposure = sum(exposure))

############################################################
# Combine SCENARIOS                                        #
############################################################

#summary plots for HDPA ###
infected_over_time_runs.scenarios.min<- rbind(cbind(readRDS(file = "./output/infected_over_time_runs_novaccination_hdpa_min.RDS"), scenario = "1. Baseline"),
                                          cbind(readRDS( file = "./output/infected_over_time_runs_clinprot_hpda_min.RDS"), scenario = "3. Clinical protection"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_60vaccination_hpda_min.RDS"), scenario = "2. Reduce Transmission"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_wane_hpda_min.RDS"), scenario = "4. Waning \n homologous"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_wane_ds_hpda_min.RDS"), scenario = "5. Waning \n heterologous"))

plot.data.scenarios.min <- infected_over_time_runs.scenarios.min%>%reframe(.by = c(run,scenario),
                                                                   Detected = max(C)+max(PC))%>%reshape2::melt(id.vars=c("run","scenario"))
ggplot(plot.data.scenarios.min, aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
  scale_fill_manual(values =c("Detected" = "red", "Screened" = "gray"), 
                    labels =c("Detected" = "Infected", "Screened" = "Screened"),
                    name = "Culling")+
  labs(x = "Number of farms",y = "Proportion of runs")+
  facet_grid(.~scenario)+
  xlim(0,60)+
  theme(legend.position = "none")#+
#  ggtitle("Infected farms")
ggsave("./output/figures/InfectedFarmshpda_min.png")

############################################################
# Combine human exposure with the no vaccination scenario   #
############################################################
hpda.exposure.total.min <- rbind(cbind(exposure.novaccin.tot, vac = 0, scenario = "1. Baseline"),
                            cbind(exposure.60vaccin.tot, vac = 0.6, scenario = "2. Reduced transmission"),
                             cbind(exposure.clinprot.tot, vac = 0.6, scenario = "3. Clinical protection"),
                             cbind(exposure.wane.tot, vac = 0.6, scenario = "4. Waning immunity \n homologous"),
                             cbind(exposure.wane.ds.tot, vac = 0.6, scenario = "5. Waning immunity \n heterologous"))
hpda.exposure.total.min$ratio <- hpda.exposure.total.min$tot.exposure/mean.exposure.novaccin.min
hpda.exposure.total.min$ratio.baseline <- hpda.exposure.total.min$tot.exposure/mean.exposure.novaccin

# Visualize Human exposure ####

ggplot(hpda.exposure.total.min, aes(ratio))+ 
  geom_histogram(aes(y =after_stat(count/sum(count))),position = "identity",binwidth = 0.5)+
  labs(x = "log10 ratio of exposure",y = "Proportion of runs")+
  geom_vline(xintercept = 1)+
  scale_x_log10()+
  facet_grid(.~scenario)+
  #ggtitle("Exposure during an outbreak of multiple farms")+
  theme(legend.position = "none")
ggsave("./output/figures/ExposureSpatial_min.png", scale = 1.23)


ggplot(hpda.exposure.total.min, aes(ratio.baseline))+ 
  geom_histogram(aes(y =after_stat(count/sum(count))),position = "identity",binwidth = 0.5)+
  labs(x = "log10 ratio of exposure",y = "Proportion of runs")+
  geom_vline(xintercept = 1)+
  scale_x_log10()+
  facet_grid(.~scenario)+
  #ggtitle("Exposure during an outbreak of multiple farms")+
  theme(legend.position = "none")
ggsave("./output/figures/ExposureSpatial_min_base.png", scale = 1.23)





############################################
#Spatial plotting                       ####
############################################
require(sf)
netherlands.contours <- st_read( "./input/Nederland/Nederland")
plot(st_geometry(netherlands.contours))

ggplot()+
  geom_sf(data= netherlands.contours,colour = "black",alpha = 0.1)+
  geom_point(data= spatial.input,aes(X,Y), colour = "darkgrey")+
  geom_point(data= spatial.input[spatial.input$host_id%in%hpda,],aes(X,Y), colour = "red")+
  geom_point(data= spatial.input[spatial.input$host_id%in%lpda,],aes(X,Y), colour = "black")+theme_bw()



netherlands.contours$geometry
coords.sims <- joined.histories%>%select(x_coord,y_coord)%>%reframe(x = x_coord,y = y_coord)%>%as.matrix
SpatialMultiPointsDataFrame(coords = coords.sims,
                            data = joined.histories%>%select(host_id,Type_event,Event_time,run),
                            proj4string = CRS(netherlands.contours))

ggplot()+
  geom_sf(data= netherlands.contours,colour = "black")+
  geom_point(aes(X,Y), colour = "black",alpha = 0.5,spatial.input)+
  geom_point(aes(x_coord*1000,y_coord*1000, colour =  as.factor(Type_event)),data = joined.histories%>%select(host_id,x_coord,y_coord,Type_event,Event_time,run))+
  facet_grid(run~.)+theme_bw(legend.position = "none" )



