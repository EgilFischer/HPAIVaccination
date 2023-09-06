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

#function to set the infectious period distribution depending on the vaccination status
#load fitted distributions
rate.fit <- readRDS("./output/pasdetRate_pasdet.RDS")
rate.function <- function(vac){with(rate.fit,
                                    coefficients[1]+coefficients[2]*vac*100 + coefficients[3]*(vac*100)^2 +coefficients[4]*(vac*100)^3 )}

shape.fit <- readRDS("./output/pasdetShape_pasdet.RDS")
shape.function <- function(vac){with(shape.fit,
                                     coefficients[1]+coefficients[2]*vac*100 + coefficients[3]*(vac*100)^2 +coefficients[4]*(vac*100)^3 )}

#function to set the infectious period distribution depending on the vaccination status
Tdist = function(vacstat,type,fadeout){
  shape = shape.function(vacstat);
  rate = rate.function(vacstat);
  return(rgamma(1, shape= shape, rate = rate))
}

#generic fade out function always return 0
fadeout_func = function(vacstat){
  return(0);
}


#introduction locations
set.seed(23011977)
spatial.input$density1KM <- pointdensity(spatial.input[,c("Xkm","Ykm")], eps = 1) 

hpda <- sample(c(spatial.input%>%
                   filter(spatial.input$density1KM>= quantile(spatial.input$density1KM,0.95))%>%select("nr"))$nr, 100, replace = TRUE)
#spatial.input%>%filter(nr %in% hpda)

lpda <- sample(c(spatial.input%>%
                   filter(spatial.input$density1KM<= quantile(spatial.input$density1KM,0.25))%>%select("nr"))$nr, 100, replace = TRUE)
#spatial.input%>%filter(nr %in% lpda)


#Simulations ####
######################################################
#  High Poultry Density Area                         #
######################################################
index.farm <- hpda
max.runs <- length(index.farm)

#no vaccination ####
#set functions ####
vac.func <- function(vac,type){0}
set.seed(23011977)

#run the simulations ####
srm <- Sys.time()
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
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$nr,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joinedhistories_novaccination_hdpa.RDS")
saveRDS(joined.vacstatus, file = "./output/joinedvacstatus_novaccination_hdpa.RDS")
Sys.time()-srm

#read saved simulations
infected_over_time_runs<- readRDS( file = "./output/infected_over_time_runs_novaccination_hdpa.RDS")
histories <- readRDS(file = "./output/histories_novaccination_hdpa.RDS")
V_stat_matrix <- readRDS(file = "./output/V_stat_matrix_hdpa.RDS")
joined.histories<- readRDS(file = "./output/joinedhistories_novaccination_hdpa.RDS")
joined.vacstatus<- readRDS( file = "./output/joinedvacstatus_novaccination_hdpa.RDS")


#visualize  & summarize runds ####
ggplot(infected_over_time_runs%>%reframe(.by = run,
                                  Detected = max(C),
                                  Screen = max(PC))%>%reshape2::melt(id.vars="run"),aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
  scale_fill_manual(values =c("Detected" = "red", "Screen" = "gray"), 
                    labels =c("Detected" = "Detected", "Screen" = "Screen"),
                    name = "Culling")+
  labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(variable~.)
#ggsave("./output/figures/spatialHPDAculledfarms.png", scale =1.23)

#length of infectious periods
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)
#add vaccination status, farn type and size 
inf.vac.novaccin<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))

#determine the exposure ###
exposure.novaccin.fit <- readRDS(file = "./output/exposurefitLayerPas.RDS")


inf.vac.novaccin$exposure <- c(predict(exposure.novaccin.fit, 
                                       newdata =data.frame(size = inf.vac.novaccin$SIZE, detection.time =inf.vac.novaccin$Tinf , vaccination = c(0.0))))
inf.vac.novaccin$exposure[inf.vac.novaccin$exposure<0]<-500
exposure.novaccin.tot <- inf.vac.novaccin%>%reframe(.by = run,
                           tot.exposure = sum(exposure))

mean.exposure.novaccin <- mean(exposure.novaccin.tot$tot.exposure)

#60% high titre -> exposure per farm is slightly less, but detection times are increased ####
set.seed(23011977)
srm <- Sys.time()
#change the vaccination function
vac.func <- function(type){ifelse(type == "LAYER",0.6,0)}

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
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$nr,
                                                   spatial.input$TYPE,
                                                   spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")
saveRDS(joined.histories, file = "./output/joinedhistories_60vaccination_hdpa.RDS")
saveRDS(joined.vacstatus, file = "./output/joinedvacstatus_60vaccination_hdpa.RDS")

Sys.time()-srm

#read saved simulations
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
exposure.vaccin.fit <- readRDS( file = "./output/exposurefitLayerPasVaccin.RDS")
exposure.function.60vac <- function(size, vac){if(vac == 0)
{max(0,predict(exposure.novaccin.fit, newdata =data.frame(size = size)))}
  else{max(0,predict(exposure.vaccin.fit, newdata =data.frame(size = size, vaccination =vac)))}}

inf.vac.60vac$exposure <- mapply(exposure.function.60vac, inf.vac.60vac$SIZE, inf.vac.60vac$vac)
exposure.60vaccin.tot <- inf.vac.60vac%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))

#Clinical protection Vaccination has 50% high titre with reduced transmission rate and 50% death of low titre birds ####
rate.fit <- readRDS( "./output/pasdetRate_clinprot.RDS")
shape.fit <- readRDS( "./output/pasdetShape_clinprot.RDS")

Tdist = function(vacstat,type,fadeout){
  if(vacstat==0)
  {
    shape = shape.fit[1,2]
    rate = rate.fit[1,2]
  }else {
    #select scenario 3 = 50% clinical protection but no protection against transmission
    shape = shape.fit[3,2]
    rate = rate.fit[3,2]
  }
    return(rgamma(1, shape= shape, rate = rate))
}

#sims
srm <- Sys.time()
vac.func <- function(x){ifelse(x == "LAYER",10^-6 #does not effect pmajor but changes the infectious period
                               ,0)} 
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
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$nr,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
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

exposure.clinprot.fit <- readRDS( file = "./output/exposurefitLayerVaccinClinprot.RDS")

inf.vac.clinprot$exposure<- predict(exposure.clinprot.fit, newdata = data.frame(detection.time =inf.vac.clinprot$Tinf))
exposure.clinprot.tot <- inf.vac.clinprot%>%reframe(.by = run,
                                                   tot.exposure = sum(exposure))
# 
# exposure.tot <- rbind(cbind(exposure.novaccin.tot, scenario = "Baseline"),
#                        cbind(exposure.60vaccin.tot, scenario = "Reduced transmission"),
#                        cbind(exposure.clinprot.tot, scenario = "Clinical protection"))
# mean.exposure.novaccin <- mean(exposure.novaccin.tot$tot.exposure)
# exposure.tot$ratio <- exposure.tot$tot.exposure/mean.exposure.novaccin
# ggplot(exposure.tot)+
#   geom_histogram(aes(x = ratio), fill = "red", alpha = 0.5)+
#   geom_vline(aes(xintercept = 1))+
#   facet_grid(.~scenario)
# ggsave("./output/figures/spatialexposure.png")
# exposure.tot%>%reframe(.by = scenario,
#                        median = median(ratio),
#                        low = quantile(ratio,.25),
#                        high = quantile(ratio,.75),
#                        above1 = sum(ratio>1)/100   )






#Waning of immunity ####
rate.fit <- readRDS( "./output/pasdetRate_pasdet.RDS")
shape.fit <- readRDS( "./output/pasdetShape_pasdet.RDS")

Tdist = function(vacstat,type,fadeout){
  if(fadeout==1)return(5)
  shape = shape.function(vacstat);
  rate = rate.function(vacstat);
  if(rate <=0|shape<=0 ) return(0);
  return(rgamma(1, shape= shape, rate = rate))
}

fadeout_func = function(vacstat){
  if(vacstat> 0.9999)return(1)else return(0)
}


srm <- Sys.time()
vac.func <- function(type){ifelse(type == "LAYER",
                                  1-pgamma(runif(1,min = 0,max = 540), shape = 36,rate =0.07),
                                  0)}
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_wane_hdpa.RDS")
saveRDS(histories, file = "./output/histories_wane_hdpa.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_wane_hdpa.RDS")
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$nr,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
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


#length of infectious periods and exposure###
infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.wane<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))

exposure.wane.fit <- readRDS("C:/Surfdrive/Projecten/AI/LNVVaccinatie/HPAIVaccination/output/exposurefitLayerVaccinwaneDifStrainPAS.RDS")
#no fit 

inf.vac.wane$exposure<- predict(exposure.wane.fit, newdata = data.frame(detection.time =inf.vac.wane$Tinf,phigh = inf.vac.wane$vac))
exposure.wane.tot <- inf.vac.wane%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))

############################################################
# Combine SCENARIOS                                        #
############################################################

#summary plots for HDPA ###
infected_over_time_runs.scenarios<- rbind(cbind(readRDS(file = "./output/infected_over_time_runs_novaccination_hdpa.RDS"), scenario = "Baseline"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_clinprot_hdpa.RDS"), scenario = "Clinical protection"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_60vaccination_hdpa.RDS"), scenario = "Reduce Transmission"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_wane_hdpa.RDS"), scenario = "Waning"))

plot.data.scenarios <- infected_over_time_runs.scenarios%>%reframe(.by = c(run,scenario),
                                            Detected = max(C)+max(PC))%>%reshape2::melt(id.vars=c("run","scenario"))
ggplot(plot.data.scenarios, aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
                   scale_fill_manual(values =c("Detected" = "red", "Screened" = "gray"), 
                                     labels =c("Detected" = "Infected", "Screened" = "Screened"),
                                     name = "Culling")+
                   labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(.~scenario)+
  ggtitle("Infected farms")
ggsave("./output/figures/InfectedFarmsHPDA.png")

############################################################
# Combine human exposure with the no vaccination scenario   #
############################################################
hpda.exposure.total <- rbind(cbind(exposure.novaccin.tot, vac = 0, scenario = "Baseline"),
                             cbind(exposure.60vaccin.tot, vac = 0.6, scenario = "60% Vaccination"),
                             cbind(exposure.clinprot.tot, vac = 0.6, scenario = "Clinical protection"),
                             cbind(exposure.wane.tot, vac = 0.6, scenario = "Waning immunity"))
hpda.exposure.total$ratio <- exposure.total$tot.exposure/mean.exposure.novaccin

# Visualize Human exposure ####
histories.scenarios<- rbind(cbind(readRDS(file = "./output/histories_novaccination_hdpa.RDS")%>%set_names(seq_along(.)) %>% 
                                    enframe %>%
                                    unnest(), data.frame(scenario = "Baseline")),
                            cbind(readRDS(file = "./output/histories_60vaccination_hdpa.RDS")%>%set_names(seq_along(.)) %>% 
                                    enframe %>%
                                    unnest(), scenario = "Reduced transmission"),
                            cbind(readRDS(file = "./output/histories_clinprot_hdpa.RDS")%>%set_names(seq_along(.)) %>% 
                                    enframe %>%
                                    unnest(), scenario = "Clinical protection"),
                            cbind(readRDS(file = "./output/histories_wane_hdpa.RDS")%>%set_names(seq_along(.)) %>% 
                                    enframe %>%
                                    unnest(), scenario = "Waning Immunity"))
names(histories.scenarios)[1]<- "run"
V_stat_matrix.base<-  readRDS(file = "./output/V_stat_matrix_hdpa.RDS")*0
V_stat_matrix.vaccin  <-   readRDS(file = "./output/V_stat_matrix_hdpa.RDS")
V_stat_matrix.clinprot    <-readRDS(file = "./output/V_stat_matrix_clinprot_hdpa.RDS")
V_stat_matrix.wane      <-  readRDS(file = "./output/V_stat_matrix_wane_hdpa.RDS")
  

exposure.baseline.fit <- readRDS( file = "./output/exposurefitLayerPas.RDS")
exposure.clinprot.fit <- readRDS( file = "./output/exposurefitLayerVaccinClinprot.RDS")
exposure.vaccin.fit <- readRDS( file = "./output/exposurefitLayerVaccin.RDS")
exposure.wane.fit  <- readRDS( file = "./output/exposurefitLayerWanePas.RDS")

infectious.periods<- histories.scenarios%>%reframe(.by = c(run, host_id, scenario),
                                                                   Tinf = diff(Event_time),
                                                                   Culled = max(Type_event)==4,
                                                                   Infected = min(Type_event)==2)
V_stat_vector.wane <- mapply(function(r,i){ V_stat_matrix.wane[as.numeric(r),as.numeric(i)]}, infectious.periods[infectious.periods$scenario == "Waning Immunity",]$run, infectious.periods[infectious.periods$scenario == "Waning Immunity",]$host_id)
V_stat_vector.vaccin <- mapply(function(r,i){ V_stat_matrix.vaccin[as.numeric(r),as.numeric(i)]}, infectious.periods[infectious.periods$scenario == "Reduced transmission",]$run, infectious.periods[infectious.periods$scenario == "Reduced transmission",]$host_id)                              

infectious.periods$exposure <-0
infectious.periods$exposure[infectious.periods$scenario =="Baseline"]<- predict(exposure.baseline.fit, newdata =data.frame(size= spatial.input$SIZE[infectious.periods$host_id[infectious.periods$scenario =="Baseline"]]))
infectious.periods$exposure[infectious.periods$scenario =="Reduced transmission"]<- predict(exposure.vaccin.fit, newdata =data.frame(size= spatial.input$SIZE[infectious.periods$host_id[infectious.periods$scenario =="Reduced transmission"]],
                                                                                                                                 vaccination = V_stat_vector.vaccin))
infectious.periods$exposure[infectious.periods$scenario =="Clinical protection"]<- predict(exposure.clinprot.fit, newdata =data.frame(detection.time = infectious.periods$Tinf[infectious.periods$scenario =="Clinical protection"]))
infectious.periods$exposure[infectious.periods$scenario =="Waning Immunity"]<- predict(exposure.baseline.fit, newdata =data.frame(size= spatial.input$SIZE[infectious.periods$host_id[infectious.periods$scenario =="Waning Immunity"]],
                                                                                                                                  vaccination =V_stat_vector.wane ))
infectious.periods[infectious.periods$exposure<0,]$exposure <-0
infectious.periods$exposure[infectious.periods$scenario =="Waning Immunity"& infectious.periods$Tinf==5]<-35 
exposure.by.run <-  infectious.periods%>%reframe(.by = c("scenario","run"),
                                                 exposure.tot = sum(exposure))
mean.exposure.baseline <- infectious.periods%>%filter(scenario == "Baseline")%>%reframe(mean(exposure))%>%unlist
exposure.by.run$ratio <- (exposure.by.run$exposure.tot)/mean.exposure.baseline

ggplot(exposure.by.run, aes(ratio))+ 
  geom_histogram(aes(y =after_stat(count/sum(count))),position = "identity",binwidth = 0.5)+
  labs(x = "log10 ratio of exposure",y = "Proportion of runs")+
  geom_vline(xintercept = 1)+
  scale_x_log10()+
  facet_grid(.~scenario)+
  ggtitle("Exposure during an outbreak of multiple farms")
ggsave("./output/figures/ExposureSpatial.png", scale = 1.23)






#low poultry density area####
#no vaccination ####
#function to set the infectious period distribution depending on the vaccination status
#load fitted distributions
rate.fit <- readRDS("./output/pasdetRate_pasdet.RDS")
rate.function <- function(vac){with(rate.fit,
                                    coefficients[1]+coefficients[2]*vac*100 + coefficients[3]*(vac*100)^2 +coefficients[4]*(vac*100)^3 )}

shape.fit <- readRDS("./output/pasdetShape_pasdet.RDS")
shape.function <- function(vac){with(shape.fit,
                                     coefficients[1]+coefficients[2]*vac*100 + coefficients[3]*(vac*100)^2 +coefficients[4]*(vac*100)^3 )}

#function to set the infectious period distribution depending on the vaccination status
Tdist = function(vacstat,type,fadeout){
  shape = shape.function(vacstat);
  rate = rate.function(vacstat);
  return(rgamma(1, shape= shape, rate = rate))
}

#generic fade out function always return 0
fadeout_func = function(vacstat){
  return(0);
}

srm <- Sys.time()
index.farm <- lpda
max.runs <- length(index.farm)
vac.func <- function(vac,type){0}
#{ifelse(x == "LAYER",runif(1,min = 0.5, max = 1.0),0)}
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_novaccination_lpda.RDS")
saveRDS(histories, file = "./output/histories_novaccination_lpda.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_lpda.RDS")
Sys.time()-srm


infected_over_time_runs<- readRDS( file = "./output/infected_over_time_runs_novaccination_lpda.RDS")
histories <- readRDS(file = "./output/histories_novaccination_lpda.RDS")

#visualize  & summarize

# ggplot(infected_over_time_runs)+
#   geom_path(aes(Time,I,colour = as.factor(run)))+
#   theme(legend.position = "none")

ggplot(infected_over_time_runs%>%reframe(.by = run,
                                         Detected = max(C),
                                         Screen = max(PC))%>%reshape2::melt(id.vars="run"),aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
  scale_fill_manual(values =c("Detected" = "red", "Screen" = "gray"), 
                    labels =c("Detected" = "Detected", "Screen" = "Screen"),
                    name = "Culling")+
  labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(variable~.)
ggsave("./output/figures/spatialLPDAculledfarms.png", scale =1.23)

#length of infectious periods
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$nr,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")

infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.novaccin<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))


exposure.novaccin.fit <- readRDS(file = "./output/exposurefitLayerPas.RDS")


inf.vac.novaccin$exposure <- c(predict(exposure.novaccin.fit, 
                                       newdata =data.frame(size = inf.vac.novaccin$SIZE, detection.time =inf.vac.novaccin$Tinf , vaccination = c(0.0))))
inf.vac.novaccin$exposure[inf.vac.novaccin$exposure<0]<-500
exposure.novaccin.tot <- inf.vac.novaccin%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))



#60% high titre -> exposure per farm is slightly less, but detection times are increased ####
srm <- Sys.time()
index.farm <- lpda
max.runs <- length(index.farm)
vac.func <- function(type){ifelse(type == "LAYER",0.6,0)}
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_60vaccination_lpda.RDS")
saveRDS(histories, file = "./output/histories_60vaccination_lpda.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_lpda.RDS")
Sys.time()-srm

infected_over_time_runs<- readRDS(file = "./output/infected_over_time_runs_60vaccination_lpda.RDS")


#summarize
# ggplot(infected_over_time_runs%>%reframe(.by = run,
#                                          C = max(C),
#                                          PC = max(PC)))+
#   geom_histogram(aes(C,after_stat(count/sum(count))))+labs(x="Infected culled",y= "Proportion of runs")
# ggplot(infected_over_time_runs%>%reframe(.by = run,
#                                          C = max(C),
#                                          PC = max(PC)))+geom_histogram(aes(PC))

#visualize
#ggplot(infected_over_time_runs%>%filter(Time>0))+geom_path(aes(Time,I,colour = as.factor(run)))+theme(legend.position = "none")
ggplot(infected_over_time_runs%>%reframe(.by = run,
                                         Detected = max(C),
                                         Screened = max(PC))%>%reshape2::melt(id.vars="run"),aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
  scale_fill_manual(values =c("Detected" = "red", "Screened" = "gray"), 
                    labels =c("Detected" = "Detected", "Screened" = "Screened"),
                    name = "Culling")+
  labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(variable~.)+ggtitle("Vaccination 60% High Titre")
ggsave("./output/figures/spatialLPDAVaccin60culledfarms.png")


#length of infectious periods
histories <- readRDS(histories, file = "./output/histories_60vaccination_lpda.RDS")

joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$nr,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")

infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.60vac<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))


exposure.vaccin.fit <- readRDS( file = "./output/exposurefitLayerPasVaccin.RDS")
exposure.function.60vac <- function(size, vac){if(vac == 0)
{max(0,predict(exposure.novaccin.fit, newdata =data.frame(size = size)))}
  else{max(0,predict(exposure.vaccin.fit, newdata =data.frame(size = size, vaccination =vac)))}}

inf.vac.60vac$exposure <- mapply(exposure.function.60vac, inf.vac.60vac$SIZE, inf.vac.60vac$vac)
exposure.60vaccin.tot <- inf.vac.60vac%>%reframe(.by = run,
                                                 tot.exposure = sum(exposure))
mean.exposure.novaccin <- mean(exposure.novaccin.tot$tot.exposure)
mean.exposure.vaccin <- mean(exposure.60vaccin.tot$tot.exposure)
exposure.total <- rbind(cbind(exposure.novaccin.tot, vac = 0),cbind(exposure.60vaccin.tot, vac = 0.6))
exposure.total$ratio <- exposure.total$tot.exposure/mean.exposure.novaccin
ggplot(exposure.total)+geom_histogram(aes(ratio, fill = as.factor(vac)),
                                      binwidth = 0.5,position = 'identity',colour = "black",alpha = 0.5)+
  geom_vline(aes(xintercept = 1))+
  scale_fill_manual(values =c("0" = "red","0.6"= "black"), name = "High titre",labels =c("0" = "0%","0.6"= "60%") )+
  labs(x = "Ratio", y = "Count")
ggsave("./output/figures/exposure.spatial.vac60%.LPDA.png")



#Clinical protection Vaccination has 50% high titre with reduced transmission rate and 50% death of low titre birds ####
rate.fit <- readRDS( "./output/pasdetRate_clinprot.RDS")
shape.fit <- readRDS( "./output/pasdetShape_clinprot.RDS")

Tdist = function(vacstat,type,fadeout){
  if(vacstat==0)
  {
    shape = shape.fit[1,2]
    rate = rate.fit[1,2]
  }else {
    #select scenario 3 = 50% clinical protection but no protection against transmission
    shape = shape.fit[3,2]
    rate = rate.fit[3,2]
  }
  return(rgamma(1, shape= shape, rate = rate))
}
#sims
srm <- Sys.time()
index.farm <- lpda
max.runs <- length(index.farm)
vac.func <- function(x){ifelse(x == "LAYER",10^-6 #does not effect pmajor but changes the infectious period
                               ,0)} 

Tdist
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_clinprot_lpda.RDS")
saveRDS(histories, file = "./output/histories_clinprot_lpda.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_clinprot_lpda.RDS")
Sys.time()-srm

#summarize
# ggplot(infected_over_time_runs%>%reframe(.by = run,
#                                          C = max(C),
#                                          PC = max(PC)))+geom_histogram(aes(C,after_stat(count/sum(count))))+labs(x="Infected culled",y= "Proportion of runs")
# ggplot(infected_over_time_runs%>%reframe(.by = run,
#                                          C = max(C),
#                                          PC = max(PC)))+geom_histogram(aes(PC))

infected_over_time_runs<-readRDS( file = "./output/infected_over_time_runs_clinprot_lpda.RDS")
histories<- readRDS(file = "./output/histories_clinprot_lpda.RDS")
V_stat_matrix<- readRDS(file = "./output/V_stat_matrix_clinprot_lpda.RDS")

#visualize
ggplot(infected_over_time_runs%>%filter(Time>0))+geom_path(aes(Time,I,colour = as.factor(run)))+theme(legend.position = "none")

ggplot(infected_over_time_runs%>%reframe(.by = run,
                                         Detected = max(C),
                                         Preemptive = max(PC))%>%reshape2::melt(id.vars="run"),aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
  scale_fill_manual(values =c("Detected" = "red", "Preemptive" = "gray"), 
                    labels =c("Detected" = "Detected", "Preemptive" = "Screened"),
                    name = "Culling")+
  labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(variable~., labeller = labeller(variable =c("Detected" = "Detected", "Preemptive" = "Screened")))
ggsave("./output/figures/spatialHPDAVaccinClinProtculledfarmsLPDA.png")
#length of infectious periods
joined.histories <- NULL;joined.vacstatus<- NULL;
for(it in c(1:length(histories))){
  joined.histories<- rbind(joined.histories, cbind(histories[[it]], data.frame(run = it)))
  joined.vacstatus <- rbind(joined.vacstatus,cbind(spatial.input$nr,spatial.input$TYPE,spatial.input$SIZE, V_stat_matrix[it,], data.frame(run = it)))
}
names(joined.vacstatus) <- c("host_id","TYPE","SIZE","vac","run")

infectious.periods<- joined.histories%>%reframe(.by = c(run, host_id),
                                                Tinf = diff(Event_time),
                                                Culled = max(Type_event)==4,
                                                Infected = min(Type_event)==2)

inf.vac.clinprot<- left_join(infectious.periods,joined.vacstatus,by = c("host_id","run"))

exposure.clinprot.fit <- readRDS( file = "./output/exposurefitLayerVaccinClinprot.RDS")

inf.vac.clinprot$exposure<- predict(exposure.clinprot.fit, newdata = data.frame(detection.time =inf.vac.clinprot$Tinf))
exposure.clinprot.tot <- inf.vac.clinprot%>%reframe(.by = run,
                                                    tot.exposure = sum(exposure))

exposure.tot <- rbind(cbind(exposure.novaccin.tot, scenario = "Baseline"),
                      cbind(exposure.60vaccin.tot, scenario = "Reduced transmission"),
                      cbind(exposure.clinprot.tot, scenario = "Clinical protection"))
mean.exposure.novaccin <- mean(exposure.novaccin.tot$tot.exposure)
exposure.tot$ratio <- exposure.tot$tot.exposure/mean.exposure.novaccin
ggplot(exposure.tot)+
  geom_histogram(aes(x = ratio), fill = "red", alpha = 0.5)+
  geom_vline(aes(xintercept = 1))+
  facet_grid(.~scenario)
ggsave("./output/figures/spatialexposureLPDA.png")
exposure.tot%>%reframe(.by = scenario,
                       median = median(ratio),
                       low = quantile(ratio,.25),
                       high = quantile(ratio,.75),
                       above1 = sum(ratio>1)/100   )





#Waning of immunity ####
rate.fit <- readRDS( "./output/pasdetRate_pasdet.RDS")
shape.fit <- readRDS( "./output/pasdetShape_pasdet.RDS")

Tdist = function(vacstat,type,fadeout){
  if(fadeout==1)return(5)
  shape = shape.function(vacstat);
  rate = rate.function(vacstat);
  if(rate <=0|shape<=0 ) return(0);
  return(rgamma(1, shape= shape, rate = rate))
}

fadeout_func = function(vacstat){
  if(vacstat> 0.80)return(1)else return(0)
}


srm <- Sys.time()
index.farm <- lpda
max.runs <- length(index.farm)
vac.func <- function(type){ifelse(type == "LAYER",
                                  1-pgamma(runif(1,min = 0,max = 540), shape = 36,rate =0.07),
                                  0)}
source("./src/SpatialSellkeModel_init_simuls.R")
#run the model
source("./src/SpatialSellkeModel_main.R")
#save output
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs_wane_lpda.RDS")
saveRDS(histories, file = "./output/histories_wane_lpda.RDS")
saveRDS(V_stat_matrix,file = "./output/V_stat_matrix_wane_lpda.RDS")
Sys.time()-srm

infected_over_time_runs<- readRDS(file = "./output/infected_over_time_runs_wane_lpda.RDS")

ggplot(data = infected_over_time_runs)+
  geom_path(aes(Time, I,group = run), colour = "red")
ggplot(data = infected_over_time_runs)+
  geom_path(aes(Time, C,group = run), colour = "darkgreen")

ggplot(infected_over_time_runs%>%reframe(.by = run,
                                         Detected = max(C),
                                         Screened = max(PC))%>%reshape2::melt(id.vars="run"),aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
  scale_fill_manual(values =c("Detected" = "red", "Screened" = "gray"), 
                    labels =c("Detected" = "Detected", "Screened" = "Screened"),
                    name = "Culling")+
  labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(variable~.)+ggtitle("Vaccination waning")
ggsave("./output/figures/spatialHPDAVaccinWane.png")

#summary plots for lpda ###
infected_over_time_runs.scenarios<- rbind(cbind(readRDS(file = "./output/infected_over_time_runs_novaccination_lpda.RDS"), scenario = "Baseline"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_clinprot_lpda.RDS"), scenario = "Clinical protection"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_60vaccination_lpda.RDS"), scenario = "Reduce Transmission"),
                                          cbind(readRDS(file = "./output/infected_over_time_runs_wane_lpda.RDS"), scenario = "Waning"))

plot.data.scenarios <- infected_over_time_runs.scenarios%>%reframe(.by = c(run,scenario),
                                                                   Detected = max(C)+max(PC))%>%reshape2::melt(id.vars=c("run","scenario"))
ggplot(plot.data.scenarios, aes(value))+
  geom_histogram(aes(y =after_stat(count/sum(count)),fill = variable),position = "identity",binwidth = 1)+
  scale_fill_manual(values =c("Detected" = "red", "Screened" = "gray"), 
                    labels =c("Detected" = "Infected", "Screened" = "Screened"),
                    name = "Culling")+
  labs(x = "Number of farms",y = "Proportion of runs")+facet_grid(.~scenario)+
  ggtitle("Infected farms")
ggsave("./output/figures/InfectedFarmsLPDA.png")
# Human exposure ####
histories.scenarios<- rbind(cbind(readRDS(file = "./output/histories_novaccination_lpda.RDS")%>%set_names(seq_along(.)) %>% 
                                    enframe %>%
                                    unnest(), data.frame(scenario = "Baseline")),
                            cbind(readRDS(file = "./output/histories_60vaccination_lpda.RDS")%>%set_names(seq_along(.)) %>% 
                                    enframe %>%
                                    unnest(), scenario = "Reduced transmission"),
                            cbind(readRDS(file = "./output/histories_clinprot_lpda.RDS")%>%set_names(seq_along(.)) %>% 
                                    enframe %>%
                                    unnest(), scenario = "Clinical protection"),
                            cbind(readRDS(file = "./output/histories_wane_lpda.RDS")%>%set_names(seq_along(.)) %>% 
                                    enframe %>%
                                    unnest(), scenario = "Waning Immunity"))
names(histories.scenarios)[1]<- "run"
V_stat_matrix.base<-  readRDS(file = "./output/V_stat_matrix_lpda.RDS")*0
V_stat_matrix.vaccin  <-   readRDS(file = "./output/V_stat_matrix_lpda.RDS")
V_stat_matrix.clinprot    <-readRDS(file = "./output/V_stat_matrix_clinprot_lpda.RDS")
V_stat_matrix.wane      <-  readRDS(file = "./output/V_stat_matrix_wane_lpda.RDS")

exposure.baseline.fit <- readRDS( file = "./output/exposurefitLayerPas.RDS")
exposure.clinprot.fit <- readRDS( file = "./output/exposurefitLayerVaccinClinprotPAS.RDS")
exposure.vaccin.fit <- readRDS( file = "./output/exposurefitLayerVaccinPAS.RDS")
exposure.wane.fit  <- readRDS( file = "./output/exposurefitLayerWanePas.RDS")

infectious.periods<- histories.scenarios%>%reframe(.by = c(run, host_id, scenario),
                                                   Tinf = diff(Event_time),
                                                   Culled = max(Type_event)==4,
                                                   Infected = min(Type_event)==2)
V_stat_vector.wane <- mapply(function(r,i){ V_stat_matrix.wane[as.numeric(r),as.numeric(i)]}, infectious.periods[infectious.periods$scenario == "Waning Immunity",]$run, infectious.periods[infectious.periods$scenario == "Waning Immunity",]$host_id)
V_stat_vector.vaccin <- mapply(function(r,i){ V_stat_matrix.vaccin[as.numeric(r),as.numeric(i)]}, infectious.periods[infectious.periods$scenario == "Reduced transmission",]$run, infectious.periods[infectious.periods$scenario == "Reduced transmission",]$host_id)                              

infectious.periods$exposure <-0
infectious.periods$exposure[infectious.periods$scenario =="Baseline"]<- predict(exposure.baseline.fit, newdata =data.frame(size= spatial.input$SIZE[infectious.periods$host_id[infectious.periods$scenario =="Baseline"]]))
infectious.periods$exposure[infectious.periods$scenario =="Reduced transmission"]<- predict(exposure.vaccin.fit, newdata =data.frame(size= spatial.input$SIZE[infectious.periods$host_id[infectious.periods$scenario =="Reduced transmission"]],
                                                                                                                                     vaccination = V_stat_vector.vaccin))
infectious.periods$exposure[infectious.periods$scenario =="Clinical protection"]<- predict(exposure.clinprot.fit, newdata =data.frame(detection.time = infectious.periods$Tinf[infectious.periods$scenario =="Clinical protection"]))
infectious.periods$exposure[infectious.periods$scenario =="Waning Immunity"]<- predict(exposure.baseline.fit, newdata =data.frame(size= spatial.input$SIZE[infectious.periods$host_id[infectious.periods$scenario =="Waning Immunity"]],
                                                                                                                                  vaccination =V_stat_vector.wane))
                                                                                                                                                                                                                                                               
#fade out in waning immunity scenario -> 35 = approx 1.13*10*3
infectious.periods$exposure[infectious.periods$scenario =="Waning Immunity"& infectious.periods$Tinf==5]<-35 
infectious.periods[infectious.periods$exposure<0,]$exposure <-0

exposure.by.run <-  infectious.periods%>%reframe(.by = c("scenario","run"),
                                                 infected.farms = length(exposure),
                                                 exposure.tot = sum(exposure))
mean.exposure.baseline <- infectious.periods%>%filter(scenario == "Baseline")%>%reframe(mean(exposure))%>%unlist
exposure.by.run$ratio <- (exposure.by.run$exposure.tot)/mean.exposure.baseline

ggplot(exposure.by.run, aes(ratio))+ 
  geom_histogram(aes(y =after_stat(count/sum(count))),position = "identity",binwidth = 0.5)+
  labs(x = "log10 ratio of exposure",y = "Proportion of runs")+
  geom_vline(xintercept = 1)+
  scale_x_log10()+
  facet_grid(.~scenario)+
  ggtitle("Exposure during an outbreak of multiple farms")
ggsave("./output/figures/ExposureSpatialLPDA.png", scale = 1.23)




#Spatial plotting####
require(sf)
netherlands.contours <- st_read( "./input/Nederland/Nederland")
plot(st_geometry(netherlands.contours))

ggplot()+
  geom_sf(data= netherlands.contours,colour = "black",alpha = 0.1)+
  geom_point(data= spatial.input,aes(X,Y), colour = "darkgrey")+
  geom_point(data= spatial.input[spatial.input$nr%in%hpda,],aes(X,Y), colour = "red")+
  geom_point(data= spatial.input[spatial.input$nr%in%lpda,],aes(X,Y), colour = "black")+theme_bw()



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



