# check the simulation of BI 
source("./src/loadLibraries.R")
#source("./src/pre_processing_detModule.R")
#source("./src/detectionModule.R")

#load data of BI and BI booster layer farm
BIBOOSTER_data <- readRDS("C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_output/2025_07_24_BIBOOSTER_layer_farm.RDS")
BI_data <- readRDS("C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_output/2025_07_25_BI_layer_farm.RDS") 
#load daily data and select BI and BIBOOSTER
BIBOOSTER_daily_data <- readRDS("C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_output/DailyData/DailyData.2025_07_24_BIBOOSTER_layer_farm.RDS")
BI_daily_data <- readRDS("C:/Users/fisch106/surfdrive/Projecten/AI/PPS_vaccinatie/Code/HPAI_vaccination/HPC_output/DailyData/DailyData.2025_07_25_BI_layer_farm.RDS")

#plot curves 
ggplot(BIBOOSTER_daily_data$out)+
  geom_point(aes(x=time, y=S.1/N),colour = "red")+
  geom_point(aes(x=time, y=S.2/N))

ggplot(BIBOOSTER_daily_data$out%>%filter(time <20))+
#  geom_point(aes(x=time, y=I.1/N),colour = "red")+
  geom_path(aes(x=time, y=I.2, colour = as.factor(run)))

ggplot(BI_daily_data$out%>%filter(time <20))+
  #  geom_point(aes(x=time, y=I.1/N),colour = "red")+
  geom_path(aes(x=time, y=I.2, colour = as.factor(run)))

ggplot(as.data.frame(BI_data$out)%>%filter(time <20))+
  #  geom_point(aes(x=time, y=I.1/N),colour = "red")+
  geom_path(aes(x=time, y=I.2, colour = as.factor(run)))


#detection module
ggplot(BI_daily_data$out%>%filter(run == 18))+geom_path(aes(x = time, y = I.1+I.2))
ggplot(BI_daily_data$out%>%filter(run == 18))+geom_path(aes(x = time, y = DR.2))

S1.results%>%filter(run ==18, scenario == "BI_layer_farm")
