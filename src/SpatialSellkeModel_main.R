########## start the simulation #####################################
### Start the simulation
# I first call the R file where I initialize the variables ("InitSim.R")
# Then I call the file ("SimLoop.R") that calls in a loop the function Event() defined above
set.seed(1234)

param.list.baseline.layer <- list(
  scenario = "baseline_Layer", #scenario
  runs = 10, #number of runs
  max.time = 17*30,#length of the run
  itypes = 2, #type
  N0 = 45000, #population size - agramatie website
  initial= 10 , #initially infected - choosen value
  p.hightitre = 0,#proportion initially protected by vaccination
  beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = 2),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976 and Gemeraard et al 2023 #Type 1  = not protected by vaccination and type 2 = protected by vaccination
  infectious.period = c(3.0,4.0),#Duration infectious period 
  variance.infectious.period = c(3.0,4.0)^2/20, #Variance infectious period - Hobbelen et al uses shape parameter of 20 -> variance = m^2 / shape 
  transRate = matrix(c(0,0.0,0.0,0), nrow = 2), #value based on https://nvaidya.sdsu.edu/SIAP2015.pdf #transition rates should be of size itypes x itypes with 0 on diagonal
  pdie = c(0.99,0.01),#probability of dying at end of infectious period
  mortRate = 0.0005/7 #per capita death rate #Mortality events - based on performance reports and pers. com. mieke matthijs 0.5% per week +Gonzales & elbers figure 1 
)

source("./src/loadLibraries.R")
source("./src/SpatialScenarios.R")
source("./src/SpatialSellkeModel_init_simuls.R")
source("./src/SpatialSellkeModel_event.R")

histories <- list()
infected_over_time_runs <- data.frame()
ptm <- proc.time()
for(repetition in c(1:100))
 {
print(repetition);
K<- 1943;#sample(c(1:totpoints),1) # define the first infected
T_inf <- T_inf_matrix[repetition,] # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
Q_init <- Q_init_matrix[repetition,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
hazardmatrix <- t(t(hazardmatrix)*P_maj_matrix[repetition,])

source("./src/SpatialSellkeModel_init.R") # initialization and setting of the first infected
while(nrow(Queue)!=0){
  source("./src/SpatialSellkeModel_loop.R")
}
# Show the simulation outputs
# The outputs are stored in the matrix History
histories[[repetition]] <- History
infected_over_time_runs <- rbind(infected_over_time_runs,cbind(infected_over_time,data.frame(run = repetition)))

}
names(infected_over_time_runs )<- c("Time","I","run")
ggplot(infected_over_time_runs)+geom_path(aes(Time,I,colour = as.factor(run)))+theme(legend.position = "none")
#saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs.RDS")
#saveRDS(histories, file = "./output/histories.RDS")


ggplot(aes(x_coord,y_coord, colour =  as.factor(Type_event)),data = History%>%select(host_id,x_coord,y_coord,Type_event,Event_time))+
  geom_point(aes(X,Y), colour = "black",spatial.input)+
    geom_point()+ 
  facet_grid((Event_time< 10)~.)
                             


    
