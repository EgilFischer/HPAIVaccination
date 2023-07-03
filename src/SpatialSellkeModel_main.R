########## start the simulation #####################################
### Start the simulation
# I first call the R file where I initialize the variables ("InitSim.R")
# Then I call the file ("SimLoop.R") that calls in a loop the function Event() defined above
source("./src/SpatialScenarios.R")
histories <- list()
infected_over_time_runs <- data.frame()
ptm <- proc.time()
for(repetition in c(1:10))
 {
K<- 1943;#sample(c(1:totpoints),1) # define the first infected
source("./src/SpatialSellkeModel_init.R") # initialization and setting of the first infected
T_inf <- T_inf_matrix[repetition,] # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
Q_init <- Q_init_matrix[repetition,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1

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

ggplot()+
  geom_point(aes(X,Y), data = spatial.input)+
  geom_point(aes(x_coord,y_coord), data = History%>%select(host_id,x_coord,y_coord)%>%unique, color = "red")+
  geom_point(aes(X,Y), data = spatial.input[K,], color = "lightblue")

