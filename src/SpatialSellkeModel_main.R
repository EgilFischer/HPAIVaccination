########## start the simulation #####################################
### Start the simulation
# I first call the R file where I initialize the variables ("InitSim.R")
# Then I call the file ("SimLoop.R") that calls in a loop the function Event() defined above

histories <- list()
infected_over_time_runs <- data.frame()
#ptm <- proc.time()
for(i in c(1:2000))
 {
K<-1943 # define the first infected
source("./src/SpatialSellkeModel_init.R") # initialization and setting of the first infected
T_inf <- T_inf_matrix[i,] # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
Q_init <- Q_init_matrix[i,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1

while(nrow(Queue)!=0){
  source("./src/SpatialSellkeModel_loop.R")
}
# Show the simulation outputs
# The outputs are stored in the matrix History
histories[i] <- History
infected_over_time_runs <- rbind(infected_over_time_runs,cbind(infected_over_time,data.frame(run = i)))

}
names(infected_over_time_runs )<- c("Time","I","run")
plot(infected_over_time, type= "l")
saveRDS(infected_over_time_runs, file = "./output/infected_over_time_runs.RDS")
saveRDS(histories, file = "./output/histories.RDS")

ggplot()+
  geom_point(aes(X,Y), data = spatial.input)+
  geom_point(aes(x_coord,y_coord), data = History%>%select(host_id,x_coord,y_coord)%>%unique, color = "red")+
  geom_point(aes(X,Y), data = spatial.input[K,], color = "lightblue")

