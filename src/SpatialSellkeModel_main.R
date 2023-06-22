########## start the simulation #####################################
### Start the simulation
# I first call the R file where I initialize the variables ("InitSim.R")
# Then I call the file ("SimLoop.R") that calls in a loop the function Event() defined above


#ptm <- proc.time()
K<-1943 # define the first infected
T_inf <- T_inf_matrix[K,] # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
Q_init <- Q_init_matrix[K,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1

source("./src/SpatialSellkeModel_init.R") # initialization and setting of the first infected
while(nrow(Queue)!=0){
  source("./src/SpatialSellkeModel_loop.R")
}
# Show the simulation outputs
# The outputs are stored in the matrix History
head(History)
tail(History)
dim(History)

plot(infected_over_time, type= "l")

ggplot()+
  geom_point(aes(X,Y), data = spatial.input)+
  geom_point(aes(x_coord,y_coord), data = History%>%select(host_id,x_coord,y_coord)%>%unique, color = "red")+
  geom_point(aes(X,Y), data = spatial.input[K,], color = "lightblue")

