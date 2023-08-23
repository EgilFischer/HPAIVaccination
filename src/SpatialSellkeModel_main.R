########## start the simulation #####################################
### Start the simulation
# I first call the R file where I initialize the variables ("InitSim.R")
# Then I call the file ("SimLoop.R") that calls in a loop the function Event() defined above
set.seed(1234)


source("./src/loadLibraries.R")
source("./src/SpatialSellkeModel_event.R")

histories <- list()
infected_over_time_runs <- data.frame()
if(!exists("max.runs")){max.runs <- 100}
ptm <- proc.time()
for(repetition in c(1:max.runs))
 {
  print(paste("Repetition ", repetition,"of",max.runs));
  K<- index.farm[repetition]  # define the first infected
  F_vector <- F_matrix[repetition,]#fade out = 1
  T_inf <- T_inf_matrix[repetition,] # rgamma(totpoints,10, scale=7/10) # mean=7, std=2
  Q_init <- Q_init_matrix[repetition,] # rexp(totpoints, rate = 1) # these are the thresholds (exposure to infection) picked from an exponential distribution of mean 1
  hazardmatrix <- t(t(hazardmatrix)*P_maj_matrix[repetition,])

  source("./src/SpatialSellkeModel_init.R") # initialization and setting of the first infected
  while(nrow(Queue)!=0){
  source("./src/SpatialSellkeModel_loop.R")
  }
  
  # The outputs are stored in the matrix History
  histories[[repetition]] <- History
  infected_over_time_runs <- rbind(infected_over_time_runs,cbind(infected_over_time,data.frame(run = repetition)))
}
names(infected_over_time_runs )<- c("Time","I","C","PC", "run")
