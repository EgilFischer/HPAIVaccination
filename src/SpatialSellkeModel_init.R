#########################################################
#                                                        
# Spatial SIR transmission model with the Sellke construction
# The model simulates infections spread between immobile hosts (e.g. poultry, farms, plants) using the Sellke construction.                         
# The model assumes only removal of the detected hosts (i.e. no preventive removal)                                                        
#                  Author:  Elisa Benincá    -  adapted by Egil Fischer                           
#                  Creation date  https://doi.org/10.1371/journal.pcbi.1008009
#########################################################


# Set input point pattern
matrix_points <- spatial.input[,c("X","Y")] # I select one point pattern as example
totpoints <- nrow(matrix_points) # total number of points
colnames(matrix_points) <- c("xcoord","ycoord")
# Add a column for the index
Index_points <- c(1:totpoints)
# Calculate the matrix of distances between points
# Express the coordinates as complex numbers for fast calculation of the euclidean distance
Coord <-
  (complex(length.out=2,real=matrix_points$xcoord,imaginary=matrix_points$ycoord));
distancematrix <- as.matrix(abs(outer(Coord,Coord,"-")))

#define a matrix for preemptive culling
culling.radius <- 3
culling.delay <- 1
cullingmatrix <- 1*(distancematrix<= culling.radius) #(times one keeps the matrix as it is and transforms boolean to numeric)
diag(cullingmatrix) <- matrix(0,nrow=totpoints); #no culling because of detection by itself
  
# Define the transmission kernel and calculate the hazard matrix
# Rescale the parameters h0 used in Boender et al. to account for the change in size and number of farms (see main text)
h0 <- 0.002*5360/totpoints;
alpha <- 2.1;
r0 <- 1.9;
h_kernel <- function(r){h0/(1 + (r/r0)^alpha)} ; # transmission kernel as a function of r
beta<-1;

# Create an hazard matrix evaluating for each host j the chance to be infected by host i as a function of distance
hazardmatrix <- as.matrix(apply(distancematrix,MARGIN=c(1,2),FUN=h_kernel));
#discount hazard by vaccination
diag(hazardmatrix) <- matrix(0,nrow=totpoints); # because the chance of infecting itself is 0






# code of InitSim.R sourced in the main code
# initialize vector status with 1 (all susceptibles)
Status <- matrix(1,nrow=totpoints) # status is a vector recording the state of each host
colnames(Status) <- "status" #1= susceptible, 2=infectious, 3 culled
Queue <- {}
History <- {}
Current <- {}
infected_over_time <- {}
time_vector <- {} # Initialize the Cumulative Force of Infection at the  beginning of the epidemic
CFI <- matrix(0,nrow=totpoints)
CFI_matrix <- matrix(0,ncol=totpoints,nrow=10000)
index_new_event <- 0
counter <- 0
# create vectors for the infected and the susceptibles, to use to compute the  cumulative force of infection
# These are index used for the calculation of the cumulative force of infection.
# They are used to keep track which hosts are infected or susceptible.
List_to_remove <- {}
List_to_infect <- {}
indexI <- matrix(0,nrow=totpoints) #indexI==0, not yet infected, indexI==1 infected, indexI==3 culled
indexS <- matrix(0,nrow=totpoints) # indexS==0, not yet infected, indexS==1 infected, indexS=3 culled
t_infection <- matrix(0,nrow=totpoints)
# coefficient of increase
bb <- matrix(0,nrow=totpoints)
tt <- 0 # start at time 0
next_infection_time <- 0
next_infection_host <- 0
####### # initialize with the first one to be infected
firstone <- K
Status[firstone] <- 2 # now Status= 2 (infected)
indexI[firstone] <- 1 # infected
indexS[firstone] <- 1 # not susceptible anymore
#calculate coefficient of increase
bb[which(indexS==0)] <-
  beta*apply(matrix(hazardmatrix[which(indexI==1),which(indexS==0)],nrow=length(which(indexI==1)),ncol=length(which(indexS==0))), MARGIN=2,FUN=sum) # In the rows (i) the infected. in the j the supectibles
t_infection[which(indexS==0)] <- (Q_init[which(indexS==0)]-  CFI[which(indexS==0)])/bb[which(indexS==0)]
t_infection[which(indexS==1)] <- 10000000 # I set an extremely high number,because it cannot infect itself
next_infection_time <- min(t_infection)
next_infection_host <- which.min(t_infection) # add this infected host to the history vector
History <-rbind(History,data.frame(Event_time=tt,Type_event=2,host_id=firstone,x_coord=Re(Coord[as.numeric(firstone)]),y_coord=Im(Coord[as.numeric(firstone)]))) # record it in the history vector. In the history vector add the coord
# update the list_to_infect and the list_to remove
List_to_infect <- rbind(List_to_infect,data.frame(Event_time=tt
                                                  +next_infection_time,Type_event=2,id_host = next_infection_host))
List_to_infect <- List_to_infect[order(List_to_infect[,1]),]
List_to_remove <- rbind(List_to_remove,data.frame(Event_time=tt
                                                  +T_inf[firstone],Type_event=3,id_host=firstone))
List_to_remove <- List_to_remove[order(List_to_remove[,1]),]
next_events <- rbind(List_to_infect[1,],List_to_remove[1,])
index_next_event <- which.min(next_events[,1])
Queue <- rbind(Queue,next_events[index_next_event,])
# now remove this event from the list_to_infect or list_to_remove
if(all.equal(cbind(List_to_infect[1,1],List_to_infect[1,2],List_to_infect[1,3]),cbind(Queue[1,1],Queue[1,2],Queue[1,3]))==TRUE){
  List_to_infect <- List_to_infect[-c(1),]
} else
  if(all.equal(cbind(List_to_remove[1,1],List_to_remove[1,2],List_to_remove[1,3]),cbind(Queue[1,1],Queue[1,2],Queue[1,3]))==TRUE){
    List_to_remove <- List_to_remove[-c(1),]
  }
# I add the current CFI (which is 0 at the beginning of the CFI_matrix)
index_new_event <- index_new_event+1
CFI_matrix[index_new_event,] <- CFI
timevector <- rbind(time_vector,tt)
infected_over_time <- rbind(infected_over_time,c(0,length(indexI[indexI==1])))

