# code of SimLoop.R sourced in the main code
########### main loop between-herd model#############################
Current <- Queue[1,]
Queue <- Queue[-(1),] # remove this event from the Queue
# update the status of that host

# call the function Event. If the function event returns 1, then save it in the history vector
status_id <- Status[Current[,3]]
if (eval(event(Current[,1],Current[,2],status_id,Current[,3]))==1){ # If the function event returns a new event save it in the history, and draw it on the   map
  History <-
    rbind(History,data.frame(Event_time=Current[,1],Type_event=Current[,2],host_id=Current[,3],x_coord=Re(Coord[as.numeric(Current[,3])]),y_coord=Im(Coord[as.numeric(Current[,3])]))) # record it in the history vector. In the history vector add the  coord
}
# to decide which is the next event compare the time of the 2 vectors (infected,removed)
next_events <- rbind(List_to_infect[1,],List_to_remove[1,])
index_next_event <- which.min(next_events[,1])
Queue <- rbind(Queue,next_events[index_next_event,])
if(all.equal(cbind(List_to_infect[1,2],List_to_infect[1,3]),cbind(Queue[1,2],Queue[1,3]))==TRUE){
  List_to_infect <- List_to_infect[-c(1),]
} else
  if(all.equal(cbind(List_to_remove[1,2],List_to_remove[1,3]),cbind(Queue[1,2],Queue[1,3]))==TRUE){
    List_to_remove <- List_to_remove[-c(1),]
  }
#####################################################################