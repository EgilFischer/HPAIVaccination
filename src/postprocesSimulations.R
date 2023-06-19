#########################################################
#                                                        
#                 Post-proces simulations                                
#                                                        
#                  Author: Egil Fischer                               
#                  Contact: e.a.j.fischer@uu.nl                             
#                  Creation date: 28-3-2023                         
#########################################################

source("./src/loadLibraries.R") 

#load simulations
load.sims <- function(path, interval = NULL, params = TRUE){
  sims <- list.files(path, pattern = ".RDS");if(length(sims) == 0) stop("No simulations in this folder");
  #recode runs
  run.counter =1;
   for(i in sims){
     if(is.null(interval)){
     if(!exists("output"))
        {
          tmp.output <-as.data.frame(readRDS(paste0(path,"/",i))$out)
          tmp.output$run <- run.counter
          run.counter = run.counter+1}else
        {
          tmp <- as.data.frame(readRDS(paste0(path,"/",i))$out)
          tmp$run <- run.counter
          run.counter = run.counter+1
          tmp.output <- rbind(tmp.output,tmp)}
        }else #merge into intervals
        {
          tmp <- as.data.frame(readRDS(paste0(path,"/",i))$out);
          tmp$interval.index <- tmp$time%/%interval;
          tmp$tround <- interval*(tmp$interval.index+sign(tmp$time))
          tmp <- tmp%>%group_by(tround)%>%slice(n())%>%ungroup
            if(!exists("tmp.output")) {
              tmp.output <-tmp
              }else{ 
                tmp.output <- rbind(tmp.output,tmp)
            }
          
        }
     if(params){if(!exists("param_lists")){param_lists <- NULL;param_lists[[1]] <- readRDS(paste0(path,"/",i))$pars}else{param_lists[[length(param_lists)+1]]<- readRDS(paste0(path,"/",i))$pars}}
     #recode run
     tmp.output
   }
  return(list(output = tmp.output,pars = if(params){param_lists}else{NULL}))
}


# ##plot the output
plot.output <- function(output,vars,title = NULL ){
  ggplot(data =
           data.frame(output)%>%select(all_of(c("time", "run",vars)))%>%reshape2::melt(id.vars = c("time","run"),value.name = "prevalence",variable.name=c("itype")))+
        geom_step(aes(x = time, y = prevalence,colour = itype, group =run))+
      ylab("#number of birds")+facet_grid(.~itype, scales = "free")+ggtitle(title)
  
}

plot.output.sparse <- function(output,vars,title = NULL, frac = 0.5){
  out <- data.frame(output)%>%sample_n(round(frac*length(output$time)));
  return(plot.output(out, vars, title))
}

# 
# #human exposure ####
# #Assumption of exposure is that the rate at which humans get exposed to infection is proportionate to the infectivity towards chickens
# a <-0.001; #factor scaling towards human exposure
# beta.human <- a*beta
# 
# 
# output.df<- data.frame(output)
# exposure.human <- output.df%>%group_by(run)%>%reframe(time = time,
#                                                       I.1 = I.1,
#                                                       I.2= I.2,
#                                                       psurvdt = exp(-beta.human[1,1]*I.1*dt-beta.human[1,2]* I.2*dt))
# exposure.human <- exposure.human%>%group_by(run)%>%mutate(
#   pinf = 1 - cumprod(psurvdt))
# ggplot(data = exposure.human) + 
#   geom_path(aes(x = time, y = psurvdt, colour = factor(run)))
# ggplot(data = exposure.human) + 
#   geom_path(aes(x = time, y = pinf, colour = factor(run)))



human.exposure.timeseries.multiple.runs <- function(output, beta.human)
  {
  exposure<- data.frame(output)%>% group_by(run)%>%reframe(time = time,
                                            I.1 = I.1,
                                            I.2= I.2,
                                            cum.exposure = cumsum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  return(exposure)
}


human.exposure.timeseries <- function(output, beta.human){
  exposure<- data.frame(output)%>% reframe(time = time,
                                               I.1 = I.1,
                                               I.2= I.2,
                                               cum.exposure = cumsum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  
}

human.exposure.total.multiple.runs<- function(output, beta.human,detection.time,var = "min.det.time", reps = 1)
{
  if(is.null(detection.time$rep))detection.time$rep <- 1;
  exposure<- data.frame(output)%>% group_by(run)%>%reframe(total.exposure = sum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  exp.det<-c();
  for(i in exposure$run){
    for(j in c(1:reps))
    {
      exp.det<-rbind(exp.det, data.frame(output)%>%filter(run == i & time<= detection.time[detection.time$run == i & detection.time$rep == j, var])%>%reframe(run = i,
                                                                                                                                                          rep = j,
                                                                                                                                                          detection.time = detection.time[detection.time$run == i & detection.time$rep == j, var],
        detection.exposure = sum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt),
        cum.I.1 = sum(I.1*dt),
        cum.I.2 = sum(I.2*dt)))
    }
  }
  exposure<- cbind( data.frame(exp.det),data.frame(total.exposure = rep(exposure$total.exposure, each = reps)))
  return(exposure)
}
# ggplot(data = exposure.human) + geom_path(aes(x = time, y = psurvdt, group = run))
# ggplot(data = exposure.human) + geom_path(aes(x = time, y = pinf, group = run))


