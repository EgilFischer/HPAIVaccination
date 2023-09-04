#########################################################
#                                                        
#                 Post-proces simulations                                
#                                                        
#                  Author: Egil Fischer                               
#                  Contact: e.a.j.fischer@uu.nl                             
#                  Creation date: 28-3-2023                         
#########################################################

source("./src/loadLibraries.R") 

#load simulations ####
load.sims <- function(path, interval = NULL, params = TRUE){
  sims <- list.files(path, pattern = ".RDS");if(length(sims) == 0) stop("No simulations in this folder");
  #recode runs
  run.counter =1;
   for(i in sims){
     if(is.null(interval)){
     if(!exists("tmp.output"))
        {
          tmp.output <-as.data.frame(readRDS(paste0(path,"/",i))$out)
          tmp.output$run <- run.counter
          run.counter = run.counter+1
          }else
        {
          tmp <- as.data.frame(readRDS(paste0(path,"/",i))$out)
          tmp$run <- run.counter
          run.counter = run.counter+1
          tmp.output <- rbind(tmp.output,tmp)}
        }else #merge into intervals and reset dt
        {
          tmp <- as.data.frame(readRDS(paste0(path,"/",i))$out);
          tmp$run <- run.counter;
          run.counter = run.counter+1;
          tmp$interval.index <- tmp$time%/%interval;
          tmp$tround <- interval*(tmp$interval.index+sign(tmp$time))
          tmp <- tmp%>%group_by(tround)%>%dplyr::slice(n())%>%ungroup
          tmp$dt <- c(0,tail(tmp$time,-1)-head(tmp$time,-1) )
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

#Visualization of output ####
#plot the output
plot.output <- function(output,vars,title = NULL, frac = NULL ){
  if(!is.null(frac)) return(plot.output.sparse(output, vars,title,frac))
  ggplot(data =
           data.frame(output)%>%select(all_of(c("time", "run",vars)))%>%reshape2::melt(id.vars = c("time","run"),value.name = "prevalence",variable.name=c("itype")))+
        geom_step(aes(x = time, y = prevalence,colour = itype, group =run))+
      ylab("#number of birds")+facet_grid(.~itype, scales = "free")+ggtitle(title)
  
}

#plot a fraction of the output
plot.output.sparse <- function(output,vars,title = NULL, frac = 0.5){
  out <- data.frame(output)%>%sample_n(round(frac*length(output$time)));
  return(plot.output(out, vars, title))
}

#plot a grid with outputs 
plot.output.grid <- function(output,vars,title = NULL, frac = NULL , scales = "fixed", scenario.label = NULL, scenario.levels = NULL,itype.label = NULL, itype.levels =NULL){
    if(is.null(frac)){out = output}else{
      out <- data.frame(output)%>%sample_n(round(frac*length(output$time)))
    }
  out <- data.frame(out)%>%select(all_of(c(vars,"time", "run","scenario")))%>%reshape2::melt(id.vars = c("time","run","scenario"),value.name = "prevalence",variable.name=c("itype"))
  scenario.levels <- if(is.null(scenario.levels)){unique(out$scenario)}else scenario.levels;
  itype.levels <- if(is.null(itype.levels)){unique(out$itype)}else itypes.levels;
  out$scenario <- factor(out$scenario,levels = scenario.levels)
  out$itype <- factor(out$itype,levels = itype.levels);
  ggplot(data =out)+
    geom_step(aes(x = time, y = prevalence,colour = itype, group =run))+
    ylab("#number of birds")+
    facet_grid(scenario~itype, 
                scales = scales, 
                labeller = labeller(scenario = scenario.label, itype = itype.label))+
    ggtitle(title)
  
}

 
#Calculate human exposure ####
#one time series
human.exposure.timeseries <- function(output, beta.human){
  exposure<- data.frame(output)%>% reframe(time = time,
                                           I.1 = I.1,
                                           I.2= I.2,
                                           cum.exposure = cumsum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  
}

#multiple time series
human.exposure.timeseries.multiple.runs <- function(output, beta.human)
  {
  exposure<- data.frame(output)%>% group_by(run)%>%reframe(time = time,
                                                           cum.I.1 = sum(I.1*dt),
                                                           cum.I.2 = sum(I.2*dt),
                                            cum.exposure = cumsum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  return(exposure)
}


#multiple time series and detection times
human.exposure.detection.multiple.runs<- function(output, beta.human,detection.time,var.det = c("pas.det.time", "ac.det.time", "min.det.time"))
{
  exposure <- NULL;
  for(det.method in var.det)
  {
  if(is.null(detection.time$rep))detection.time$rep <- 1;
  tot.exposure<- data.frame(output)%>% group_by(run)%>%reframe(total.exposure = sum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
  exp.det<-c();
  for(i in tot.exposure$run){
    reps <- detection.time%>%filter(run == i)%>%select("rep")%>%max
    for(j in c(1:reps))
    {
      dettime <-  detection.time[detection.time$run == i & detection.time$rep == j, det.method];
      exp.det<-rbind(exp.det, data.frame(output)%>%filter(run == i & time<=dettime)%>%reframe(run = i,
                                rep = j,
                                detection.time = dettime,
        detection.exposure = sum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt),
        cum.I.1 = sum(I.1*dt),
        cum.I.2 = sum(I.2*dt)))
    }
  }
  exposure<- rbind(exposure,cbind(data.frame(exp.det),
                                  data.frame(total.exposure = rep(c(tot.exposure$total.exposure), each = reps)), 
                                  data.frame(detection.method = c(det.method))))
  }
  return(exposure)
}

