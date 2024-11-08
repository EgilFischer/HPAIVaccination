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
load.sims<- function(path){
  sim.files <- list.files(path, pattern = ".RDS");if(length(sim.files) == 0) stop("No simulations in this folder");
  output <- lapply(sim.files,function(i){readRDS(paste0(path,i))$out})
  pars<- lapply(sim.files,function(i){readRDS(paste0(path,i))$pars})
  
  return(list(output = output,pars =pars))
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
plot.output.grid <- function(output,vars,
                             title = NULL, 
                             frac = NULL , 
                             scales = "fixed", 
                             scenario.label = NULL, 
                             scenario.levels = NULL,
                             itype.label = NULL, 
                             itype.levels =NULL,
                             legend.position = "none"){
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
    theme(legend.position = legend.position )+
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
  tot.exposure<- data.frame(output)%>% 
    group_by(run)%>%reframe(total.exposure = sum(beta.human[1,1]*I.1*dt+beta.human[1,2]* I.2*dt))
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

#fit distributions
fit.distribution <-function(data, 
                            distr = "gamma", 
                            method = c("pas.det.time","ac.det.time","min.det.time"),
                            scale  = 1, #prevent numerical problems with very large numbers
                            long.form = FALSE, 
                            long.var = "detection.exposure",
                            added.vars =NULL){
  dists <- NULL;
  for(its in unique(data$scenario)){
    det.data <- data%>%filter(scenario == its)
    for(m in method)
     { 

      # #select unique values for passive detection
       if(m == "pas.det.time"){
          #number of simulations
          n <- max(det.data$run)
          if(!long.form){ tmp <-  det.data%>%filter(det.data[,m]!=Inf)%>%filter(rep ==1) %>%dplyr::select(all_of(m)) }
          else {
            tmp <-  det.data%>%filter(det.data[, long.var]!=Inf)%>%filter(rep ==1) %>%dplyr::filter(detection.method==m)%>%select(all_of(long.var))
          }
        }else{
         #number of simulations
         n <- max(det.data$run)*max(det.data$rep)
         if(!long.form){
           tmp <- det.data%>%dplyr::select(all_of(m))%>%filter(det.data[,m]!=Inf)
         }
          else
          {
            tmp <-  det.data%>%filter(det.data[, long.var]!=Inf)%>%dplyr::filter(detection.method==m)%>%select(all_of(long.var))
          }
        }

      if(length(tmp[,1])>1)
      {
        fit <- fitdistrplus::fitdist(data = tmp[,1]/scale, distr = distr, method = "mle")
        summary.fit <-summary(fit)
        tmp.dist<- data.frame(shape = as.numeric(summary.fit$estimate[1]),
                              rate = as.numeric(summary.fit$estimate[2])/scale,
                              scenario = its,
                              detection.method = m,
                              pdetect = length(tmp[,1])/n)
        for(v in c(1:length(add.vars)))
          {
            tmp.dist[,names(added.vars)[v]]<- added.vars[[v]][its]
          }
      }else{
        tmp.dist <- data.frame(shape =c(0), rate =c(10^-5),scenario = its, detection.method = m, pdetect = 0)
        for(v in c(1:length(add.vars)))
        {
          tmp.dist[,names(added.vars)[v]]<- added.vars[[v]][its]
        }
        }
      dists <-rbind(dists,tmp.dist)
       }
    }
  
  dists$mean <- dists$shape/dists$rate
  dists$var <- dists$shape/dists$rate^2
  return(dists)
  
}
