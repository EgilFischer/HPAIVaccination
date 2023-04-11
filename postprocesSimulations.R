#########################################################
#                                                        
#                 Post-proces simulations                                
#                                                        
#                  Author: Egil Fischer                               
#                  Contact: e.a.j.fischer@uu.nl                             
#                  Creation date: 28-3-2023                         
#########################################################


#include libraries ####
packages <- c("ggplot2","dplyr","tidyverse")


## Now load or install&load all
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

#import simulations
for(i in c(1:10)){
  
  if(i == 1)output <- as.data.frame(readRDS(paste0("20230407outputlayerT1001.RDS"))$out)
  output <- rbind(output,as.data.frame(readRDS(paste0("20230407outputlayerT100",i,".RDS"))$out))
}

##plot the output
ggplot(data = 
         data.frame(output)%>%select(time, run,I.1,I.2)%>%reshape2::melt(id.vars = c("time","run"),value.name = "prevalence",variable.name=c("itype")))+
  #geom_path(aes(x = time, y = S.1, group = run), color = "magenta", linetype = 1)+
  #geom_path(aes(x = time, y = S.2, group = run), color = "magenta", linetype = 2)+
  # geom_step(aes(x = time, y = L.1, group = run), color = "blue", linetype = 1)+
  # geom_step(aes(x = time, y = L.2, group = run), color = "blue", linetype = 2)+
  geom_step(aes(x = time, y = prevalence,colour = itype))+
  #geom_step(aes(x = time, y = I.2, group = run), color = "red", linetype = 2)+
  #geom_step(aes(x = time, y = R.1, group = run), color = "darkgreen")+
  # geom_step(aes(x = time, y = R.2, group = run), color = "darkgreen", linetype = 2)+
  # geom_step(aes(x = time, y = DS.1 + DL.1 +DI.1+DR.1, group = run), color = "black")+
  # geom_step(aes(x = time, y = DS.2 + DL.2 +DI.2+DR.2, group = run), color = "black", linetype = 2)+
  ylab("#number of birds")+facet_grid(run~.)


#human exposure ####
#Assumption of exposure is that the rate at which humans get exposed to infection is proportionate to the infectivity towards chickens
a <-0.001; #factor scaling towards human exposure
beta.human <- a*beta


output.df<- data.frame(output)
exposure.human <- output.df%>%group_by(run)%>%reframe(time = time,
                                                      I.1 = I.1,
                                                      I.2= I.2,
                                                      psurvdt = exp(-beta.human[1,1]*I.1*dt-beta.human[1,2]* I.2*dt))
exposure.human <- exposure.human%>%group_by(run)%>%mutate(
  pinf = 1 - cumprod(psurvdt))
ggplot(data = exposure.human) + 
  geom_path(aes(x = time, y = psurvdt, colour = factor(run)))
ggplot(data = exposure.human) + 
  geom_path(aes(x = time, y = pinf, colour = factor(run)))

#detection time ####
#detect if a certain number of R dead are found within a certain time interval
detection.time.threshold.interval <- function(times,Deaths, time.interval,threshold){
  #determine the number of deaths at the end of the interval
  Dinterval <- data.frame(t(sapply(FUN = function(t){c(time = t,deaths = max(Deaths[times < t]))}, X = seq(1, ceiling(max(times)),time.interval))));
  #define the increment during the interval
  Dinterval$deaths <- c(0, tail(Dinterval$deaths,-1) -  head(Dinterval$deaths,-1))
  #return the first time the threshold is passed /  return Inf if never
  if(max(Dinterval$deaths)<threshold){return(Inf)}
  return(min(Dinterval$time[Dinterval$deaths>=threshold], max(Dinterval$time)));
  
  
}

det.time <- detection.time.threshold.interval(times = output$time[output$run ==1],D = output$DR.1[output$run ==1]+output$DR.2[output$run ==1], time.interval = 1,0.05*20000)
det.time

# #detect if a certain number of R animals are found 
# detection.time.threshold <- function(times,cumD,threshold){
#   return(head(times[cumD >= threshold],1))
# }
# 
# detection.time.threshold(output$time[output$run ==1],output$DR.1[output$run ==1]+output$DR.2[output$run ==1], 0.05*20000)

output.df<- data.frame(output)%>%filter(time <= det.time)
exposure.human <- output.df%>%group_by(run)%>%summarise(time = time,
                                                        I.1 = I.1,
                                                        I.2= I.2,psurvdt = exp(-beta.human[1,1]*I.1*dt-beta.human[1,2]* I.2*dt))
exposure.human <-exposure.human%>%group_by(run) %>% mutate(pinf = 1 - cumprod(psurvdt))
ggplot(data = exposure.human) + geom_path(aes(x = time, y = psurvdt, group = run))
ggplot(data = exposure.human) + geom_path(aes(x = time, y = pinf, group = run))
