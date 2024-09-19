#########################################################
#                                                        
#                  Probability of a major outbreak                               
#                  Description                           
#                                                        
#                  Author:  E.A.J. Fischer                             
#                  Contact:   e.a.j.fischer@uu.nl                           
#                  Creation date: 4-4-2023                         
#########################################################

#load libraries
source("./src/loadLibraries.R") 


#Function to find the root of q1,q2 = pgf (probability generating functions) for q1, q2 for offspring of type 1 or  2, given ultimate extinction e.g. I1(\[Infinity]) = 0;  I2(\[Infinity]) =0
#In this function f1 = -u1 + pgf1 and f2 = -u2 + pgf2 
#exponential distributed infectious period
model.exp <- function(u, param.list) {with((param.list),
  return(c(f1 =-u[1] + (gamma[1] + mu + beta[1,1]* (N[1]/(N[1] +N[2]) )* u[1]^2 + beta[2,1]* (N[2] /(N[1] + N[2]))* u[1]* u[2])/(gamma[1]+ mu + beta[1,1]* (N[1]/(N[1] + N[2]))+ beta[2,1]* N[2] /(N[1] + N[2])),
                       f2 = -u[2] + (gamma[2] + mu + beta[1,2]* (N[1]/(N[1] +N[2]) )* u[1] * u[2] + beta[2,2]* (N[2] /(N[1] + N[2]))* u[2]^2)/(gamma[2]+ mu + beta[1,2]* (N[1]/(N[1] + N[2]))+ beta[2,2]* N[2] /(N[1] + N[2])))))}

#gamma distributed infectious period 
#From Clancy and Pearce 2013
model.gamma <- function(u, param.list, mortality = TRUE) {
  with((param.list),{
  #shape parameter for infectious period
  a = infectious.period^2 / variance.infectious.period;
    
  #scale parameters to mean infectious period (if mortality is true, discounted for background mortality)
  R = beta*infectious.period * (1 * as.numeric(!mortality) + as.numeric(mortality)*(1 + (infectious.period/a)*mortRate)^(-1-a) );
 
  #determine fraction in each population
  p = c(N[1]/sum(N), N[2]/sum(N));

  return(c(f1 =-u[1] + (1 + ((1 - u[1]) * p[1] * R[1, 1] + (1 - u[2]) * p[2] * R[2, 1])/a[1])^-a[1], #type 1
          f2 = -u[2] + (1 + ((1 - u[1]) * p[1] * R[1, 2] + (1 - u[2]) * p[2] * R[2, 2])/a[2])^-a[2])) #type 2
       }
    )}




#select the model to use
model <- model.exp
  

#
Rmodel.exp <- function(param.list, p){with(param.list,
  { N <- c(1-p,p)*N0;
    return(0.5*(beta[1,1]*(N[1]/(N[1]+N[2]))/(gamma[1]+mu) + 
                  beta[2,2]*(N[2]/(N[1]+N[2])) /(gamma[2]+mu) + 
                  sqrt((beta[2,2]*(N[2]/(N[1]+N[2]))/(gamma[2]+mu) - beta[1,1]*(N[1]/(N[1]+N[2]))/(gamma[1]+mu))^2 + 
                       4 * (beta[1,2]*(N[2]/(N[1]+N[2]))*beta[2,1]*(N[1]/(N[1]+N[2])))/((gamma[2]+mu)*(gamma[1]+mu)))))
  })}

Rmodel.gamma <- function(param.list, p, mortality = TRUE){with(param.list,
                                       { N <- c(1-p,p)*N0;
                                       #shape parameter for infectious period
                                       a = infectious.period^2 / variance.infectious.period;
                                       
                                       #scale parameters to mean infectious period (if mortality is true, discounted for background mortality)
                                       infectious.period.disc <- infectious.period * (1 * as.numeric(!mortality) + as.numeric(mortality)*(1 + (infectious.period/a)*mortRate)^(-1-a) );
                                       return(0.5*(beta[1,1]*(N[1]/(N[1]+N[2])) * infectious.period.disc[1] + 
                                                     beta[2,2]*(N[2]/(N[1]+N[2]))* infectious.period.disc[2] + 
                                                     sqrt((beta[2,2]*(N[2]/(N[1]+N[2]))* infectious.period.disc[2] - beta[1,1]*(N[1]/(N[1]+N[2]))* infectious.period.disc[1])^2 + 
                                                            4 * (beta[1,2]*(N[2]/(N[1]+N[2]))*infectious.period.disc[1]*beta[2,1]*(N[1]/(N[1]+N[2])))*infectious.period.disc[2])))
                                       })}

q1q2 <- function(param.list,...){
  dat<- data.frame(t(sapply(FUN = function(p){c(p,multiroot(f = model, start = c(.01,.01),
                                                      rtol = 1E-19,
                                                      parms = c(param.list,N = list(param.list$N0*c(1-p,p))),...)$root)}, X = seq(0,1,.01))))
  dat$Rv <- sapply(FUN = Rmodel, param.list = param.list,X = seq(0,1,.01),... )
  names(dat)<- c("p","q1","q2","Rv");
  return(dat)
}

pmajor <- function(param.list,p,n,...){
  dat<- data.frame(p = p)
  dat <- cbind(dat,data.frame(t( multiroot(f = model, start = c(.1,.1),
                                                            rtol = 1E-19,
                                                            parms = c(param.list,N = list(param.list$N0*c(1-p,p))))$root)))
  dat$Rv <- Rmodel(param.list = param.list,p = p);
  names(dat)<- c("p","q1","q2","Rv");
  dat$pmajor <- max(1-dat$q1^round(n*(1-p)) * dat$q2^round(n*p),0);
  return(dat)
}


#test the exponential and gamma are equal if the variance is the square of the mean (e.g. shape has value 1) ####
param.list <- list(beta = matrix(c(1.13,1.13,0.05,0.05),nrow =2),  
                   infectious.period = c(3,4), 
                   mortRate = c(7.14E-05),
                   N0 = 64000)
param.list$gamma <- 1/param.list$infectious.period
param.list$variance.infectious.period <- param.list$infectious.period^2
param.list$mu <- param.list$mortRate
#exponential distributed infectious period
model <- model.exp
Rmodel <- Rmodel.exp
q1q2.exp <- q1q2(param.list)
#exponential distributed infectious period using gamma distribution with shape = 1

model <- model.gamma
Rmodel <- Rmodel.gamma
q1q2.gamma <-q1q2(param.list)
#gamma distributed infectious period using gamma distribution with shape = 20
param.list$variance.infectious.period <- (param.list$infectious.period^2) /20
q1q2.gamma.var <-q1q2(param.list, mortality = TRUE)


#combine 
q1q2.df <- rbind(cbind(q1q2.exp,data.frame(infectiousperiod = "exponential")),
                    cbind(q1q2.gamma,data.frame(infectiousperiod = "gamma.exp")),
              cbind(q1q2.gamma.var,data.frame(infectiousperiod = "gamma.var")))

 
#initial number of infected animals per group
 i1 = .5;
 i2 = .5;
 
 ggplot(data =q1q2.df)+
   geom_path(aes(p, q1^(i1+i2), colour = "1",linetype = infectiousperiod), linewidth = 1.5)+
   geom_path(aes(p, q2^(i1+i2), colour = "2",linetype = infectiousperiod), linewidth = 1.5)+
   geom_path(aes(p,(q1^i1* q2^i2), colour = "1&2",linetype = infectiousperiod), linewidth = 1.5)+
   geom_path(aes(p,Rv, colour = "R",linetype = infectiousperiod), linewidth = 1.5)+
   xlab("Vaccinated coverage")+
   ylab("Probability of a minor outbreak")+ 
     scale_colour_manual(name = paste("Introduction by",i1+i2,"birds"),
                         labels = c("Type 1 only","Type 1 & 2","Type 2 only","R"),
                         values = c("red","blue","green","black"))+
  facet_grid(infectiousperiod~.)

   



