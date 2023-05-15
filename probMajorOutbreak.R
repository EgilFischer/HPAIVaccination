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
source("loadLibraries.R") 


#probability generating functions q1, q2 for offspring of type 1 or  2, given ultimate extinction e.g. I1(\[Infinity]) = 0;  I2(\[Infinity]) =0
model <- function(u, param.list) {with((param.list),
  return(c(f1 =-u[1] + (gamma[1] + mu + beta[1,1]* (N[1]/(N[1] +N[2]) )* u[1]^2 + beta[2,1]* (N[2] /(N[1] + N[2]))* u[1]* u[2])/(gamma[1]+ mu + beta[1,1]* (N[1]/(N[1] + N[2]))+ beta[2,1]* N[2] /(N[1] + N[2])),
                       f2 = -u[2] + (gamma[2] + mu + beta[1,2]* (N[1]/(N[1] +N[2]) )* u[1] * u[2] + beta[2,2]* (N[2] /(N[1] + N[2]))* u[2]^2)/(gamma[2]+ mu + beta[1,2]* (N[1]/(N[1] + N[2]))+ beta[2,2]* N[2] /(N[1] + N[2])))))}

Rmodel <- function(param.list, p){with(param.list,
  { N <- c(1-p,p)*N0;
    return(0.5*(beta[1,1]*(N[1]/(N[1]+N[2]))/(gamma[1]+mu) + 
                  beta[2,2]*(N[2]/(N[1]+N[2])) /(gamma[2]+mu) + 
                  sqrt((beta[2,2]*(N[2]/(N[1]+N[2]))/(gamma[2]+mu) - beta[1,1]*(N[1]/(N[1]+N[2]))/(gamma[1]+mu))^2 + 
                       4 * (beta[1,2]*(N[2]/(N[1]+N[2]))*beta[2,1]*(N[1]/(N[1]+N[2])))/((gamma[2]+mu)*(gamma[1]+mu)))))
  })}

q1q2 <- function(param.list){
  dat<- data.frame(t(sapply(FUN = function(p){c(p,multiroot(f = model, start = c(.01,.1),
                                                      rtol = 1E-15,
                                                      parms = c(param.list,N = list(param.list$N0*c(1-p,p))))$root)}, X = seq(0,1,.01))))
  dat$Rv <- sapply(FUN = Rmodel, param.list = param.list,X = seq(0,1,.01) )
  names(dat)<- c("p","q1","q2","Rv");
  return(dat)
}

# q1q2 <- data.frame(t(sapply(FUN = function(p){c(p,multiroot(f = model, start = c(.01,
#                                                                                 0.1),
#                                                       rtol = 1E-15,parms = list(N = c(10000 * (1-p),10000 * p),
#                                                                                gamma = c(1/1.47,1/1.47), beta = matrix(c(3.84,3.84,0.058,0.058),nrow = 2), mu = 0.0005))$root)}, X = seq(0,1,.05))))
# q1q2
# names(q1q2)<- c("p","q1","q2")
# 
# #initial number of infected animals per group
# i1 = 5;
# i2 = 5;
# 
# ggplot(data =q1q2)+
#   geom_path(aes(p,1 - q1^(i1+i2), colour = "1"))+
#   geom_path(aes(p,1-  q2^(i1+i2), colour = "2"))+
#   geom_path(aes(p,1-(q1^i1* q2^i2), colour = "1&2"))+
#   xlab("Vaccinated coverage")+
#   ylab("Probability of a major outbreak")+ 
#   scale_colour_manual(name = paste("Introduction by",i1+i2,"birds"),
#                       labels = c("Type 1 only","Type 2 only","Type 1 & 2"),
#                       values = c("red","blue","black"))
