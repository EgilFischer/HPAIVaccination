#########################################################
#                                                        
#                  Deterministic multitype transitions model                                 
#                  
#                                                        
#                  Author: Egil Fischer                              
#                  Contact: e.a.j.fischer@uu.nl                             
#                  Creation date                         
#########################################################
#include libraries ####
packages <- c("ggplot2","deSolve","tidyverse")


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

# SEIR transition model ####
SLIRtransmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dS.1 <- S.1 * (-(beta11*I.1+beta21*I.2)/N - mu  - T12) + T21*S.2
    dS.2 <- S.2 * (-(beta12*I.1+beta22*I.2)/N - mu  - T21) + T12*S.1
    dL.1 <- S.1 * (beta11*I.1+beta21*I.2)/N - mu*L.1 - sigma1*L.1 
    dL.2 <- S.2 * (beta12*I.1+beta22*I.2)/N - mu*L.2 - sigma2*L.2 
    dI.1 <- sigma1 * L.1 - (gamma1 + mu) * I.1
    dI.2 <- sigma2 * L.2 - (gamma2 + mu) * I.2
    dR.1 <- gamma1 * I.1 - mu * R.1
    dR.2 <- gamma2 * I.2 - mu * R.2
    dN <- -mu * N
    return(list(c(dS.1,dS.2,dL.1,dL.2,dI.1,dI.2,dR.1,dR.2,dN)))
  })
}
#parameters ####
# params<- c(beta11 = param.list$beta[1,1],
#            beta12 = param.list$beta[2,1],
#            beta21 = param.list$beta[1,2],
#            beta22 = param.list$beta[2,2],
#            mu = param.list$mortRate,
#            T12 = param.list$transRate[1,2],
#            T21 = param.list$transRate[2,1],
#            sigma1 = 1/param.list$latency.period[1],
#            sigma2 = 1/param.list$latency.period[2],
#            gamma1 = 1/param.list$infectious.period[1],
#            gamma2 = 1/param.list$infectious.period[2])

# #initial values ####
# N0 = 75000;
# p.vaccinated = 0.9;
# L0 = 10;
# dt = 0.1;
# max.time  = 42;param.list$max.time;
# init <- c(S.1 = floor((N0-L0)*(1-p.vaccinated)),
#           S.2 = ceiling((N0-L0)*(p.vaccinated)),
#           L.1 = floor(L0*(1-p.vaccinated)),
#           L.2 = ceiling(L0*(1-p.vaccinated)),
#           I.1 =0,
#           I.2 = 0,
#           R.1 = 0,
#           R.2 =0,
#           N = N0)
# 
# #run the model ####
# ode.out <- ode(init,seq(0, max.time, dt), SLIRtransmod, params)
# ggplot(data = as.data.frame(ode.out)) + 
#   geom_path(aes(x = time, y = S.1, colour = "S",linetype = "1"))+
#   geom_path(aes(x = time, y = S.2, colour = "S",linetype = "2")) +
#   geom_path(aes(x = time, y = L.1, colour = "L",linetype = "1"))+
#   geom_path(aes(x = time, y = L.2, colour = "L",linetype = "2")) +
#   geom_path(aes(x = time, y = I.1, colour = "I",linetype = "1"))+
#   geom_path(aes(x = time, y = I.2, colour = "I",linetype = "2")) +
#   geom_path(aes(x = time, y = R.1, colour = "R",linetype = "1"))+
#   geom_path(aes(x = time, y = R.2, colour = "R",linetype = "2")) +
#   ylab("#number of birds")
#            