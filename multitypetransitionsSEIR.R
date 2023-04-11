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
SEIRtransmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dS1 <- S1 * (-(beta11*I1+beta21*I2)/N - mu  - T12) + T21*S2
    dS2 <- S2 * (-(beta12*I1+beta22*I2)/N - mu  - T21) + T12*S1
    dE1 <- S1 * (beta11*I1+beta21*I2)/N - mu*E1 - sigma*E1 
    dE2 <- S2 * (beta12*I1+beta22*I2)/N - mu*E2 - sigma*E2 
    dI1 <- sigma * E1 - (gamma + mu) * I1
    dI2 <- sigma * E2 - (gamma + mu) * I2
    dR1 <- gamma * I1 - mu * R1
    dR2 <- gamma * I2 - mu * R2
    dN <- mu * N
    return(list(c(dS1,dS2,dE1,dE2,dI1,dI2,dR1,dR2,dN)))
  })
}