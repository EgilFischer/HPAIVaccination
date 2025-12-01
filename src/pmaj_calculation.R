#########################################################
#                                                        
#                  Probability of a major outbreak                               
#                  Description                           
#                                                        
#                  Author:  E.A.J. Fischer                             
#                  Contact:   e.a.j.fischer@uu.nl                           
#                  Creation date: 4-4-2023                         
#########################################################
# Uses code from Probability of a major outbreak (probMajorOutbreak) but simplified to use any k or p.protect input

#source("./src/loadLibraries.R") 

 itypes = 2;
 param.list.baseline.layer <- list(
   N0 = 46000, # Population size
   itypes = itypes, # Number of types
   beta = matrix(c(1.13, 1.13,0.05,0.05),ncol = itypes),#,#transmission coefficient matrix for a 2x2 matrix (1 -> 1, 1->2, 2-> 1, 2-> 2)#Use values for infectivity and infectious periods from Sitaris et al 2016 https://doi.org/10.1098/rsif.2015.0976#Type 1  = not protected by vaccination and type 2 = protected by vaccination
   infectious.period = c(3.0, 4.0), # Duration infectious period
   #variance.infectious.period = c(0.45, 0.80),
   k.infectious = c(2,2), # k parameter Erlang distribution
   pdie = c(1.0, 0.01), # Probability of dying at end of infectious period
   mortRate = 0 #0.0005 # Per capita death rate
 )
# 
# variance.infectious.period <- param.list.baseline.layer$infectious.period^2 / param.list.baseline.layer$k.infectious
# 
 param.list <- param.list.baseline.layer

model.gamma <- function(u, param.list, mortality = FALSE) {
  with(param.list, {
    # Shape parameter for infectious period
    if(exists("k.infectious")){if(length(k.infectious)==1){a<- c(k.infectious, k.infectious)}else{a<-k.infectious}}else{
    a <- infectious.period^2 / variance.infectious.period}
    
    # Scale parameters to mean infectious period
    R <- beta * infectious.period * 
      (1 * as.numeric(!mortality) + as.numeric(mortality) * (1 + (infectious.period / a) * mortRate)^(-1 - a))
    
    # determine fraction in each population
    p <- c(N[1] / sum(N), N[2] / sum(N))
    
    # calculate the probability generating functions
    return(c(
      f1 = -u[1] + (1 + ((1 - u[1]) * p[1] * R[1, 1] + (1 - u[2]) * p[2] * R[2, 1]) / a[1])^-a[1],
      f2 = -u[2] + (1 + ((1 - u[1]) * p[1] * R[1, 2] + (1 - u[2]) * p[2] * R[2, 2]) / a[2])^-a[2]
    ))
  })
}


Rmodel.gamma <- function(param.list, p, mortality = TRUE) {
  with(param.list, {
    N <- c((1 - p) * N0, p * N0)
    # Shape parameter for infectious period
    if(exists("variance.infectious.period")){
    a <- infectious.period^2 / variance.infectious.period}else
      {a<- ifelse(length(k.infectious)==1,c(k.infectious, k.infectious) ,k.infectious)}
    
    # Scale parameters to mean infectious period
    infectious.period.disc <- infectious.period *
      (1 * as.numeric(!mortality) + as.numeric(mortality) * (1 + (infectious.period / a) * mortRate)^(-1 - a))
    
    # Calculate Rv
    return(0.5 * (beta[1, 1] * (N[1] / sum(N)) * infectious.period.disc[1] +
                    beta[2, 2] * (N[2] / sum(N)) * infectious.period.disc[2] +
                    sqrt((beta[2, 2] * (N[2] / sum(N)) * infectious.period.disc[2] -
                            beta[1, 1] * (N[1] / sum(N)) * infectious.period.disc[1])^2 +
                           4 * (beta[2, 1] * (N[2] / sum(N)) * infectious.period.disc[1] *
                                  beta[1,2] * (N[1] / sum(N))) * infectious.period.disc[2])))
  })
}

#critical vaccination coverage
pc <- function(param.list, mortality = TRUE) {
  with(param.list, {
    #throw warning if calculation is incorrect
    if(beta[1,1]!= beta[2,1] | beta[2,2]!= beta[1,2]) {warning("The calculation in this function only holds when the infectiousness of types differs and not if also susceptibility is reduced.")}
    # Calculate critical vaccination coverage when only infectiousness is reduced
    return((1 - beta[1,1] * infectious.period[1])/(beta[2,2] * infectious.period[2] - beta[1,1] * infectious.period[1]))
  })
}

#calculate the probability of a major outbreak
model <- model.gamma

pmajor <- function(param.list,p,n, prefer.low.titre = TRUE,average_by_p = FALSE,...){
  param.list$N <- c((1 - p) * param.list$N0, p * param.list$N0)
  
  dat<- data.frame(p = p)
  dat <- cbind(dat,data.frame(t(multiroot(f = model, start = c(.5,.5),
                                           rtol = 1E-19,
                                           parms = c(param.list,N = list(param.list$N0*c(1-p,p))))$root)))
  dat$R0 <- Rmodel.gamma(param.list = param.list,p = 0);
  dat$R_hightitre <- Rmodel.gamma(param.list = param.list,p = 1);
  dat$Rv <- Rmodel.gamma(param.list = param.list,p = p);
  dat$pc <- pc(param.list = param.list);
  names(dat)<- c("p","q1","q2","R0","R_high","Rv","pc");
  
  if(length(n) ==2){n1 = n[1];n2=n[2]}else{
  if(prefer.low.titre ==TRUE){
    n1 <- ceiling(n*(1-p)); 
    n2 <- floor(n*p);}else{
      n1 <- floor(n*(1-p)); 
      n2 <- ceiling(n*p);
      }};
  #if n == 1 it is possible to average by taking the probability that either one or the other based on the proportion low titre is the first to get infected
  dat$pmajor <- if(sum(n)==1 & average_by_p){max((1-p)*(1-dat$q1) + p*(1-dat$q2),0)}else{
  max(1-dat$q1^n1 * dat$q2^n2,0)}
  return(dat)
  
}

#Rmodel.gamma(param.list , 0)

# 
# 
# # Simulate results for p = 0 to 1, given introduction of n infectious birds
#  n <- 1
#  p_values <- seq(0, 1, by = 0.050)
# #p_values <- 0.50
#  p_results <- list()
#  for (p in p_values) {
#    p_results[[as.character(p)]] <- pmajor(param.list, p = p, n = n)
#  }
#  results_combined <- do.call(rbind, p_results)
#  print(results_combined)
# 
#  param.list$N <-c(0.5*param.list$N0 ,0.5*param.list$N0 )
# 
# model.gamma(c(.1,.1), param.list, TRUE)
# 
# param.list <- param.list.baseline.layer
# param.list
# p <- 0.5
# param.list$N <-c(p*param.list$N0 ,(1-p)*param.list$N0 )
# param.list$beta <- matrix(c(5,5,0.5,0.5), ncol =2)
# param.list$infectious.period <- c(1,1)
# param.list$k.infectious <- c(2,2)
# 
# fu1<- function(u1){ifelse(param.list$N[1] == 0, 1, uniroot(f = function(u2){model.gamma(c(u1,u2), param.list, TRUE)[1]},c(0,2))$root)}
# fu2 <- function(u1){uniroot(f = function(u2){model.gamma(c(u1,u2), param.list, TRUE)[2]},c(0,1))$root}
# model.gamma(c(.2,0.1), param.list, TRUE)
# seq_u1 <- seq(0.0,1.0,.01)
# rm(u1u2)
# u1u2 <- data.frame(u1 = seq_u1,
#            u1_iso = sapply(seq_u1,function(u){tryCatch(fu1(u), error = function(e){NA})}),
#            u2_iso = sapply(seq_u1,fu2))
# ggplot(data = u1u2)+
#   geom_path(aes(u1,u1_iso))+
#   geom_path(aes(u1,u2_iso))+ylab("u2")
