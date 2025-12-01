##### Test for the simulation model ########
source("./src/faster_multitypetransitionsSEIR_tleap.R")
source("./src/default_parameter_list.R")
source("./src/pmaj_calculation.R")

#Test 1: R0 = 0####
#Expectation no transmission #

#Changes to parameters:
test.param.list <- param_list_default
test.param.list$beta<- 0*test.param.list$beta
test.param.list$runs <- 10
test.param.list$intro.time<- 1
test.param.list$length.round <- 5

if(exists("x")){rm(x)};x <-sim.multitypeSEIR_tleap(test.param.list,
                                                   seed = NULL)

#expected result maximum number of infections equals the test.param.list$no.introduction value
test.param.list$no.introduction
data.frame(x)%>%reframe(.by = run, 
                        max.inf = max(L.1+L.2+I.1+ I.2))

#Results: Correct: no more than 1 infected





# FIXED!  
#Test 2: Waning shape and timing
#Expectation: Fraction susceptible should follow the distribution of waning

#Changes to parameters:
test.param.list <- param_list_default
test.param.list$runs <- 10
test.param.list$beta<- 0*test.param.list$beta #check waning with R0 = 0 so no interference of the infection#

if(exists("x")){rm(x)};x <-sim.multitypeSEIR_tleap(test.param.list,
                                                   seed = NULL)

plot(x$time, 1-x$S[,1]/rowSums(x$S))

#plot expectation
test.param.list$age_at_d0
test.param.list$trans.mean.wane # = scale * shape
test.param.list$trans.var.wane # = scale^2 * shape
scale = test.param.list$trans.var.wane /test.param.list$trans.mean.wane #var/mean
shape = test.param.list$trans.mean.wane /scale #mean/scale

lines(seq(0,560,1),1-sapply(seq(0,560,1)+test.param.list$age_at_d0,function(t){pgamma(t,shape = shape,
       scale =scale)}), col = "red" )
#Result: Incorrect: timing of waning is shifted by age_at_d0



#FIXED! 
#Test 3: Probability of a major outbreak when no immune waning and  several levels of protection ###
#See also: pmaj_calculation.R for more in-depth pmaj calculations
#Expectation: Higher protection (p.protect) should decrease the probability of a major outbreak

#Changes to parameters:
test.param.list <- param_list_default
test.param.list$k.infectious<- c(2,2)
test.param.list$trans.mean.wane = Inf
test.param.list$trans.mean.buildup = Inf
test.param.list$no.introduction <- c(1,0)
test.param.list$N0 <- 1000
test.param.list$length.round <-50
test.param.list$deltat <- 0.005
test.param.list$runs <- 50
#expected fraction minor outbreaks ####
test.param.list$variance.infectious.period <- (test.param.list$infectious.period^2) /test.param.list$k.infectious
res.test <- NULL
for(i in seq(0,1,0.05)){
  res.test<- rbind(res.test,pmajor(test.param.list,
                                   p = i,test.param.list$no.introduction))
}
res.test

res.test_high <- NULL
test.param.list$k.infectious<- c(2,2)
for(i in seq(0,1,0.05)){
  res.test_high<- rbind(res.test_high,pmajor(test.param.list,
                                   p = i,c(0,1)))
}
res.test_high

#test p.protect 0, 0.5,0.7, 1.0
if(exists("x")){rm(x)};x <- NULL
test.param.list$k.latency <- 2
test.param.list$k.infectious <- 2
for(p in c(0,0.1,0.2,0.3,0.40,0.50,0.6,0.7,0.75,0.8,0.85,0.9,0.95,1.0)){
  
 test.param.list$p.protect <- p; 
 x <- rbind(x,cbind(data.frame(sim.multitypeSEIR_tleap(test.param.list,
                                                     seed = NULL)), p.protect = p))
}
if(exists("results")){rm(results)};results <- NULL
results <- x %>%
  reframe(.by = c(run, p.protect), 
          max.inf = max(DL.1 + DL.2 + DI.1 + DI.2 + DR.1 + DR.2 +R.1+R.2)  > (0.05 * test.param.list$N0)) %>%
  reframe(.by = p.protect,
          pmajor = mean(max.inf))
results <- results[order(results$p.protect),]  
p_protect_values <- results$p.protect
p_major_results <- results$pmajor
plot(
  p_protect_values, p_major_results,
  type = "o" ,
  ylim = c(0,1)
  
)
lines(x = res.test$p, y = res.test$pmajor, col = "red")
lines(x = res.test_high$p, y = res.test_high$pmajor, col = "red", lt = "dashed")

#Correct pmajor as expected

hist_dat <- x %>%filter(p.protect == 1.)%>%
  reframe(.by = c(run, p.protect), 
          max.inf = max(DL.1 + DL.2 + DI.1 + DI.2 + DR.1 + DR.2 +L.1 + L.2 + I.1+I.2+R.1+R.2)-1 )
hist(hist_dat$max.inf/1000,breaks = 10)


#Test 4: R0 < 1 ###
#Expectation: no major outbreaks

#Changes to parameters
test.param.list <- param_list_default
test.param.list$runs <- 10
test.param.list$beta[,1]<- 0.9/test.param.list$infectious.period
test.param.list$beta[,2]<- 0


if(exists("x")){rm(x)};x <-sim.multitypeSEIR_tleap(test.param.list,
                                                   seed = NULL)

#expected result maximum number of infections equals the test.param.list$no.introduction value
data.frame(x)%>%reframe(.by = run, 
                        max.inf = max(L.1+L.2+I.1+ I.2))
#Output:
#  run    max.inf
#    1       1
#    2       1
#    3       1
#    4       1
#    5       1
#    6       1
#    7       1
#    8       1
#    9       1
#   10       1
#Result: Correct 





#Test 5: When beta is 0.05 for all birds and zero protection:
#Expectation: no major outbreak (pmaj= 1-1/0.2) 

#Changes to parameters:
test.param.list <- param_list_default
test.param.list$runs <- 100
test.param.list$beta <- matrix(0.05, ncol = test.param.list$itypes, nrow = test.param.list$itypes)
test.param.list$p.protect <- 0.7

if(exists("x")){rm(x)};x <-sim.multitypeSEIR_tleap(test.param.list, seed = NULL)

data.frame(x) %>%
  reframe(.by = run, max.inf = max(L.1 + L.2 + I.1 + I.2))

#Output:
# run   max.inf
#    1      20
#    2      19
#    3       8
#    4       2
#    5       4
#    6       1
#    7       1
#    8       1
#    9       1
#   10       1
#Results: Correct - no major outbreaks





#Test 6: When beta is 1.13 for all birds and zero protection:
#Expectation: major outbreak ~70% (pmaj= 1-1/3.39)

#Changes to parameters:
test.param.list <- param_list_default
test.param.list$runs <- 10
test.param.list$beta <- matrix(c(1.13, 1.13,1.13,1.13), ncol = test.param.list$itypes, nrow = test.param.list$itypes)
test.param.list$p.protect <- 1

if(exists("x")){rm(x)};x <-sim.multitypeSEIR_tleap(test.param.list, seed = NULL)

results <- data.frame(x)
max_inf <- tapply(rowSums(results[, c("L.1", "L.2", "I.1", "I.2")]), results$run, max)
threshold <- 0.1 * test.param.list$N0
major_outbreak <- max_inf > threshold
p_major <- mean(major_outbreak)
print(p_major)

#Output:
#  run  max.inf
#    1   45797
#    2   45804
#    3   45789
#    4   45787
#    5   45789
#    6   45798
#    7   45810
#    8   45783
#    9   45795
#   10   45781
#Results: Incorrect, too many major outbreaks? Even with 100 runs, p.maj is 1 





#Test X: Flock production 
#Expectation: when age_at_d0 is 0, no eggs until ~120

#Changes to parameters:
age_at_d0 <- 0
flock_production <- function(t, age_at_d0) {
  t_shift <- (t + age_at_d0)/7 #function assumes week 0 = birth, model uses day 0 = first day of maturity
  a <- 126.65
  b <- 0.01251
  c <- 1.164
  d <- 23.82
  ((a * exp(-b * t_shift)) / (1 + exp(-c * (t_shift - d))))/100
}
sim.time <- seq(0,550, by = 1)
production_values <- sapply(sim.time, flock_production, age_at_d0 = age_at_d0)
plot(sim.time, production_values, type = "l")

#Results: Correct 





