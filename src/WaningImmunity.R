#load libraries
source("./src/loadLibraries.R") 

TiterWaning <- read_csv("TiterWaning.csv")
TiterWaning <- reshape2::melt(TiterWaning,id.vars = c("Time"))
mod <- lm(value ~ Time * variable, data = TiterWaning)
drop1(mod)
mod.uni <- lm(value ~ Time, data = TiterWaning)

TiterWaning<- cbind(TiterWaning, pred = predict(mod), preduni = predict(mod.uni))
ggplot(data = TiterWaning)+ 
  geom_point(aes(x = Time,
                 y = value,
                 colour = variable))+
  geom_path(aes(x = Time,
                 y = pred,
                 colour = variable))+
  geom_path(aes(x = Time,
                y = preduni),
                colour = "black")
#rate of decline
runvac<- mod$coefficients["Time"]
rvac <- mod$coefficients["Time"]+mod$coefficients["Time:variableVaccinated"]
runi <- mod.uni$coefficients["Time"]
#time until dropping below a certain titer
starttiter <- 8
lowtiter <- 2
av.time.runvac <- (lowtiter-starttiter)/runvac
av.time.rvac <- (lowtiter-starttiter)/rvac
av.time.runi <- (lowtiter-starttiter)/runi
1/av.time.runvac
1/av.time.rvac
1/av.time.runi
