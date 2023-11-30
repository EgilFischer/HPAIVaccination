#load libraries ####
source("./src/loadLibraries.R") 


#Estimate waning rate ####
TiterWaning <- read_csv("./input/TiterWaning.csv")
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


#plot waning rates #
times <- seq(0,200,.1)
wane <- function(t,rate){exp(-rate*t)}
plot.data.waning <- data.frame(time = rep(times,3),
                               Immunity = c(rep(1,length(times)),
                                            wane(times,0.012),
                                            wane(times,0.038)),
                               Waning = c(rep("none",length(times)),
                                          rep("0.012 per day",length(times)),
                                          rep("0.038 per day",length(times))))
ggplot(data = plot.data.waning)+
  geom_path(aes(time,Immunity,colour = Waning), linewidth = 1)+
  geom_vline(xintercept = c(0,1,10,100),linetype = "dotted")+
  geom_hline(yintercept = c(0.75),linetype = "solid")+
  ylim(0, 1.001)
ggsave("./output/figures/Waning.png")


#Get farm sizes
#load spatial input
spatial.input<- read_excel("./input/20230404_AI_AnimalLocations_SizeType_v02.xlsx", sheet = 2)

spatial.input%>%reframe(.by = TYPE,
                        mean = mean(SIZE),
                        median = median(SIZE),
                        Prec01 = quantile(SIZE,0.01),
                        Prec25 = quantile(SIZE,0.25),
                        Prec75 = quantile(SIZE,0.75),
                        Prec99 = quantile(SIZE,0.99))
