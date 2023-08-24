#load libraries
source("./src/loadLibraries.R") 

#paper Mirzaie
TiterWaningM <- read_csv("./input/TiterWaningMirzaie.csv")
TiterWaningM <- reshape2::melt(TiterWaning,id.vars = c("Time"))

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


#plot waning rates ####
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


#Paper Rudolf et al
TitreWaningR <- read_delim("./input/TitreWaningRudolf2010.csv", delim =";")
library(bbmle)

LL <- function(mu, vaccin.status = "BI"){
  dataToFit <- TitreWaningR%>%filter(VaccinStatus ==vaccin.status)
  -log(sum(dbinom(x = dataToFit$PosH5N2, size = dataToFit$NH5N2, p = f(-mu*dataToFit$TimeWeeks))))
}
curve<-sapply(seq(0,0.05,0.001),LL)
plot(seq(0.,0.05,0.001),curve)
with(TitreWaningR%>%filter(VaccinStatus =="BI"),plot(x = TimeWeeks,y = PosH5N2/NH5N2))
fit<- mle2(LL,c(mu =0.02))
fit
fit@fullcoef/7

#fit generalized logistic model to the data 
genlog <- function(mu,f,nu,t){
  1 - (1+exp(-mu*(t-f)))^(-1/nu)
}
sapply(TitreWaningR$TimeWeeks,genlog, f = 55, mu = 0.1, nu = 0.1 )

LL <- function(mu,f,nu, vaccin.status = "BI"){
  dataToFit <- TitreWaningR%>%filter(VaccinStatus ==vaccin.status)
  -log(sum(dbinom(x = dataToFit$PosH5N2, size = dataToFit$NH5N2, p = sapply(TitreWaningR$TimeWeeks,FUN = genlog, f = f, mu = mu, nu = nu))))
}
sapply(TitreWaningR$TimeWeeks,genlog, f = 1, mu = 1, nu = 1)

fit<- mle2(LL,start = list(f =55, mu = 2),fixed = list(nu =1))
fit
genlogpar <- function(t){with(as.list(fit@fullcoef),
                              1 - (1+exp(-mu*(t-f)))^(-1/nu))
}          
plot((TitreWaningR%>%filter(VaccinStatus =="BI"))$TimeWeeks ,
     sapply((TitreWaningR%>%filter(VaccinStatus =="BI"))$TimeWeeks,genlogpar),"l", ylim=c(0,1))
points((TitreWaningR%>%filter(VaccinStatus =="BI"))$TimeWeeks ,(TitreWaningR%>%filter(VaccinStatus =="BI"))$PosH5N2/60 )

LL.gamma <- function(shape,rate, vaccin.status = "BI"){
  dataToFit <- TitreWaningR%>%filter(VaccinStatus ==vaccin.status)
  with(TitreWaningR%>%filter(VaccinStatus =="BI"),sum((NH5N2-PosH5N2)*log(pgamma(TimeWeeks,shape,rate))+PosH5N2*log((1-pgamma(TimeWeeks,shape,rate)))))
}

LL.gamma(50,.10)
pgamma(77,60,.5)

fit.gamma<- mle2(LL.gamma,start = list(shape = 60,rate =1 ),fixed = list(),method = "Nelder-Mead")
fit.gamma
with((TitreWaningR%>%filter(VaccinStatus =="BI")),
     plot(TimeWeeks,1-pgamma(TimeWeeks,fit.gamma@fullcoef["shape"], fit.gamma@fullcoef["rate"]),"l", ylim = c(0,1)))
with((TitreWaningR%>%filter(VaccinStatus =="BI")),
     points(TimeWeeks,PosH5N2/NH5N2))

