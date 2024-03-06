#load libraries
source("./src/loadLibraries.R") 

#paper Mirzaie####
TiterWaningM <- read_csv("./input/TiterWaningMirzaie2020.csv")
TiterWaningM <- reshape2::melt(TiterWaningM,id.vars = c("Time"))

mod <- lm(value ~ Time * variable, data = TiterWaningM)
drop1(mod)
mod.uni <- lm(value ~ Time, data = TiterWaningM)

TiterWaningM<- cbind(TiterWaningM, pred = predict(mod), preduni = predict(mod.uni))
ggplot(data = TiterWaningM)+ 
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
#ggsave("./output/figures/Waning.png")



#Paper Rudolf et al ####
TitreWaningR <- read_delim("./input/TitreWaningRudolf2010.csv", delim =";")
library(bbmle)

#H5N2
LL <- function(mu, vaccin.status = "BI"){
  dataToFit <- TitreWaningR%>%filter(VaccinStatus ==vaccin.status)
  -log(sum(dbinom(x = dataToFit$PosH5N2, size = dataToFit$NH5N2, p = exp(-mu*dataToFit$TimeWeeks))))
}
curve<-sapply(seq(0,0.05,0.001),LL)
plot(seq(0.,0.05,0.001),curve)
fit<- mle2(LL,c(mu =0.02))
fit
fit@fullcoef/7
with(TitreWaningR%>%filter(VaccinStatus =="BI"),plot(x = TimeWeeks,y = PosH5N2/NH5N2))
lines(seq(0.,max(TitreWaningR$TimeWeeks),1),exp(-fit@fullcoef*seq(0.,max(TitreWaningR$TimeWeeks),1)))

#H5N1
LL <- function(mu, vaccin.status = "BI"){
  dataToFit <- TitreWaningR%>%filter(VaccinStatus ==vaccin.status)
  -log(sum(dbinom(x = dataToFit$PosH5N1, size = dataToFit$NH5N1, p = exp(-mu*dataToFit$TimeWeeks))))
}
curve<-sapply(seq(0,0.05,0.001),LL)
plot(seq(0.,0.05,0.001),curve)
fit<- mle2(LL,c(mu =0.02))
fit
fit@fullcoef/7
with(TitreWaningR%>%filter(VaccinStatus =="BI"),plot(x = TimeWeeks,y = PosH5N1/NH5N1))
lines(seq(0.,max(TitreWaningR$TimeWeeks),1),exp(-fit@fullcoef*seq(0.,max(TitreWaningR$TimeWeeks),1)))

#fit generalized logistic model to the data ####
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
plot(seq(0, max((TitreWaningR%>%filter(VaccinStatus =="BI"))$TimeWeeks)) ,
     sapply(seq(0, max((TitreWaningR%>%filter(VaccinStatus =="BI"))$TimeWeeks)),genlogpar),"l", ylim=c(0,1))
points((TitreWaningR%>%filter(VaccinStatus =="BI"))$TimeWeeks ,(TitreWaningR%>%filter(VaccinStatus =="BI"))$PosH5N2/60 )

#fit gamma distribution ####
LL.gamma <- function(shape,rate, vaccin.status = "BI"){
  dataToFit <- TitreWaningR%>%filter(VaccinStatus ==vaccin.status&TimeWeeks>3)
  with(TitreWaningR%>%filter(VaccinStatus =="BI"),sum((NH5N2-PosH5N2)*log(pgamma(TimeWeeks,shape,rate))+PosH5N2*log((1-pgamma(TimeWeeks,shape,rate)))))
}

LL.gamma(50,.10)
pgamma(77,60,.5)

fit.gamma<- mle2(LL.gamma,start = list(shape = 70,rate =1.0),fixed = list(),method = "Nelder-Mead")
fit.gamma
fit.gamma@fullcoef[1]/fit.gamma@fullcoef[2]
fit.gamma@fullcoef[1]/fit.gamma@fullcoef[2]^2
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     plot(seq(0,max(TimeWeeks),1),1-pgamma(seq(0,max(TimeWeeks),1),fit.gamma@fullcoef["shape"], fit.gamma@fullcoef["rate"]),"l", ylim = c(0,1)))
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     points(TimeWeeks,PosH5N2/NH5N2))

#H5N1###
LL.gamma <- function(shape,rate, vaccin.status = "BI"){
  dataToFit <- TitreWaningR%>%filter(VaccinStatus ==vaccin.status&TimeWeeks>3)
  with(TitreWaningR%>%filter(VaccinStatus =="BI"),sum((NH5N1-PosH5N1)*log(pgamma(TimeWeeks,shape,rate))+PosH5N1*log((1-pgamma(TimeWeeks,shape,rate)))))
}

LL.gamma(50,.10)
pgamma(77,60,.5)

fit.gamma<- mle2(LL.gamma,start = list(shape = 70,rate =1.0),fixed = list(),method = "Nelder-Mead")
fit.gamma
fit.gamma@fullcoef[1]/fit.gamma@fullcoef[2]
fit.gamma@fullcoef[1]/fit.gamma@fullcoef[2]^2
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     plot(seq(0,max(TimeWeeks),1),1-pgamma(seq(0,max(TimeWeeks),1),fit.gamma@fullcoef["shape"], fit.gamma@fullcoef["rate"]),"l", ylim = c(0,1)))
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     points(TimeWeeks,PosH5N1/NH5N1))

#bayesian fit####
library(rstan)

scode <- "
data {
  int N;
  vector[N] Times;
  vector[N] left_censor;
  vector[N] right_censor;
}

parameters {
  real<lower=0> shape;
  real<lower=0> rate;
  }

model {
  shape ~ normal(60,10);
  rate ~ uniform(0.1,2.);
  target += right_censor*gamma_lccdf(Times | shape, rate) + left_censor*gamma_lcdf(Times | shape, rate);
}
"
dat<- with(TitreWaningR,list(N =length(TimeWeeks),Times = TimeWeeks, right_censor = PosH5N2,left_censor = NH5N2-PosH5N2)) 
fit <- stan(model_code=scode, data=dat, iter=1000, chains=3, verbose=F)
fit.shape <- extract(fit)$shape%>%median
fit.rate <-extract(fit)$rate%>%median
fit.shape
fit.rate
fit.shape/fit.rate
plot(fit)
stan_hist(fit)
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     plot(seq(0,max(TimeWeeks),1),1-pgamma(seq(0,max(TimeWeeks),1),fit.shape, fit.shape),"l", ylim = c(0,1)))
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     points(TimeWeeks,PosH5N2/NH5N2))


library(rstan)

scode <- "
data {
  int N;
  vector[N] Times;
  vector[N] left_censor;
  vector[N] right_censor;
}

parameters {
  real<lower=0, upper = 100> mu;
 real<lower=0> sigma;
  }

model {
  mu ~ normal(60,10);
  sigma ~ normal(2,2);
  target += right_censor*normal_lccdf(Times | mu, sigma) + left_censor*normal_lcdf(Times | mu, sigma);
}
"
dat<- with(TitreWaningR,list(N =length(TimeWeeks),Times = TimeWeeks, right_censor = PosH5N2,left_censor = NH5N2-PosH5N2)) 
fit <- stan(model_code=scode, data=dat,init = list(list(mu = 70)), iter=1000, chains=3, verbose=F)
fit.mu <- extract(fit)$mu%>%median
fit.sigma <-extract(fit)$sigma%>%median

plot(fit)
stan_hist(fit)
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     plot(seq(0,max(TimeWeeks),1),1-pnorm(seq(0,max(TimeWeeks),1),fit.mu, fit.sigma),"l", ylim = c(0,1)))
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     points(TimeWeeks,PosH5N2/NH5N2))

#all does not work - easier approach
#We see that at time 77 weeks approximately half of the birds are positive
TitreWaningR$fH5N2 <- TitreWaningR$PosH5N2/TitreWaningR$NH5N2
TitreWaningR$fH5N1 <- TitreWaningR$PosH5N1/TitreWaningR$NH5N1
#given an approximate normal distribution of the infectious periods and tweak it
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     plot(seq(0,max(TimeWeeks),1),1-pgamma(seq(0,max(TimeWeeks),1),30, 0.4),"l", ylim = c(0,1)))
with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     points(TimeWeeks,PosH5N2/NH5N2))


with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     plot(seq(0,max(TimeWeeks),1)*7,(1-pgamma(seq(0,max(TimeWeeks),1)*7,4, .1/7)),"l", ylim = c(0,1)))

with((TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)),
     points(TimeWeeks*7,(1/0.7) *PosH5N1/NH5N1))


#scale to days 
tmp <- TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)
data.fit <- rbind(data.frame(time = seq(0,max(tmp$TimeWeeks*7),1),prob = 1-pgamma(seq(0,max(tmp$TimeWeeks*7),1),36, 0.07), strain = "H5N2 - homologous"),
                  data.frame(time = seq(0,max(tmp$TimeWeeks*7),1),prob = .7*(1-pgamma(seq(0,max(tmp$TimeWeeks*7),1),4, 0.1/7)), strain = "H5N1 - heterologous"))
data.exp <- rbind(data.frame(time = tmp$TimeWeeks*7,prob = tmp$PosH5N2/tmp$NH5N2, strain = "H5N2 - homologous"),
                  data.frame(time = tmp$TimeWeeks*7,prob = tmp$PosH5N1/tmp$NH5N1, strain = "H5N1 - heterologous"))
ggplot()+ 
  geom_path(data =data.fit , aes(time, prob,colour = strain, linetype = strain, group = strain), linewidth = 1.5 )+
  geom_point(data = data.exp, aes(time, prob,colour = strain, shape = strain, group = strain),size = 4)+
  scale_colour_manual(values = c("H5N2 - homologous"="#F3965E", "H5N1 - heterologous" = "#AA1555"))+
  #scale_linetype_discrete(values = c("dash",  "solid"))+
  labs(x = "Time (days post full vaccination)", y = "Prob. high titre")+
  ylim(0,NA)
  # +ggtitle("Prob. of having a high titre")
ggsave("./output/figures/WaningCurve.png",scale =1.2)

#scale to days 
tmp <- TitreWaningR%>%filter(VaccinStatus =="BI"&TimeWeeks>3)
data.fit <- rbind(data.frame(time = seq(0,max(tmp$TimeWeeks*7),1),prob = 1-pgamma(seq(0,max(tmp$TimeWeeks*7),1),36, 0.07), strain = "Homologous"),
                  data.frame(time = seq(0,max(tmp$TimeWeeks*7),1),prob = .7*(1-pgamma(seq(0,max(tmp$TimeWeeks*7),1),4, 0.1/7)), strain = "Heterologous"))
data.exp <- rbind(data.frame(time = tmp$TimeWeeks*7,prob = tmp$PosH5N2/tmp$NH5N2, strain = "Homologous"),
                  data.frame(time = tmp$TimeWeeks*7,prob = tmp$PosH5N1/tmp$NH5N1, strain = "Heterologous"))
ggplot()+ 
  geom_path(data =data.fit , aes(time, prob,colour = strain, linetype = strain, group = strain), linewidth = 1.5 )+
  geom_point(data = data.exp, aes(time, prob,colour = strain, shape = strain, group = strain),size = 4)+
  scale_colour_manual(values = c("Homologous"="#F3965E", "Heterologous" = "#AA1555"))+
  #scale_linetype_discrete(values = c("dash",  "solid"))+
  labs(x = "Time (days post full vaccination)", y = "Prob. high titre")+
  ylim(0,NA)
# +ggtitle("Prob. of having a high titre")
ggsave("./output/figures/WaningCurve_poster.png",scale =1.2)


test<- rgamma(10000, 36,0.07)
mean(test)/7
var(test/7)
hist(test/7)
