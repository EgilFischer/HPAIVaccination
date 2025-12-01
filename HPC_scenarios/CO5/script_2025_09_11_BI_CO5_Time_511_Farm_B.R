# Scenario:BI_CO5_Time_511_Farm_B
 param_list <- list(scenario="BI_CO5_Time_511_Farm_B",
runs=100,
deltat=0.005,
length.round=90,
age_at_d0=511,
N0=46000,
itypes=2,
p.protect=0.95,
no.introduction=1,
introduction.at.random=TRUE,
intro.time=0,
beta=matrix(c(2.13, 2.13, 0.21, 0.21), ncol = 2),
latency.period=c(1, 1),
k.latency=2,
infectious.period=c(4.42, 3.4),
k.infectious=2,
trans.mean.wane=Inf,
trans.var.wane=Inf,
trans.mean.buildup=Inf,
trans.var.buildup=Inf,
pdie=c(0.2, 0.01),
mortRate=2e-04,
eh=0.57,
ei=0.285,
disfigured=0.1,
pickup_time=3) 
# Source the simulator
source("/home/uu_vet_te/efischer/HPAI_vaccination/src/faster_multitypetransitionsSEIR_tleap.R")
# Run the simulator with the parameter set
sim.out <- sim.multitypeSEIR_tleap(param_list, seed = 12345)
op <- list(out = sim.out, pars = param_list)
saveRDS(op, file = paste0("/home/uu_vet_te/efischer/HPAI_vaccination/output/", format(Sys.Date(), "%Y_%m_%d_"), param_list$scenario, ".RDS"))
