#correct error with first line in simulation having wrong run number
folder <- "./output/baseline_broiler"
sims <- list.files(folder, pattern = ".RDS")
for(s in sims){
  rds <- readRDS(paste0(folder,"/",s))
  rds$out$run <-matrix(rep(max(rds$out$run),length(rds$out$time[,1])),ncol =1);
  saveRDS(rds,file = paste0(folder,"/",s)) 
}

