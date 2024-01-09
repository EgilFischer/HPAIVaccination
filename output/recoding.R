setwd("./output")
#replace file names
ofs<- list.files()
nfn<- data.frame(oldname =subset(ofs,grepl("20000",ofs)|grepl("38000",ofs)|grepl("73000",ofs)),
                 newname = gsub("layer","broiler",subset(ofs,grepl("20000",ofs)|grepl("38000",ofs)|grepl("73000",ofs))))
nfn
file.rename(nfn$oldname,nfn$newname)

#replace scenario names in the files

scenarios.tmp.rec <-  data.frame(
  type = c(rep("broiler",6*3)),
  size = rep(c(20000,38000,73000),each = 6),
  vac  = rep(c("Vac0","Vac50","Vac60","Vac70","Vac80","Vac90"),3))

scenario.list.size.type.rec <-unlist( mapply(FUN = function(type, size,vac){list(data.frame(scenario  = paste0(type,"Size",size,vac)))},
                                      scenarios.tmp.rec$type,scenarios.tmp.rec$size,scenarios.tmp.rec$vac))

#for all vaccination types
for(i in c(1:(length(scenario.list.size.type.rec)))){
  rec.data <- readRDS(paste0("./output/coverage/20231210output",scenario.list.size.type.rec[i],".RDS"))
  rec.data$pars$scenario<- gsub("layer","broiler", rec.data$pars$scenario)
  saveRDS(rec.data,paste0("./output/coverage/20231210output",scenario.list.size.type.rec[i],".RDS"))
}


#put is subdirectories
dir.create("clinic")
dir.create("coverage_broiler")
dir.create("coverage_layer")
dir.create("waning_homologous")
dir.create("waning_heterologous")