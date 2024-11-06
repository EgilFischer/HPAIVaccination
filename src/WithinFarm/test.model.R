# param.list$intro.time <-0
# param.list$beta <- 0*param.list$beta
# param.list$latency.period <- 1000
# param.list$trans.mean.wane <- 14
  
rm(test.sim)
test.sim <- sim.multitypeSEIR_tleap_buildup_distwaning(param.list, 
                                           inits.gamma.buildup.wane,
                                           gamma.buildup.distribution,
                                           gamma.waning.distribution)
out <- data.frame(cbind(test.sim$time, test.sim$run,test.sim$S,test.sim$L,test.sim$I))

ggplot(out)+geom_path(aes(X1,X3,colour = factor(X2)))+ggtitle("Low titre")

ggplot(out)+geom_path(aes(X1,X4,colour = factor(X2)))+ggtitle("high titre")

ggplot(out)+geom_path(aes(X1,X5,colour = factor(X2)))

ggplot(out)+geom_path(aes(X1,log10(X6),colour = factor(X2)))
 #all infected
ggplot(out)+geom_path(aes(X1,log10(X5+X6+X7+X8),colour = factor(X2)))+geom_hline(yintercept=log10(1))

