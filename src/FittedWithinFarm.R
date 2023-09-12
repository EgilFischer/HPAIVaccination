#######################################################
#   Fitted distributions for the scenarios            #
#######################################################
#gamma distributions shape = mean^2/var, rate = mean/var

####T dist ####
#no vaccination ####
Tdist.novac.pas = function(vacstat,size, type,fadeout,intro.time){
  mean = 8.3+1.33E-05*size-0.01*as.numeric(type=="LAYER"); 
  var = (8.282-1.03E-05 * size + 9.39E-03*as.numeric(type=="LAYER"));
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

#vaccination no waning ###
Tdist.vac.pas = function(vacstat,size, type,fadeout,intro.time){
  if(type !="LAYER") return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  mean = 8.3+0.21992*vacstat+3.22E-05*size; var = (1.055-1.42E-05 * size + 4.15E-02*vacstat)^2;
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

#clinical protection ####
Tdist.clinprot.pas = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat==0 || type !="LAYER")return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  mean = 9.4; 
  var = 0.23;
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

#waning homologuous ####
Tdist.vac.wane.pas = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat == 0) return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  if(intro.time>= 500) {mean = 18.1; var = 9.47} 
  else {mean = 8; var = 1}#fade out
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

#waning heterologuous ####
Tdist.vac.wane.ds.pas = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat < 10^-10) return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  if(intro.time>= 50) {mean = 10^(1.60942-0.00154*intro.time); var = (3.781674-0.0064*intro.time)^2}
  else {mean = 8; var =1}#fade out
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}


####Human exposure ####
#no vaccination ####
exposure.function.novac.pas <- function(sim.data){
  return(-3034.545 + as.numeric(sim.data$SIZE)*1.287 - (584.923*(sim.data$TYPE == "LAYER")))
}

#vaccination no waning ####
exposure.function.vac.pas <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){
    if(sim.data$TYPE[x]!="LAYER"){return(exposure.function.novac.pas(sim.data[x,]))}else{
  if(sim.data[x,]$vac <= 0.50)return(10^(4.11 + sim.data[x,]$SIZE*1.22E-05 -1.36E-02*sim.data[x,]$vac*100)) 
  if(sim.data[x,]$vac <= 0.70)return(10^(3.347 + sim.data[x,]$SIZE*1.37E-05)) 
  if(sim.data[x,]$vac <= 1.0)return(10^(13.7 + sim.data[x,]$SIZE*7.25E-06 +-1.47E-01*sim.data[x,]$vac*100)) }})
  
}

#clinical protection ####
exposure.function.clinprot.pas <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){
  if(sim.data[x,]$vac >0)  return(52735.995) else
    return(exposure.function.novac.pas(sim.data[x,]))})
}

#waning homologuous ####
exposure.function.wane.pas <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){if(sim.data[x,]$vac>0){
    if(sim.data[x,]$intro.time <400) {return(2.41872)}else {
      -27053.2+sim.data[x,]$intro.time*67.66
    }}
    else{    return(exposure.function.novac.pas(sim.data[x,]))}
  })
}

#waning heterologuous ####
exposure.function.wane.ds.pas <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){if(sim.data[x,]$vac>0){
    10^(3.81+0.00131*sim.data[x,]$intro.time)}
    else{    
      return(exposure.function.novac.pas(sim.data[x,]))}
  })
}

####Fade out ####
#no vaccination ####
fadeout_func.novac.pas = function(vacstat){
  return(0);
}
#vaccination no waning ####
fadeout_func.vac.pas = function(vacstat){
  return(0);
}

#clinical protection ####
fadeout_func.clinprot.pas = function(vacstat){
  return(0);
}
#waning homologuous ####
fadeout_func.wane.pas = function(vacstat){
  return(0);
}
#waning heterologuous ####
fadeout_func.wane.ds.pas = function(vacstat){
  return(0);
}
#test#####
# set.seed(1)
# Tdist.novac.pas(0.5,32000, "LAYER",0,0)
# set.seed(1)
# Tdist.novac.pas(0.5,32000, "BROILER",0,0)
# set.seed(1)
# Tdist.vac.pas(0.5,32000, "BROILER",0,0)
# set.seed(1)
# Tdist.vac.pas(0.5,32000, "LAYER",0,0)
# 
# 
# exposure.function.novac.pas(data.frame(SIZE = 32000, vac =c(0.25, 0.6,0.8), TYPE = "LAYER"))
# exposure.function.vac.pas(data.frame(SIZE = 32000, vac =c(0.25, 0.6,0.8), TYPE = "LAYER"))
# exposure.function.novac.pas(data.frame(SIZE = 32000, vac =c(0.25, 0.6,0.8), TYPE = "BROILER"))
# exposure.function.vac.pas(data.frame(SIZE = 32000, vac =c(0.25, 0.6,0.8), TYPE = "BROILER"))
