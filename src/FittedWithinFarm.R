#######################################################
#   Fitted distributions for the scenarios            #
#######################################################
#gamma distributions shape = mean^2/var, rate = mean/var

####T dist ####
#no vaccination ####
Tdist.novac.pas = function(vacstat,size, type,fadeout,intro.time){
  mean = 8.3+1.33E-05*size-0.01*as.numeric(type=="LAYER"); 
  var = 0.5;#choose value that was in the middle range 
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}


Tdist.novac.min = function(vacstat,size, type,fadeout,intro.time){
  mean = 8.13+1.33E-05*size-0.01*as.numeric(type=="LAYER"); 
  var = 1.75;
  if(mean <0 )stop("Negative Tinf")
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

#vaccination no waning ###
Tdist.vac.pas = function(vacstat,size, type,fadeout,intro.time){
  if(type !="LAYER") return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  mean = 8.3+0.21992*vacstat+3.22E-05*size; var = (1.055-1.42E-05 * size + 4.15E-02*vacstat)^2;
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

Tdist.vac.min = function(vacstat,size, type,fadeout,intro.time){
  if(type !="LAYER") return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  mean = 8.91 +0.21992*vacstat-3.22E-05*size; var = (1.055-1.42E-05 * size + 4.15E-02*vacstat)^2;
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}


#clinical protection ####
Tdist.clinprot.pas = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat==0 || type !="LAYER")return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  mean = 9.4+0.21992*vacstat+3.22E-05*size; 
  var = 0.23;
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

Tdist.clinprot.min = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat==0 || type !="LAYER")return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  mean = 8.97+0.21992*vacstat+3.22E-05*size; 
  var = 0.64;
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

#waning homologuous ####
Tdist.vac.wane.pas = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat == 0) return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  if(intro.time>= 500) {mean = 18.1; var = 9.47} 
  else {mean = 8; var = 1}#fade out
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

Tdist.vac.wane.min = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat == 0) return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  if(intro.time>= 400) {mean = 10.202914; var = 3.6781475} 
  else {mean = 8; var = 1}#fade out
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}


#waning heterologuous ####
Tdist.vac.wane.ds.pas = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat == 0) return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  if(intro.time>= 50) {mean = 10^(1.60942-0.00154*intro.time)+0.21992*vacstat+3.22E-05*size; var = (3.781674-0.0064*intro.time)^2}
  else {mean = 8; var =1}#fade out
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}

Tdist.vac.wane.ds.min = function(vacstat,size, type,fadeout,intro.time){
  if(vacstat == 0) return(Tdist.novac.pas(vacstat,size, type,fadeout,intro.time))
  if(intro.time>= 50) {mean =11.1315912-0.0052078*intro.time+0.21992*vacstat+3.22E-05*size;
var = (2.5542855-0.0037109*intro.time)^2}
  else {mean = 8; var =1}#fade out
  return(rgamma(1, shape= mean^2/var, rate = mean/var))
}


####Human exposure ####
#no vaccination ####
exposure.function.novac.pas <- function(sim.data){if(sim.data$SIZE <= 2357.844) {return(return(2.41872))}else{
  return(-3034.545 + as.numeric(sim.data$SIZE)*1.287 - (584.923*(sim.data$TYPE == "LAYER")))}
}

exposure.function.novac.min <- function(sim.data){
  return(-825.742 + as.numeric(sim.data$SIZE)*1.055 - (-70.609*(sim.data$TYPE == "LAYER")))
}

#vaccination no waning ####
exposure.function.vac.pas <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){
    if(sim.data$TYPE[x]!="LAYER"){return(exposure.function.novac.pas(sim.data[x,]))}else{
  if(sim.data[x,]$vac <= 0.50)return(10^(4.11 + sim.data[x,]$SIZE*1.22E-05 -1.36E-02*sim.data[x,]$vac*100)) 
  if(sim.data[x,]$vac <= 0.70)return(10^(3.347 + sim.data[x,]$SIZE*1.37E-05)) 
  if(sim.data[x,]$vac <= 1.0)return(10^(13.7 + sim.data[x,]$SIZE*7.25E-06 +-1.47E-01*sim.data[x,]$vac*100)) }})
  
}

exposure.function.vac.min <- function(sim.data){
  0.64 * 10^(4.53+sim.data$size*5.62E-06 -3.96E-02 * sim.data$vac*100) #0.64 is scaling to match no vaccination
  }

#clinical protection ####
exposure.function.clinprot.pas <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){
  if(sim.data[x,]$vac >0)  return(52735.995) else
    return(exposure.function.novac.pas(sim.data[x,]))})
}


exposure.function.clinprot.min <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){
    if(sim.data[x,]$vac >0)  return(31983.922) else
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

exposure.function.wane.min <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){if(sim.data[x,]$vac>0){
    if(sim.data[x,]$intro.time <400) {return(2.41872)}else {
      -4072.2 +sim.data[x,]$intro.time*10.2
    }}
    else{    return(exposure.function.novac.pas(sim.data[x,]))}
  })
}

#waning heterologuous ####
exposure.function.wane.ds.pas <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){if(sim.data[x,]$vac>0){
    10^(3.81+0.00131*sim.data[x,]$intro.time)}
    else{    return(exposure.function.novac.pas(sim.data[x,]))}
  })
}


exposure.function.wane.ds.min <- function(sim.data){
  sapply(c(1:length(sim.data[,1])),FUN = function(x){if(sim.data[x,]$vac>0){
    if(sim.data[x,]$intro.time>50){
    max(-9218.703    +70.254*sim.data[x,]$intro.time, 2.41872)}else{2.41872}}
    else{    return(exposure.function.novac.pas(sim.data[x,]))}
  })
}

####Fade out ####
fadeout_func.generic = function(vacstat){
  return(0);
}

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
