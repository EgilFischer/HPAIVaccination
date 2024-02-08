#Labels for plotting #

itype.labels <- c("I.1" ="I.1",
                  "I.2" = "I.2",
                  "R.1"="R.1",
                  "R.2"="R.2",
                  "DR.1"= "Dead 1",
                  "DR.2" = "Dead 2")
SIR_labels <- c("S.1" ="S",
                "I.1"="I",
                "R.1"="R",
                "DR.1"= "Dead",
                "DI.1" = "Dead")
labels.type <-c("46"  = "Broiler", "510" = "Layer")
labels.size <-c("15000"  = "Small", "32000" = "Medium", "64000" = "Large", "20000"  = "Small", "38000" = "Medium", "73000" = "Large")
scenario.label.clin <- c(layerSize32000Vac0 = "Baseline",
                         layerClinic50Phigh50 = "Clinic 50% \n High titre 50%",
                         layerClinic50Phigh0 = "Clinic 50% \n High titre 0%",
                         layerClinic0.1Phigh0 = "Clinic 0.1% \n High titre 0%")


scenario.label.wane <- c(layerSize32000Vac0 = "Baseline",
                         layerStartTime0 = "T0 = 0",
                         layerStartTime50 = "T0 =  50",
                         layerStartTime100 = "T0 = 100",
                         layerStartTime200 = "T0 = 200",
                         layerStartTime300 = "T0 = 300",
                         layerStartTime400 = "T0 = 400",
                         layerStartTime500 = "T0 = 500")

scenario.levels.label.wane <- names(scenario.label.wane)

scenario.label.waneDifStrain <- c(layerSize32000Vac0 = "Baseline",
                                  layerStartTime0DifStrain = "T0 = 0",
                                  layerStartTime50DifStrain = "T0 =  50",
                                  layerStartTime100DifStrain = "T0 = 100",
                                  layerStartTime200DifStrain = "T0 = 200",
                                  layerStartTime300DifStrain = "T0 = 300",
                                  layerStartTime400DifStrain = "T0 = 400",
                                  layerStartTime500DifStrain = "T0 = 500")

scenario.levels.label.waneDifStrain <- names(scenario.label.waneDifStrain)



detection.method.label <- c(pas.det.time = "Passive",
                            ac.det.time = "Active",
                            min.det.time = "Minimum")

