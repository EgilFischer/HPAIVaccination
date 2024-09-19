# HPAIVaccination

Code for assessing the risk of zoonotic transmission towards humans after vaccination of poultry flocks

# Code
- TechReport.RMD: R markdown file with explanation of methods and models

Code for simulations can be found under ./src
There are two folders 
- WithinFarm: consists of files for within farm modelling
Most important files

* MultitypetransitionsSEIR_tleap.R: Stochasttic multitype model SEIR with transitions between types (currently used)
* probMajorOutbreak.R: Calculation of probability of major outbreak given vaccination coverage
* postprocessSimulations.R: process output of stochastic models to visualize, determine detection moment and human exposure

- Spatial: files for spatial modelling




