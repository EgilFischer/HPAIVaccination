# The model used by Elisa Beninca is no longer available on CRAN. 

#just stupidly simulated x and y coordinates
n = 50
results <- data.frame(x = sample(1:200,n, replace = F),
                      y = sample(1:200,n, replace = F))

plot(results)
save(results,file = "configuration_points_200X200_whittle.Rdata")

n.points <- 2000
scalepar_values <- c(0.1, 2, 4, 8, 16)
varpar_values   <- c(1, 2, 4, 8, 16)

# Initalize a matrix of zeros where you store the number of points
points_conf <- as.data.frame(
  matrix(0,
         nrow = n.points,
         ncol = 2))#length(scalepar_values)*length(varpar_values)*2)) # Multiply for 2 because you have x and y coordinates
points_conf[1]<- runif(n.points, min = 0, max = 140)
points_conf[2]<- runif(n.points, min = 300, max = 325)
plot(points_conf)
# # Load packages
# library(sp)
# library(RandomFields) #not available any more
# library(gridExtra)
# library(rprojroot)
# library(here)
# library(RColorBrewer)
# 
# # Set spConform to FALSE to enable faster simulation
# # Set seed to NA to make seed arbitrary
# # See help(RFoptions)
# RlFoptions(spConform = FALSE, seed = NA)
# 
# # Set the total number of points
# n.points <- 2000
# 
# # First choose a covariance model. See help(RMmodel) for the options
# # We choose the RMstable model here. See help(RMstable)
# # Parameters:
# #   var = variance of the random field
# #   scale = correlation distance. Note: the effective correlation distance depends on the model!
# #   alpha = 'smoothness' parameter (0, 2]. alpha = 1 is exponential decay
# scalepar_values <- c(0.1, 2, 4, 8, 16)
# varpar_values   <- c(1, 2, 4, 8, 16)
# 
# # Initalize a matrix of zeros where you store the number of points
# points_conf <- as.data.frame(
#   matrix(0,
#          nrow = n.points,
#          ncol = length(scalepar_values)*length(varpar_values)*2)) # Multiply for 2 because you have x and y coordinates
# 
# # Choose one set of scale and variance parameters
# scalepar <- scalepar_values[4]
# varpar   <- varpar_values[4]
# 
# # Covariance model
# cov.model <- RMwhittle(var = varpar, scale = scalepar, nu = 1)
# 
# # Set the locations in the x and y direction
# # We make a rectangular grid in the unit square
# # (also used for plotting the results later on)
# x.seq <- seq(from = 0, to = 200, by = 1)
# y.seq <- seq(from = 0, to = 200, by = 1)
# 
# # Calculate the mean intensity
# # This is the mean number of points per grid cell: n.points/n.cells
# # However, because of the exponent in the rf term, we should divide by the smearing factor: exp(0.5*var)
# n.cells <- length(x.seq)*length(y.seq)
# mean.intensity <- n.points/n.cells/exp(0.5*cov.model@par.general$var)
# 
# # Set up a dataframe for all results
# # We use the expand.grid function from the base package
# result.data <- expand.grid(x = x.seq, y = y.seq)
# 
# # Fix seed
# set.seed(1)
# 
# # Add this to result.data
# result.data <- within(result.data, {
#   # Generate a realisation of the random field
#   rf <- RFsimulate(model = cov.model, x = x, y = y)
#   # Calculate log intensity
#   log.intensity <- log(mean.intensity) + rf
#   # Intensity (needs to be positive)
#   intensity <- exp(log.intensity)
# })
# 
# # Add this to result.data
# result.data <- within(result.data, {
#   # Generate a realisation of the random field
#   rf <- RFsimulate(model = cov.model, x = x, y = y)
#   # Calculate log intensity
#   log.intensity <- log(mean.intensity) + rf
#   # Intensity (needs to be positive)
#   intensity <- exp(log.intensity)
# })
# 
# # The sum of all intensities is approximately equal to n.points
# # We take care of this later
# sum(result.data$intensity)
# 
# # Plot log(intensity)
# z.breaks <- pretty(result.data$log.intensity, n = 10)
# cols <- colorRampPalette(colors = brewer.pal(n = 9, name = "Blues"))(n = length(z.breaks) - 1)
# par(mar = c(4.5, 4.5, 0.5, 0.5))
# image(x = x.seq, y = y.seq,
#       z = matrix(result.data$log.intensity, nrow = length(x.seq), ncol = length(y.seq)),
#       breaks = z.breaks, col = cols)
# 
# # Simulate inhomogeneous Poisson process
# # For each grid cell, draw a number of points from a Multinomial distribution with a given probability vector
# # We use the link between the Poisson distribution and the Multinomial distribution to condition on the total n.points
# # The probability = intensity/sum(intensity).
# # Within the grid cell, the location of the points is uniform in space
# result.data <- within(result.data, {
#   # The number of points within each cell
#   n <- rmultinom(n = 1, size = n.points, prob = intensity/sum(intensity))
# })
# 
# # What is the size of cells?
# x.dim <- diff(x.seq)[1]
# y.dim <- diff(y.seq)[1]
# 
# # Draw points
# points.list <- with(result.data, mapply(
#   x = x, y = y, n = n,
#   FUN = function(x, y, n) {
#     # Only draw points if n > 0
#     if (n > 0) {
#       return(cbind(
#         x = runif(n = n, min = x - x.dim/2, max = x + x.dim/2),
#         y = runif(n = n, min = y - y.dim/2, max = y + y.dim/2)))
#     }
#   }))
# 
# points_Z.list <- with(result.data, mapply(
#   x = x, y = y, n = n,
#   FUN = function(x, y, n) {
#     # Only draw points if n > 0
#     if (n > 0) {
#       return(cbind(
#         x_Z = runif(n = n, min = x - x.dim/2, max = x + x.dim/2),
#         y_Z = runif(n = n, min = y - y.dim/2, max = y + y.dim/2)))
#     }
#   }))
# points.mat <- do.call(what = "rbind", args = points.list)
# 
# # Plot the results
# plot(points.mat, pch = ".") 