##################################################################
#  input weibull parameters and convert to erlang parameters     #
##################################################################

#Minimize the difference between the mean and variance of the Weibull distribution and Erlang distribution
brute_force_method <- function(shape, scale,kmin =1, kmax = 100) {
  # Calculate the mean and variance of the Weibull distribution
  mean_weibull <- scale * gamma(1 + 1 / shape)
  var_weibull <- scale^2 * (gamma(1 + 2 / shape) - (gamma(1 + 1 / shape))^2)
  
  # Iterate through possible Erlang shape parameters to search for minimal difference with Weibull
  # Erlang distribution has parameter shape k and parameter scale beta
  dif_vars = Inf;
  best_pars = list(shape = kmin, scale = mean_weibull)
  for(k in kmin:kmax) {
    # Calculate the Erlang parameters
    erlang_scale <- mean_weibull/k
    var_erlang <- k*erlang_scale^2
    
    # Check if the difference of the variances is smaller than of a previous set of parameters
    if (abs(var_weibull - var_erlang) < dif_vars) {
      dif_vars = abs(var_weibull - var_erlang)
      best_pars = list(shape = k, scale = erlang_scale)
    }
 
  }
  # If no better parameters found, return the best found so far
  return(best_pars)
}



#random sample based method return a gamma distribution
random_sample_method <- function(shape, scale, sample_size = 1000) {
  #create a random sample from the Weibull distribution
  weibull_sample <- rweibull(sample_size, shape, scale)
  
  #fit erlang distribution
  fit <- MASS::fitdistr(weibull_sample, "gamma")
  
  #return the parameters of the fitted gamma distribution
  return(list(shape = fit$estimate[1], scale = 1/fit$estimate[2]))
}

#random sample based method returning an erlang distribution
random_sample_erlang_method <- function(shape, scale, sample_size = 1000, kmin = 1, kmax = 10, start_val =NULL) {
  #create a random sample from the Weibull distribution
  weibull_sample <- rweibull(sample_size, shape, scale)

  #define likelihood function for Erlang distribution
  log_likelihood_function <- function(params, k) {
    scale_param <- params[1]
    -sum(dgamma(weibull_sample, shape = k, scale = scale_param, log = TRUE))
  }
  #fit erlang distribution
  LL = Inf; fit = NULL;erlang_shape =kmin;
  for(k in c(kmin:kmax))
  {
    #if no start value is given, use the mean of the sample
    if (is.null(start_val)) {
      start_val_use <- scale * gamma(1 + 1 / shape)/k
    } else {
      start_val_use <- start_val
    }
    
    fit_new <- optim(c(start_val_use), 
                     log_likelihood_function, 
                     k = k,
                     method = "L-BFGS-B", lower = c(1e-6, 1e-6), upper = c(Inf, Inf))
    if (fit_new$value < LL) {
      LL = fit_new$value
      fit = fit_new
      erlang_shape = k;
    }
  }
  
  
  #round the shape parameter to the nearest integer for Erlang distribution
  return(list(shape = erlang_shape, scale = fit$par[1]))
}




# Apply to the field data
## Call:
## survreg(formula = formula, data = db, dist = "weibull")
##                 Value Std. Error     z      p
## factor(HI_5)0  1.3895     0.0345  40.3 <2e-16
## factor(HI_5)1  1.5722     0.1048  15.0 <2e-16
## Log(scale)    -1.2213     0.0888 -13.8 <2e-16
## 
## Scale= 0.295#

# Apply the brute force method to the field data parameters
weibull_shape_field <- 1.3895
weibull_scale_field <- 1/exp(-1.2213)
# Print mean, median and variance of the Weibull distribution
cat("Weibull Distribution Parameters:\n")
cat("Mean:", weibull_scale_field * gamma(1 + 1 / weibull_shape_field), "\n")
cat("Median, 2.5%, 97.5%:\n", qweibull(c(0.5,0.025,0.975), shape = weibull_shape_field, scale =  weibull_scale_field), "\n")
cat("Variance:",weibull_scale_field^2 * (gamma(1 + 2 / weibull_shape_field) - (gamma(1 + 1 / weibull_shape_field))^2), "\n")



# Using brute force method for field data
erlang_params_field_brute <- brute_force_method(weibull_shape_field, weibull_scale_field, kmin =1, kmax = 5)
# Using random sample method for field data rounded gamma
round_gamma_params_field <- random_sample_method(weibull_shape_field, weibull_scale_field)
# Using ramdom sample method for field data erlang
erlang_params_field <- random_sample_erlang_method(weibull_shape_field, weibull_scale_field, kmin = 1, kmax = 5, start_val = 0.7)

# Compare the results
cat("Comparison of results:\n")
library(ggplot2)
time_vector <- seq(0, 10, by = 1)
comparison_df <- data.frame(
  Time = time_vector,
  Weibull = pweibull(time_vector, weibull_shape_field, weibull_scale_field),
  Erlang_Brute = pgamma(time_vector,erlang_params_field_brute$shape, scale = erlang_params_field_brute$scale),
  Erlang_Random = pgamma(time_vector, erlang_params_field$shape, scale =  erlang_params_field$scale),
  Erlang_Rounded = pgamma(time_vector, round(round_gamma_params_field$shape), scale = round_gamma_params_field$scale)
)

#plot the comparison
ggplot(comparison_df, aes(x = Time)) +
  geom_line(aes(y = Weibull, color = "Weibull")) +
  geom_line(aes(y = Erlang_Brute, color = "Erlang (Brute Force)")) +
  geom_line(aes(y = Erlang_Random, color = "Erlang (Random Sample)")) +
  geom_line(aes(y = Erlang_Rounded, color = "Erlang (Rounded Gamma)")) +
  labs(title = "Comparison of Weibull and Erlang Distributions",
       y = "Cumulative Probability",
       color = "Distribution") +
  theme_minimal()
# Print the parameters
cat("Erlang Parameters (Brute Force Method):\n")
cat("Shape:", erlang_params_field_brute$shape, "\n")
cat("Scale:", erlang_params_field_brute$scale, "\n")
cat("Erlang Parameters (Random Sample Method):\n")
cat("Shape:", erlang_params_field$shape, "\n")
cat("Scale:", erlang_params_field$scale, "\n")
cat("Erlang Parameters (Rounded Gamma Method):\n")
cat("Shape:", round(round_gamma_params_field$shape), "\n")
cat("Scale:", round_gamma_params_field$scale, "\n")

#input into simulations
mean_infT_Low <-erlang_params_field_brute$shape *erlang_params_field_brute$scale


# vaccine group ###
# Apply the brute force method to the field data parameters
weibull_shape_field <- 1.5722
weibull_scale_field <- 1/exp(-1.2213)
# Print mean, median and variance of the Weibull distribution
cat("Weibull Distribution Parameters:\n")
cat("Mean:", weibull_scale_field * gamma(1 + 1 / weibull_shape_field), "\n")
cat("Median, 2.5%, 97.5%:\n", qweibull(c(0.5,0.025,0.975), shape = weibull_shape_field, scale =  weibull_scale_field), "\n")
cat("Variance:",weibull_scale_field^2 * (gamma(1 + 2 / weibull_shape_field) - (gamma(1 + 1 / weibull_shape_field))^2), "\n")



# Using brute force method for field data
erlang_params_field_brute <- brute_force_method(weibull_shape_field, weibull_scale_field, kmin =1, kmax = 5)
# Using random sample method for field data rounded gamma
round_gamma_params_field <- random_sample_method(weibull_shape_field, weibull_scale_field)
# Using ramdom sample method for field data erlang
erlang_params_field <- random_sample_erlang_method(weibull_shape_field, weibull_scale_field, kmin = 1, kmax = 5)

# Compare the results
cat("Comparison of results:\n")
library(ggplot2)
time_vector <- seq(0, 10, by = 1)
comparison_df <- data.frame(
  Time = time_vector,
  Weibull = pweibull(time_vector, weibull_shape_field, weibull_scale_field),
  Erlang_Brute = pgamma(time_vector,erlang_params_field_brute$shape, scale = erlang_params_field_brute$scale),
  Erlang_Random = pgamma(time_vector, erlang_params_field$shape, scale = erlang_params_field$scale),
  Erlang_Rounded = pgamma(time_vector, round(round_gamma_params_field$shape), scale = round_gamma_params_field$scale)
)

#plot the comparison
ggplot(comparison_df, aes(x = Time)) +
  geom_line(aes(y = Weibull, color = "Weibull")) +
  geom_line(aes(y = Erlang_Brute, color = "Erlang (Brute Force)")) +
  geom_line(aes(y = Erlang_Random, color = "Erlang (Random Sample)")) +
  geom_line(aes(y = Erlang_Rounded, color = "Erlang (Rounded Gamma)")) +
  labs(title = "Comparison of Weibull and Erlang Distributions",
       y = "Cumulative Probability",
       color = "Distribution") +
  theme_minimal()
# Print the parameters
cat("Erlang Parameters (Brute Force Method):\n")
cat("Shape:", erlang_params_field_brute$shape, "\n")
cat("Scale:", erlang_params_field_brute$scale, "\n")
cat("Erlang Parameters (Random Sample Method):\n")
cat("Shape:", erlang_params_field$shape, "\n")
cat("Scale:", erlang_params_field$scale, "\n")
cat("Erlang Parameters (Rounded Gamma Method):\n")
cat("Shape:", round(round_gamma_params_field$shape), "\n")
cat("Scale:", round_gamma_params_field$scale, "\n")

#input into simulations
mean_infT_High <-erlang_params_field_brute$shape *erlang_params_field_brute$scale
median_infT_High <- qgamma(0.5, shape = erlang_params_field_brute$shape, scale = erlang_params_field_brute$scale)

#print input
cat("Input for simulations:\n")
cat("Mean InfT Low:", mean_infT_Low, "\n")
cat("Mean InfT High:", mean_infT_High, "\n")
cat("Scale:",erlang_params_field_brute$shape)
