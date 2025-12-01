# Get waning distribution parameters from aggregate values ####
get_distribution_params <- function(agg_data, waning_model = "poisson") {
  # Check if the waning model is supported
  if (!waning_model %in% c("poisson", "negbin" ,"normal")) {
    stop("Unsupported waning model. Choose 'poisson','negbin' or 'normal'.")
  }
  # Check if the aggregate data has the required columns
  if (!all(c("Mean",	"Median",	"Min", 	"Max",	"Range",	" STDEV", 	"FracPos") %in% colnames(agg_data))) {
    warning("Aggregate data must contain the following columns: Mean, Median, Min, Max, Range, STDEV, FracPos.")
  }
  # create empty params list to return
  distribution <- list()
  # Fit the model to the aggregate data
  fit <- tryCatch({
    if (waning_model == "poisson") {
      #set parameter for poisson
    distribution$params <- list(lambda = agg_data$Mean)
      #determine fit to other values
    distribution$diagnostics <- list(median = (agg_data$Median -(agg_data$Mean + 1/3-1/(50*agg_data$Mean)))^2,
                                       min = (agg_data$Min-qpois(0.01,agg_data$Mean))^2,
                                       max = (agg_data$Max-qpois(0.99,agg_data$Mean))^2,
                                       stdev = (agg_data$STDEV^2 -agg_data$Mean)^2,
                                       frac_pos = (agg_data$FracPos - ppois(3, agg_data$Mean,lower.tail = FALSE)))
      
    } else if (waning_model == "negbin") {
      # Fit a negative binomial distribution
      p_est = agg_data$Mean/agg_data$STDEV^2
      r_est = agg_data$Mean^2 / (agg_data$STDEV^2 - agg_data$Mean)
      distribution$params <- list(p = p_est, 
                                  r = r_est)
      #determine fit to other values
      distribution$diagnostics <- list(median = NA,
                                       min = (agg_data$Min-qnbinom(0.01,r_est,p_est))^2,
                                       max = (agg_data$Max-qnbinom(0.99,r_est,p_est))^2,
                                       stdev = (agg_data$STDEV^2 - r_est * (1-p_est)/p_est^2)^2,
                                       frac_pos = (agg_data$FracPos - pnbinom(3, r_est,p_est,lower.tail = FALSE)))
      
   
    } else if (waning_model == "normal") {
      # Fit a normal distribution
      mu_est = agg_data$Mean
      sigma_est = agg_data$STDEV
      distribution$params <- list(mu = mu_est, 
                                  sigma = sigma_est)
      #determine fit to other values
      distribution$diagnostics <- list(median = (agg_data$Median - mu_est)^2,
                                       min = (agg_data$Min - qnorm(0.01, mu_est, sigma_est))^2,
                                       max = (agg_data$Max - qnorm(0.99, mu_est, sigma_est))^2,
                                       stdev = (agg_data$STDEV^2 - sigma_est^2)^2,
                                       frac_pos = (agg_data$FracPos - pnorm(3, mu_est, sigma_est, lower.tail = FALSE)))
    }
  }, error = function(e) {
    stop("Error fitting the distribution: ", e$message)
  })
  
  # Return the fitted parameters
  return(distribution)
}


data.waning = read_xlsx("input/Waning_immunity_fieldstudy.xlsx")
data.waning$FracPos <- data.waning$FracPos/100 # convert to fraction
data.waning$day <- as.numeric(sapply(c(data.waning["Dag van monstername"]),FUN = gsub,pattern ="D",replacement = ""))
#visualize data
ggplot(data.waning) +
  geom_point(aes(day, Min),colour ="lightblue") +
  geom_point(aes(day, Max),colour ="lightblue") +
  geom_point(aes(day, Median),colour ="orange") +
  geom_point(aes(day, Mean),colour = "red") +
  geom_point(aes(day, Mean - 1.96*STDEV),colour = "red") +
  geom_point(aes(day, Mean + 1.96*STDEV),colour = "red") +
  labs(title = "Waning Immunity Over Time", x = "Day", y = "Titre") +
  theme_minimal()


# get distribution parameters for different models ####
fp<- get_distribution_params(data.waning, waning_model = "poisson")
fnb<- get_distribution_params(data.waning, waning_model = "negbin")
fn<- get_distribution_params(data.waning, waning_model = "normal")

# Plot the expected fraction below 5
expected.waning <- data.frame(day =data.waning$day,
                              farm = data.waning$Locatie,
                              mu = fn$params$mu,
                              sigma = fn$params$sigma,
                              frac_pos_data = data.waning$FracPos,
                              frac_pos_3 = pnorm(3, fn$params$mu, fn$params$sigma, lower.tail = FALSE),
                              frac_pos_5 = pnorm(5, fn$params$mu, fn$params$sigma, lower.tail = FALSE),
                              frac_pos_6 = pnorm(6, fn$params$mu, fn$params$sigma, lower.tail = FALSE),
                              frac_pos_7 = pnorm(7, fn$params$mu, fn$params$sigma, lower.tail = FALSE))
ggplot(expected.waning)+
  geom_point(aes(day, frac_pos_data, colour = farm, shape = farm)) +
  geom_line(aes(day, frac_pos_3), colour = "blue", linetype = "dashed") +
  geom_line(aes(day, frac_pos_5), colour = "red", linetype = "dashed") +
  geom_line(aes(day, frac_pos_6), colour = "orange", linetype = "dashed") +
  geom_line(aes(day, frac_pos_7), colour = "pink", linetype = "dashed") + 
  scale_y_continuous(transform ="probit") 
  


# Get waning survival distribution parameters fraction positive at a certain day ####
