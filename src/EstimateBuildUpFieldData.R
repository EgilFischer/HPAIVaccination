# STAN code to estimate the parameters of a gamma function
library(rstan)

#source ProcessFieldExperiment.R
source("./src/ProcessFieldExperiment.R")

# Define the Stan model
stan_model <- "
data {
  int<lower=0> N; // Number of observations
  int<lower=0, upper=1> y[N]; // Binary outcomes
  real<lower=0> time_point[N]; // Time points
}

parameters {
  real<lower=0> shape; // Shape parameter of the gamma distribution
  real<lower=0> rate; // Rate parameter of the gamma distribution
}

model {
  vector[N] gamma_cdf;

  // Calculate the CDF of the gamma distribution for each time point
  for (n in 1:N) {
    gamma_cdf[n] = gamma_cdf(time_point[n] | shape, rate);
  }

  // Likelihood
  y ~ bernoulli(gamma_cdf);

  // Priors
  shape ~ gamma(1, 1);
  rate ~ gamma(1, 1);
}
"

# Write the Stan model to a file
writeLines(stan_model, "gamma_model.stan")

# Example data: replace this with your actual data
create_surv_data <- function(input_data){list(
  N = length(input_data$time_point),
  y = as.integer(input_data$positive),
  time_point = input_data$time_point
)}

surv_data_field_study%>%filter(Treatment == "CEVA" & time_point<= 120)%>%create_surv_data() -> input_data_CEVA
surv_data_field_study%>%filter(Treatment == "BI" & time_point<= 120&!is.na(positive))%>%create_surv_data() -> input_data_BI
surv_data_field_study%>%filter(Treatment == "BI+booster" & time_point<= 120)%>%create_surv_data() -> input_data_BIBOOSTER


# Fit the model
fit_CEVA <- stan(
  file = "gamma_model.stan",
  data = input_data_CEVA,
  iter = 2000,
  chains = 4
)

# Print the summary of the model
print(fit_CEVA)
saveRDS(fit_CEVA, "./input/FieldExperiment/fit_CEVA.rds")

# Plot the posterior distributions of the parameters
plot(fit_CEVA,plotfun = "hist", pars = c("shape", "rate")) +
  labs(title = "Posterior Distributions for CEVA")
plot(fit_CEVA,plotfun = "trace", pars = c("shape", "rate")) +
  labs(title = "Trace for CEVA")

CEVA_shape <- median(extract(fit_CEVA, "shape")[[1]])
CEVA_rate <- median(extract(fit_CEVA, "rate")[[1]])
CEVA_CDF <- function(t) {
  pgamma(t, shape = CEVA_shape, rate = CEVA_rate)
}

# Fit the model for BI
fit_BI <- stan(
  file = "gamma_model.stan",
  data = input_data_BI,
  iter = 2000,
  chains = 4
)
# Print the summary of the model  
print(fit_BI)
saveRDS(fit_BI, "./input/FieldExperiment/fit_BI.rds")

# Plot the posterior distributions of the parameters
plot(fit_BI,plotfun = "hist", pars = c("shape", "rate")) +
  labs(title = "Posterior Distributions for BI")
plot(fit_BI,plotfun = "trace", pars = c("shape", "rate")) +
  labs(title = "Trace for BI")

BI_shape <- median(extract(fit_BI, "shape")[[1]])
BI_rate <- median(extract(fit_BI, "rate")[[1]])
BI_CDF <- function(t) {
  pgamma(t, shape = BI_shape, rate = BI_rate)
}

# Fit the model for BIBOOSTER
fit_BIBOOSTER <- stan(
  file = "gamma_model.stan",
  data = input_data_BIBOOSTER,
  iter = 2000,
  chains = 4
)
# Print the summary of the model
print(fit_BIBOOSTER)
saveRDS(fit_BIBOOSTER, "./input/FieldExperiment/fit_BIBOOSTER.rds")

# Plot the posterior distributions of the parameters
plot(fit_BIBOOSTER,plotfun = "hist", pars = c("shape", "rate")) +
  labs(title = "Posterior Distributions for BIBOOSTER")
plot(fit_BIBOOSTER,plotfun = "trace", pars = c("shape", "rate")) +
  labs(title = "Trace for BIBOOSTER")

BIBOOSTER_shape <- median(extract(fit_BIBOOSTER, "shape")[[1]])
BIBOOSTER_rate <- median(extract(fit_BIBOOSTER, "rate")[[1]])
BIBOOSTER_CDF <- function(t) {
  pgamma(t, shape = BIBOOSTER_shape, rate = BIBOOSTER_rate)
}

#compare fitted models with data
comparison_data = data_field_study_cuttoffs%>%filter(Treatment %in% c("CEVA","BI","BI+booster") & time_point<=118)%>%
  mutate(CO_5_predict = case_when(
    Treatment == "CEVA" ~ CEVA_CDF(time_point),
    Treatment == "BI" ~ BI_CDF(time_point),
    Treatment == "BI+booster" ~ BIBOOSTER_CDF(time_point),
    NULL))

ggplot(comparison_data) +
  geom_point( aes(x = time_point, y = CO_5,colour= Farm) )+
  geom_line( aes(x = time_point, y = CO_5_predict),colour= "black", linetype = "dashed") +
  labs(title = "Percentage Below Cut-off",
       x = "Time Point",
       y = "Percentage Below Cut-off ",
       color = "Parameter") +
  theme_minimal()+
  facet_wrap(. ~ Group + Treatment) 


#create a table with mean and variance for each of the treatments
summary_table <- data.frame(
  Group = c("CEVA", "BI", "BIBOOSTER"),
  Mean = c(CEVA_shape/CEVA_rate, BI_shape/BI_rate, BIBOOSTER_shape/BIBOOSTER_rate),
  Variance = c(CEVA_shape/(CEVA_rate^2), BI_shape/(BI_rate^2), BIBOOSTER_shape/(BIBOOSTER_rate^2))
)

print(summary_table)

df<-data.frame(cbind(c(0:200),t(sapply(X = c(0:200), function(x){pgamma(x,summary_table$Mean^ 2/summary_table$Variance, scale = summary_table$Variance/summary_table$Mean)}))))
colnames(df) <- c("time_point", "CEVA", "BI", "BIBOOSTER")
ggplot(df, aes(x = time_point)) +
  geom_line(aes(y = CEVA, color = "CEVA")) +
  geom_line(aes(y = BI, color = "BI")) +
  geom_line(aes(y = BIBOOSTER, color = "BIBOOSTER")) +
  labs(title = "Gamma CDF for Treatments",
       x = "Time Point",
       y = "Cumulative Probability",
       color = "Treatment") +
  theme_minimal() +
  scale_color_manual(values = c("CEVA" = "blue", "BI" = "red", "BIBOOSTER" = "green"))

# Sensitivity analyses with cutt-off values 6 and 7 ####
# Cutt-off 6 ####
surv_data_field_study_co6 <- data_field_study_numeric %>%
  mutate(positive = as.numeric(Uitslag>=6)) 

surv_data_field_study_co6%>%filter(Treatment == "CEVA" & time_point<= 120)%>%create_surv_data() -> input_data_CEVA_co6
surv_data_field_study_co6%>%filter(Treatment == "BI" & time_point<= 120&!is.na(positive))%>%create_surv_data() -> input_data_BI_co6
surv_data_field_study_co6%>%filter(Treatment == "BI+booster" & time_point<= 120)%>%create_surv_data() -> input_data_BIBOOSTER_co6


# Fit the model for BI
fit_CEVA_co6 <- stan(
  file = "gamma_model.stan",
  data = input_data_CEVA_co6,
  iter = 2000,
  chains = 4
)

# Fit the model for BI
fit_BI_co6 <- stan(
  file = "gamma_model.stan",
  data = input_data_BI_co6,
  iter = 2000,
  chains = 4
)

# Fit the model for BI
fit_BIBOOSTER_co6 <- stan(
  file = "gamma_model.stan",
  data = input_data_BIBOOSTER_co6,
  iter = 2000,
  chains = 4
)

# functions with estimates

CEVA_shape_co6 <- median(extract(fit_CEVA_co6, "shape")[[1]])
CEVA_rate_co6 <- median(extract(fit_CEVA_co6, "rate")[[1]])
CEVA_CDF_co6 <- function(t) {
  pgamma(t, shape = CEVA_shape_co6, rate = CEVA_rate_co6)
}


BI_shape_co6 <- median(extract(fit_BI_co6, "shape")[[1]])
BI_rate_co6 <- median(extract(fit_BI_co6, "rate")[[1]])
BI_CDF_co6 <- function(t) {
  pgamma(t, shape = BI_shape_co6, rate = BI_rate_co6)
}


BIBOOSTER_shape_co6 <- median(extract(fit_BIBOOSTER_co6, "shape")[[1]])
BIBOOSTER_rate_co6 <- median(extract(fit_BIBOOSTER_co6, "rate")[[1]])
BIBOOSTER_CDF_co6 <- function(t) {
  pgamma(t, shape = BIBOOSTER_shape_co6, rate = BIBOOSTER_rate_co6)
}


#compare fitted models with data
comparison_data_co6 = data_field_study_cuttoffs%>%filter(Treatment %in% c("CEVA","BI","BI+booster") & time_point<=118)%>%
  mutate(CO_6_predict = case_when(
    Treatment == "CEVA" ~ CEVA_CDF_co6(time_point),
    Treatment == "BI" ~ BI_CDF_co6(time_point),
    Treatment == "BI+booster" ~ BIBOOSTER_CDF_co6(time_point),
    NULL))

ggplot(comparison_data_co6) +
  geom_point( aes(x = time_point, y = CO_6,colour= Farm) )+
  geom_line( aes(x = time_point, y = CO_6_predict),colour= "black", linetype = "dashed") +
  labs(title = "Percentage Below Cut-off",
       x = "Time Point",
       y = "Percentage Below Cut-off ",
       color = "Parameter") +
  theme_minimal()+
  facet_wrap(. ~ Group + Treatment) 


#create a table with mean and variance for each of the treatments
summary_table_co6 <- data.frame(
  Group = c("CEVA", "BI", "BIBOOSTER"),
  Mean = c(CEVA_shape_co6/CEVA_rate_co6, BI_shape_co6/BI_rate_co6, BIBOOSTER_shape_co6/BIBOOSTER_rate_co6),
  Variance = c(CEVA_shape_co6/(CEVA_rate_co6^2), BI_shape_co6/(BI_rate_co6^2), BIBOOSTER_shape_co6/(BIBOOSTER_rate_co6^2))
)

print(summary_table_co6)

df<-data.frame(cbind(c(0:200),t(sapply(X = c(0:200), function(x){pgamma(x,summary_table_co6$Mean^ 2/summary_table_co6$Variance, scale = summary_table_co6$Variance/summary_table_co6$Mean)}))))
colnames(df) <- c("time_point", "CEVA", "BI", "BIBOOSTER")
ggplot(df, aes(x = time_point)) +
  geom_line(aes(y = CEVA, color = "CEVA")) +
  geom_line(aes(y = BI, color = "BI")) +
  geom_line(aes(y = BIBOOSTER, color = "BIBOOSTER")) +
  labs(title = "Gamma CDF for Treatments",
       x = "Time Point",
       y = "Cumulative Probability",
       color = "Treatment") +
  theme_minimal() +
  scale_color_manual(values = c("CEVA" = "blue", "BI" = "red", "BIBOOSTER" = "green"))


# HI titre >=7 ####

surv_data_field_study_co7 <- data_field_study_numeric %>%
  mutate(positive = as.numeric(Uitslag>=7)) 

surv_data_field_study_co7%>%filter(Treatment == "CEVA" & time_point<= 120)%>%create_surv_data() -> input_data_CEVA_co7
surv_data_field_study_co7%>%filter(Treatment == "BI" & time_point<= 120&!is.na(positive))%>%create_surv_data() -> input_data_BI_co7
surv_data_field_study_co7%>%filter(Treatment == "BI+booster" & time_point<= 120)%>%create_surv_data() -> input_data_BIBOOSTER_co7


# Fit the model for BI
fit_CEVA_co7 <- stan(
  file = "gamma_model.stan",
  data = input_data_CEVA_co7,
  iter = 2000,
  chains = 4
)

# Fit the model for BI
fit_BI_co7 <- stan(
  file = "gamma_model.stan",
  data = input_data_BI_co7,
  iter = 2000,
  chains = 4
)

# Fit the model for BI
fit_BIBOOSTER_co7 <- stan(
  file = "gamma_model.stan",
  data = input_data_BIBOOSTER_co7,
  iter = 2000,
  chains = 4
)


# functions with estimates
CEVA_shape_co7 <- median(extract(fit_CEVA_co7, "shape")[[1]])
CEVA_rate_co7 <- median(extract(fit_CEVA_co7, "rate")[[1]])
CEVA_CDF_co7 <- function(t) {
  pgamma(t, shape = CEVA_shape_co7, rate = CEVA_rate_co7)
}


BI_shape_co7 <- median(extract(fit_BI_co7, "shape")[[1]])
BI_rate_co7 <- median(extract(fit_BI_co7, "rate")[[1]])
BI_CDF_co7 <- function(t) {
  pgamma(t, shape = BI_shape_co7, rate = BI_rate_co7)
}


BIBOOSTER_shape_co7 <- median(extract(fit_BIBOOSTER_co7, "shape")[[1]])
BIBOOSTER_rate_co7 <- median(extract(fit_BIBOOSTER_co7, "rate")[[1]])
BIBOOSTER_CDF_co7 <- function(t) {
  pgamma(t, shape = BIBOOSTER_shape_co7, rate = BIBOOSTER_rate_co7)
}

#compare fitted models with data
comparison_data_co7 = data_field_study_cuttoffs%>%filter(Treatment %in% c("CEVA","BI","BI+booster") & time_point<=118)%>%
  mutate(CO_7_predict = case_when(
    Treatment == "CEVA" ~ CEVA_CDF_co7(time_point),
    Treatment == "BI" ~ BI_CDF_co7(time_point),
    Treatment == "BI+booster" ~ BIBOOSTER_CDF_co7(time_point),
    NULL))

ggplot(comparison_data_co7) +
  geom_point( aes(x = time_point, y = CO_7,colour= Farm) )+
  geom_line( aes(x = time_point, y = CO_7_predict),colour= "black", linetype = "dashed") +
  labs(title = "Percentage Below Cut-off",
       x = "Time Point",
       y = "Percentage Below Cut-off ",
       color = "Parameter") +
  theme_minimal()+
  facet_wrap(. ~ Group + Treatment) 


#create a table with mean and variance for each of the treatments
summary_table_co7 <- data.frame(
  Group = c("CEVA", "BI", "BIBOOSTER"),
  Mean = c(CEVA_shape_co7/CEVA_rate_co7, BI_shape_co7/BI_rate_co7, BIBOOSTER_shape_co7/BIBOOSTER_rate_co7),
  Variance = c(CEVA_shape_co7/(CEVA_rate_co7^2), BI_shape_co7/(BI_rate_co7^2), BIBOOSTER_shape_co7/(BIBOOSTER_rate_co7^2))
)

print(summary_table_co7)

df<-data.frame(cbind(c(0:200),t(sapply(X = c(0:200), function(x){pgamma(x,summary_table_co7$Mean^ 2/summary_table_co7$Variance, scale = summary_table_co7$Variance/summary_table_co7$Mean)}))))
colnames(df) <- c("time_point", "CEVA", "BI", "BIBOOSTER")
ggplot(df, aes(x = time_point)) +
  geom_line(aes(y = CEVA, color = "CEVA")) +
  geom_line(aes(y = BI, color = "BI")) +
  geom_line(aes(y = BIBOOSTER, color = "BIBOOSTER")) +
  labs(title = "Gamma CDF for Treatments",
       x = "Time Point",
       y = "Cumulative Probability",
       color = "Treatment") +
  theme_minimal() +
  scale_color_manual(values = c("CEVA" = "blue", "BI" = "red", "BIBOOSTER" = "green"))



