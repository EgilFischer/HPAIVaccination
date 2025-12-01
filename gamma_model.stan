
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

