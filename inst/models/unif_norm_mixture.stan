// Simple mixture of a uniform [0,5] and N(15, 1)

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0, upper=1> theta;
}

model {
  sigma ~ normal(0, 2);
  mu ~ normal(15, 6);
  theta ~ beta(5, 5);
  for (n in 1:N)
   target += log_mix(theta,
                     uniform_lpdf(y[n] | 0, 5),
                     normal_lpdf(y[n] | mu, sigma));
}

