//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] y;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real mu;
  real<lower=0> sigma;
  real<lower=0, upper=1> theta;
}

model {
  sigma ~ normal(0, 2);
  mu ~ normal(15, 6);
  //theta ~ beta(5, 5);
  theta ~ uniform(0, 1);
  for (n in 1:N)
   target += log_mix(theta,
                     uniform_lpdf(y[n] | 0, 5),
                     normal_lpdf(y[n] | mu, sigma));
}

