//Growth function
functions{
  //Growth function for use with Runge-Kutta method
  real growth(real y, real g_max, real S_max, real k){
    return g_max * exp(-0.5 * pow(log(y / S_max) / k, 2));
  }
}

// Data structure
data {
  int n_obs;
  real y_obs[n_obs];
  real delta_obs[n_obs];
  real interval[n_obs];
}

// The parameters accepted by the model.
parameters {
  //Species level
  real<lower=0> pop_max_growth;
  real<lower=0> pop_size_at_max_growth;
  real<lower=0> pop_k;

  //Global level
  real<lower=0> global_error_sigma;
}

// The model to be estimated.
model {
  real delta_hat[n_obs];

  for(i in 1:n_obs){
    delta_hat[i] = growth(y_obs[i], pop_max_growth,
            pop_size_at_max_growth,
            pop_k) * interval[i];
  }

  //Likelihood
  delta_obs ~ normal(delta_hat, global_error_sigma);

  //Priors
  //Species level
  pop_max_growth ~lognormal(0, 5);
  pop_size_at_max_growth ~lognormal(1,5);
  pop_k ~lognormal(0, 5);

  //Global level
  global_error_sigma ~cauchy(0, 5);
}


generated quantities{
}
