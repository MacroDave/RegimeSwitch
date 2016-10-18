data {
  int T; // The number of time periods
  int S; //number of states
  int L; // the number of AR terms
  vector[T] Y; // the time series
}

parameters {
  simplex[2] row_1;
  simplex[2] row_2;
  real<lower=0> constant_high;
  real<upper=0> constant_low;
  real<lower=0> sigma;
  vector[L] beta;
}

transformed parameters {
  matrix[2,2] your_markov_transition_matrix;
  your_markov_transition_matrix[1] = row_1';
  your_markov_transition_matrix[2] = row_2';
}

model {
  
  vector[2] p_state;
  // set initial value
  p_state = rep_vector(0.5, 2);

  // prior for the transition matrix
  row_1 ~ dirichlet(rep_vector(2.0, 2));
  row_2 ~ dirichlet(rep_vector(2.0, 2));

  beta ~ normal(0, 1);
  sigma ~ cauchy(0, 5);

  for(t in 5:T) {
    p_state = your_markov_transition_matrix * p_state;
    target += log_sum_exp(log(p_state[1]) + normal_lpdf(Y[t] | constant_high + beta[1]*Y[t-1] + beta[2]*Y[t-2] + beta[3]*Y[t-3] + beta[4]*Y[t-4], sigma),
                                          log(p_state[2]) + normal_lpdf(Y[t] | constant_low + beta[1]*Y[t-1] + beta[2]*Y[t-2] + beta[3]*Y[t-3] + beta[4]*Y[t-4], sigma));
  }
}
