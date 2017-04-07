data {
  int T; // The number of time periods
  int S; //number of states
  int L; // the number of AR terms
  vector[T] Y; // the time series
}
parameters {
  real<lower = 0, upper = 1> p11;
  real<lower = 0, upper = 1> p22;
  real<lower = 0> sigma; // standard deviations
  ordered[2] alpha;
  vector[L] beta;   //lagged AR terms in GNP equation
}
transformed parameters {
  vector[T] Xi;
  vector[T-4] l_likelihood;
  matrix[T-4, 4] contributions;
  Xi[1:4] = rep_vector(log(0.05), 4);
  
  for(t in 5:T) {
    
    contributions[t-4, 1] = log(p11) + Xi[t-1] + normal_lpdf(Y[t] | alpha[1] + beta[1]*Y[t-1] + beta[2]*Y[t-2] + beta[3]*Y[t-2] + beta[4]*Y[t-4], sigma);
    contributions[t-4, 2] = log(1-p11) + Xi[t-1] + normal_lpdf(Y[t] | alpha[2] + beta[1]*Y[t-1] + beta[2]*Y[t-2] + beta[3]*Y[t-2] + beta[4]*Y[t-4], sigma);
    contributions[t-4, 3] = log(1-p22) + log(1-exp(Xi[t-1])) + normal_lpdf(Y[t] | alpha[1] + beta[1]*Y[t-1] + beta[2]*Y[t-2] + beta[3]*Y[t-2] + beta[4]*Y[t-4], sigma);
    contributions[t-4, 4] = log(p22) + log(1-exp(Xi[t-1])) + normal_lpdf(Y[t] | alpha[2] + beta[1]*Y[t-1] + beta[2]*Y[t-2] + beta[3]*Y[t-2] + beta[4]*Y[t-4], sigma);
    l_likelihood[t-4] = log_sum_exp(contributions[t-4]) ;
    Xi[t] = log_sum_exp(contributions[t-4, 1], contributions[t-4, 3]) - l_likelihood[t-4];
  }
}
model {
  // priors 
  p11 ~ beta(5, 1);
  p22 ~ beta(5, 1);
  
  sigma ~ normal(0, 1); //variance term in GNP equations
  
  to_vector(beta) ~ normal(0, 1); //coefficients on AR terms in GNP equations  
  alpha[1] ~ normal(-1,1); //constant term in high growth state GNP equation
  alpha[2] ~ normal(1,1); //constant term in high growth state GNP equation
  
  // likelihood
  for(t in 5:T) {
    target += log_sum_exp(contributions[t-4])  ;
  }
}
