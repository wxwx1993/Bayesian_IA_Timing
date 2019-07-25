// Bayesian analysis in "Optimizing Interim Analysis Timing for Bayesian Adaptive Commensurate Designs"
data {
  int<lower=1> N;    // rows of data
  vector[N]  Y;       // data
  int<lower=1,upper=2>  k[N];
}
parameters {
  real mu01; 
  real mu02; 
  real<lower=0>  s2; 
  real mu1; 
  real mu2; 
  real<lower=0> ss;
}
transformed parameters {
  real mu[2];
  real s_sqrt;
  real ss_sqrt;
  real mudiff;
  mu[1] = mu1;  
  mu[2] = mu2;
  s_sqrt = sqrt(s2);
  ss_sqrt = sqrt(ss);
  mudiff= mu[2]-mu[1];
} 
model {  
  //  Variance and variance of study-specific prior
  mu01 ~ normal(0,100);
  mu02 ~ normal(0,100);
  s2 ~ inv_gamma(0.01,1);
  mu[1] ~ normal(mu01,s_sqrt);
  mu[2] ~ normal(mu02,s_sqrt);
  ss ~ inv_gamma(0.01,1);
   //  Likelihoods 
  for (i in 1:N)
    Y[i] ~ normal(mu[k[i]],ss_sqrt);
}
