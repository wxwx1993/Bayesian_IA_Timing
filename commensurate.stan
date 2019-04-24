// Bayesian analysis using Pocock's method with piecewise constant baseline hazard
data {
  int<lower=1> N;    // rows of data
  vector[N]  Y;       // data
  int<lower=1,upper=2>  k[N];
  int<lower=0,upper=1>  hist[N];
}
parameters {
  real mu01; 
  real mu02; 
  real<lower=0>  s2; 
  real<lower=0>  s3; 
  real mu1; 
  real mu2; 
  real<lower=0> ss;
}
transformed parameters {
  real mu[4];
  real s2_sqrt;
  real s3_sqrt;
  real ss_sqrt;
  real mudiff;
  mu[1] = mu1;  
  mu[2] = mu2;
  mu[3] = mu01;  
  mu[4] = mu02;
  s2_sqrt=sqrt(s2);
  s3_sqrt=sqrt(s3);
  ss_sqrt = sqrt(ss);
  mudiff= mu[2]-mu[1];
} 
model {  
  // Variance and precision of study-specific prior
  mu[3] ~ normal(0,100);
  mu[4] ~ normal(0,100);
  s2 ~ inv_gamma(0.02,1);
  s3 ~ inv_gamma(0.02,1);
  ss ~ inv_gamma(0.01,1);
  mu[1] ~ normal(mu[3],s2_sqrt);
  mu[2] ~ normal(mu[4],s3_sqrt);
   // Likelihoods 
    for (i in 1:N)
          // define the likelihood
          Y[i] ~ normal(mu[hist[i]*2+k[i]],ss_sqrt);
}