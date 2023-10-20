data{

  int<lower = 0> x; // observed positive responses
  int<lower = x> n; // number of responses
  
  real<lower = 0> alpha; // Value of alpha for the prior distribution
  real<lower = 0> beta;  // Value of beta for the prior distribution

  int n_new;               // length of new n values you want to simulate for 
  array[n_new] int n_star; // Number of future respondents for posterior predictions
  
}
parameters{

  real<lower = 0, upper = 1> theta;

}
model{
  // Priors
  theta ~ beta(alpha, beta);
  
  // Likelihood
  x ~ binomial(n, theta);
}
generated quantities{
  
  array[n_new] int x_star = binomial_rng(n_star, theta); 
  
}

