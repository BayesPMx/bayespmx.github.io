data{
  
  int n_obs;                    // Number of observations
  real<lower = 0> dose;         // Dose amount
  array[n_obs] real time;       // Times at which we have observations
  real time_of_first_dose;      // Time of first dose
  vector[n_obs] dv;             // Observed PK data
  
  real<lower = 0> scale_cl;     // Prior Scale parameter for CL
  real<lower = 0> scale_v;      // Prior Scale parameter for V
  real<lower = 0> scale_ka;     // Prior Scale parameter for KA
  
  real<lower = 0> scale_sigma;  // Prior Scale parameter for lognormal error
  
  int n_pred;                   // Number of new times at which to make a prediction
  array[n_pred] real time_pred; // New times at which to make a prediction
 
}
