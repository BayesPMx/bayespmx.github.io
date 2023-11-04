transformed parameters{
  vector[n_obs] ipred;
  
  for(i in 1:n_obs){
    ipred[i] = depot_1cmt(dose, CL, V, KA, time_since_dose[i]);
  }
  
}
