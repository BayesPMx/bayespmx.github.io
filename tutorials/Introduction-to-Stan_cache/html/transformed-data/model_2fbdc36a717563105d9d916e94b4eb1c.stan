transformed data{ 
  
  vector[n_obs] time_since_dose = to_vector(time) - time_of_first_dose;
  vector[n_pred] time_since_dose_pred = to_vector(time_pred) - 
                                        time_of_first_dose;
  int n_cmt = 2;                                        
  
}
