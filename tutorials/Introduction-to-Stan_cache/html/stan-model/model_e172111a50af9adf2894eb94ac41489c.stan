functions{
  
  real depot_1cmt(real dose, real cl, real v, real ka, 
                  real time_since_dose){
    
    real ke = cl/v;
    
    real conc = dose/v * ka/(ka - ke) * 
              (exp(-ke*time_since_dose) - exp(-ka*time_since_dose));
    
    return conc;
    
  }
  
}
