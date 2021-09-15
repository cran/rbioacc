//functions {
int findfirst(real t, vector xt) {
  int i = 0 ;
  if(t == max(xt)){
    i = num_elements(xt) - 1 ;
    return i ;
  } else if(t < min(xt) || t > max(xt)){
    return i ;
  } else {
    while (t >= xt[i+1]){
      i = i+1 ;
    }
    return i ;
  }
}

real interpolate(real x, vector xpt, vector ypt){
  if(x >= min(xpt) && x <= max(xpt)){
    int idx = findfirst(x, xpt) ;
    return ypt[idx] + (x - xpt[idx]) * (ypt[idx+1] - ypt[idx]) / (xpt[idx+1] - xpt[idx]) ;
  } else{
    return 0.0 ;
  }
}
//}

real[] odeTK(real t,      // time
             real[] y,    // variables
             real[] theta,
             real[] x_r,
             int[] x_i) {
  
  // parameters
  
  int lentp_rmNA = x_i[1] ;
  int lentp = x_i[2] ;
  int n_exp = x_i[3] ; 
  int n_out = x_i[4] ; 
  int n_met = x_i[5] ;

  real ku[n_exp] = theta[1:n_exp] ;
  real ke[n_out] = theta[(n_exp+1):(n_exp+n_out)] ;
  real km[n_met] = theta[(n_exp+n_out+1):(n_exp+n_out+n_met)] ;
  real kem[n_met] = theta[(n_exp+n_out+n_met+1):(n_exp+n_out+n_met+n_met)] ;
   
  // vector[1+n_met] dydt ;
  real dydt[1+n_met] ;
  
  real tacc = x_r[1] ;
  // real tp_rmNA[lentp_rmNA] = x_r[2:(lentp_rmNA+1)] ;
  vector[lentp_rmNA] tp_rmNA = to_vector(x_r[2:(lentp_rmNA+1)]) ;
  // WORK ONLY FOR ONE EXPOSURE PROFILE
  
  // real Cexp_rmNA[lentp_rmNA] = x_r[(lentp_rmNA+2):(lentp_rmNA+2+lentp_rmNA)] ;
  vector[lentp_rmNA] Cexp_rmNA = to_vector(x_r[(lentp_rmNA+2):(lentp_rmNA+1+lentp_rmNA)]) ;

  // latent vairable
  real U ;
  real M ;
  real E ;
  
  if(n_met == 0){
    M = 0 ;
  } else{
    M = sum(km) ;
  }
  E = sum(ke) ;
  
  // model
  // U = ku[1] * interpolate(t, tp_rmNA, Cexp_rmNA[1:lentp_rmNA, 1]) ;
  U = ku[1] * interpolate(t, tp_rmNA, Cexp_rmNA) ;
  // if(n_exp >= 2){
  //   for(i in 2:n_exp){
  //     U = U +  ku[i] * interpolate(t, tp_rmNA, Cexp_rmNA[1:lentp_rmNA, i]) ;
  //   }
  // }
  if(t <= tacc){
    // Accumulation
    dydt[1] = U - (E + M) * y[1] ;
  } else{
    // Depuration
    dydt[1] = - (E + M) * y[1] ;
  }
  if(n_met > 0){
    for(i in 2:(n_met+1)){
      dydt[i] = km[i] * y[1] -  kem[i] * y[i] ;
    }
  }
  return(dydt) ;
}


// matrix solve_TK(
//   real[] y0, real t0, real[] ts, real[] theta, int[] n_int,
//   real tacc, real[] tp_rmNA, real[] Cexp_rmNA, real[] odeParam){
//     
//     //vector[1+size(Cexp_rmNA) + size(tp_rmNA)] vector_r ;
//     real vector_r[1+size(Cexp_rmNA) + size(tp_rmNA)]  ;
//     // real vector_r[1+size(tp_rmNA)]  ;
//     
//     vector_r[1] = tacc ;
//     // vector_r[2:(size(Cexp_rmNA)+1)] = tp_rmNA ;
//     for(i in 1:size(tp_rmNA)){
//       vector_r[i+1] = tp_rmNA[i] ;
//     }
//     for(i in 1:size(Cexp_rmNA)){
//        vector_r[i+1+size(tp_rmNA)] = Cexp_rmNA[i] ;
//     }
//     
//     // vector_r[(size(Cexp_rmNA)+2):(size(Cexp_rmNA)+size(tp_rmNA)+1)] = Cexp_rmNA ;
//     
//     // append_row(to_vector(tp_rmNA), to_vector(Cexp_rmNA));
// 
//     return(to_matrix(
//       integrate_ode_rk45(odeTK, y0, t0, ts, theta,
//                          //to_array_1d(vector_r),
//                          vector_r,
//                          // to_array_1d(append_row(to_vector(vector_r), to_vector(Cexp_rmNA))),
//                          to_array_1d(n_int),
//                          // additional control parameters for the solver: real rel_tol, real abs_tol, int max_num_steps
//                          odeParam[1], odeParam[2], odeParam[3]))) ;
// 
//   }
