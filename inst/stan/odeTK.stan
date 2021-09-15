functions{
  
#include /include/linear_interpolation.stan
// #include /include/TKodeSolver.stan
}

data {
  // Number of replicate
  int<lower=0> n_rep ;
  
  // Time points
  int<lower=0> lentp;
  real tp[lentp] ;
  
  // Exposure profiles
  int<lower=0> n_exp ;
  int<lower=0> lentp_rmNA ;
  vector[lentp_rmNA] tp_rmNA ;
  matrix[lentp_rmNA, n_exp] Cexp_rmNA ;
  
  // Internal concentraion
  // Growth
  int<lower=0> n_out ;
  real CGobs[lentp, n_out, n_rep];

  // Metabolites
  int<lower=0> n_met ;
  real Cmet[lentp, n_met, n_rep] ;

  //vector[1+n_met] y0;
  //real t0;

  real<lower=0> gmaxsup ;
  
  // TK accumulation / depuration
  int<lower=0> rankacc ;
  real<lower=0> tacc ;
  real<lower=0> C0 ;
  
  real unifMax ;
  int<lower=0> len_vt;
  vector[len_vt] vt ;
  
  // Parameters for integration of differentiol equations
  real y0[1+n_met];
  real t0 ;
  
  // real rel_tol;
  // real abs_tol;
  // int max_num_steps;
}
transformed data{
  real x_r[1+lentp_rmNA+lentp_rmNA] ;
  int x_int[5] ;
  
  x_r[1] = tacc ;
  
  for(i in 1:lentp_rmNA){
    x_r[1+i] = tp_rmNA[i] ;
  }
  for(i in 1:lentp_rmNA){
    x_r[i+lentp_rmNA+1] = Cexp_rmNA[i,1] ;
  }
  
  x_int[1] = lentp_rmNA ;
  x_int[2] = lentp ;
  x_int[3] = n_exp ;
  x_int[4] = n_out ;
  x_int[5] = n_met ;

}
parameters {
  vector[n_exp] log10ku ; // uptake
  vector[n_out] log10ke ;
  vector[n_met] log10km ;
  vector[n_met] log10kem ;
  
  real<lower=0> sigmaCGpred[n_out] ; 
  vector<lower=0>[n_met] sigmaCmetpred ;
  
  real<lower=0> gmax[n_out - 1] ;
  real<lower=0> G0[n_out -1];
  
}
transformed parameters{
  // PRIORS
  vector<lower=0>[n_exp] ku ;
  vector<lower=0>[n_out] ke ;
  vector<lower=0>[n_met] km ;
  vector<lower=0>[n_met] kem ;
  // with time-variable exposure profile,
  // exposure may have NA values
  matrix[lentp, n_exp] Cexp ;
  // little hack merging Cpred and Gpred
  matrix[lentp,n_out] CGpred ; 
  matrix[lentp,n_met] Cmetpred ;

  // int n_val[5] ;
  // real theta ;

  real y_sim[lentp, 1+n_met] ; 
  // matrix[lentp, 1+n_met] y_sim ;
  
  real theta[n_exp+n_out+n_met+n_met] ;
  
  for(i in 1:n_exp){
    ku[i] = 10 ^ log10ku[i] ;
  }
  for(i in 1:n_out){
    ke[i] = 10 ^ log10ke[i] ;
  }
  for(i in 1:n_met){
    km[i] = 10 ^ log10km[i] ;
    kem[i] = 10 ^ log10kem[i] ;
  }
  
  // n_val[1] = lentp_rmNA ;
  // n_val[2] = lentp ;
  // n_val[3] = n_exp ;
  // n_val[4] = n_out ;
  // n_val[5] = n_met ;
 
  // y_sim = ode_rk45(odeTK, y0, t0, tp,lentp_rmNA,entp, n_exp,
  // n_out,n_met,n_val, tp_rmNA, Cexp_rmNA, tacc, ku, ke, km, kem) ;
  
  for(i in 1:n_exp){
    theta[i] = ku[i] ;
  }
   for(i in 1:n_out){
    theta[n_exp+i] = ke[i] ;
  }
   for(i in 1:n_met){
    theta[n_exp+n_out+i] = km[i] ;
  }
   for(i in 1:n_met){
    theta[n_exp+n_out+n_met+i] = kem[i] ;
  }

  y_sim = integrate_ode_rk45(odeTK, y0, t0, tp, theta, x_r, x_int) ;
  // y_sim = run(y0, t0, tp,lentp_rmNA,entp, n_exp,
  // n_out,n_met,n_val, tp_rmNA, Cexp_rmNA, tacc, ku, ke, km, kem) ;

  for(t in 1:lentp){
    CGpred[t,1] = y_sim[t, 1] ;
    for(i in 1:n_met){
      Cmetpred[t,i] = y_sim[t, i+1] ;
    }
  }
  // GROWTH
  if(n_out == 2){
    for(t in 1:lentp){
      CGpred[t, 2] = (G0[1] - gmax[1]) * exp(-ke[2] * tp[t]) + gmax[1] ;
    }
  }
}
model {
  // PRIORS
  target += uniform_lpdf(log10ku | -5, 5) ;
  target += uniform_lpdf(log10ke | -5, 5) ;
  target += uniform_lpdf(log10km | -5, 5) ;
  target += uniform_lpdf(log10kem | -5, 5) ;
  target += uniform_lpdf(sigmaCGpred[1] | 0, unifMax) ;
  target += uniform_lpdf(sigmaCmetpred | 0, unifMax) ;
  if(n_out == 2){
     target +=  uniform_lpdf(sigmaCGpred[2] | 0, unifMax) ;
     target +=  uniform_lpdf(gmax[1] | gmaxsup/6, gmaxsup) ;
     target +=  uniform_lpdf(G0[1] | 0, gmaxsup) ;
  }
  
  for(rep in 1:n_rep){
    // ACCUMULATION PHASE (0 =< t =< tacc) #
    for(t in 1:rankacc){
      // Parent compound
      if(!is_inf(CGobs[t,1,rep])){
        target += normal_lpdf(CGobs[t,1,rep] | CGpred[t, 1], sigmaCGpred[1]) ;
      }
      // Metabolites
      for(i in 1:n_met){
        if(!is_inf(Cmet[t,i,rep])){
          target += normal_lpdf(Cmet[t,i, rep] | Cmetpred[t,i], sigmaCmetpred[i]) ;
        }
      }
    }
    //DEPURATION PHASE (t > tacc)
    for(t in (rankacc+1):lentp){
      // Parent compound
      if(!is_inf(CGobs[t,1,rep])){
        target += normal_lpdf(CGobs[t,1,rep] | CGpred[t, 1], sigmaCGpred[1]) ;
      }
      // Metabolites
      for(i in 1:n_met){
        if(!is_inf(Cmet[t,i,rep])){
          target += normal_lpdf(Cmet[t,i,rep] | Cmetpred[t,i], sigmaCmetpred[i]) ;
        }
      }
    }
    // GROWTH
    if(n_out == 2){
      for(t in 1:lentp){
        if(!is_inf(CGobs[t,2,rep])){
          target += normal_lpdf(CGobs[t,2,rep] | CGpred[t,2], sigmaCGpred[2]) ;
        }
      }
    }
  }
}
