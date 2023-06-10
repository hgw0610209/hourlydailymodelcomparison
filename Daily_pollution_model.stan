
data {
  int<lower=0> N; // number of observations
  int<lower=0> K; // number of periodical covariates
  matrix[N, K] x; // periodical covariates matrix, a part in X_thp in the model
  real y[N]; // response in the model,y
  int<lower=0> Nw; //=7, as there are total 7 weekdays
  int<lower=1,upper=Nw> week[N]; // weekday indicator, a colum in X_thp
  int<lower=0> Nd;// the number of how many days, e.g. two year data usually is 365+365=730, could be 731 as well if 366
  int<lower=1,upper=Nd> day[N];// day indicator, a column in X_thp
  int<lower=1> N_p; // the number of pollutants
  
  int<lower=1> trend_reference_day; // the reference day for time trend
  
  int<lower=1> pollutant[N]; // pollutant indicator 
  
  int<lower=1> N_cov; // the number of covariate using for covariance regression
  matrix[Nd, N_cov] x_cov; // covariance predictor matrix
  
}


parameters {
  vector<lower=0,upper=1>[N_p] arcoef_phi; // the AR coef of daily random effect
 
  matrix[N_p, (Nw-1)] factor_w_pre;// week factor, which weekday it belongs, one colum in X_thp
  matrix[Nd, N_p] factor_d; // daily random effects. =D_tp
  vector[N_p] miu; // mu_p
  vector[N_p-1] factor_trend_raw; // time trend, which assumes PM25 and pm10 share the same trend
  matrix[N_p, K] betas; //the regression parameters
  vector<lower=0> [N_p] alpha_pre ; // shape parameter for gamma distribution
  // vector<lower=0, upper=1>[N_p] mixture_rho[N_rho]; // mixture proportion rho in the model
  
   cholesky_factor_corr[N_p] L1; //cholesky factor of covariance L in the model
  vector<lower=0>[N_p] tau;// tau in the model
   
   //for covariance
     vector[N_p] miu_cov; // b_j in the model, is the intercept
     vector[N_cov] betas_cov[N_p]; // eta in the model, is the regression coef for varying covariance.

}

 transformed parameters
 {
  // save the needed parameters
   matrix[(N_p-1),(N_p-1)] L_save;
   matrix[N_p, Nw] factor_w;
   vector[N_p] factor_trend; // time trend
   
   vector[N_p] alpha;
  // 
  matrix[N_p,N_p] L[Nd];
  matrix[N_p,N_p] L_prior1;
  matrix[N_p,N_p] L_prior;
  vector [N_p] temp_arcoef_phi;
  
  // use wedesday as baseline
  for(i in 1:2)
  {
    factor_w[,i]=factor_w_pre[,i];
  }
  factor_w[,3]=rep_vector(0,N_p);
    for(i in 4:Nw)
  {
    factor_w[,i]=factor_w_pre[,(i-1)];
  }
  
  
  for(i in 1:N_p)
  {
    alpha[i]=1/((alpha_pre[i])^2);
  }

   factor_trend[1]=factor_trend_raw[1];

  for(i in 1:(N_p-1))
  {
   factor_trend[i+1]=factor_trend_raw[i]; //which assumes PM25 and pm10 share the same trend
  }

  for(i in 1:Nd)
  {
  L[i]=L1;
     for(j in 1:N_p)
       {
          L[i][N_p,j]=miu_cov[j] + x_cov[i,] * betas_cov[j];
        }
  }
  
  L_prior1=quad_form_diag(tcrossprod(L[1]), tau);
  
  for(i in 1:5)
  {
    for(j in 1:N_p)
    {
    temp_arcoef_phi[j]=arcoef_phi[j]^i;
      
    }
   //temp_arcoef_phi[2]=arcoef_phi[2]^i;
    //temp_arcoef_phi[3]=arcoef_phi[3]^i; 
    //temp_arcoef_phi[4]=arcoef_phi[4]^i;
    L_prior1 += quad_form_diag(tcrossprod(L[1]), tau .* temp_arcoef_phi);
  }
  
  L_prior=(L_prior1+L_prior1')./2;
  
  
  // 
  L_save = L[1, :(N_p-1), :(N_p-1)];
  
}

model {
  
  alpha_pre ~ exponential(0.5);
  miu ~ normal(0,10);
  factor_trend_raw ~ normal(0,10);
  
  
  
    for(j in 1:N_p)
    {
     factor_w_pre[j,] ~ normal(0,10);
    }
  
  
  miu_cov ~ normal(0,1);
  for(i in 1:N_p)
  {
  betas[i] ~ normal(0,10);
  betas_cov[i] ~ normal(0,0.05);
    
  }
  // Priors
  tau ~ exponential(0.5);
  arcoef_phi ~ uniform(0,1);
  
  L1 ~ lkj_corr_cholesky(1);
  
  
//##
  factor_d[1]' ~ multi_normal(rep_vector(0,N_p), L_prior);
  
  for(i in 2:Nd)
  {
  factor_d[i]' ~ multi_normal_cholesky(diag_matrix(arcoef_phi) * factor_d[(i-1)]', diag_pre_multiply(tau, L[i]));
  }
  
  //########
  // Likelihood
 
             

for(n in 1:N)
  {
     target +=  gamma_lpdf(y[n] | alpha[pollutant[n]], (alpha[pollutant[n]] / exp(miu[pollutant[n]] + x[n,] * betas[pollutant[n]]' + 
                                 factor_trend[pollutant[n]]*((day[n]-trend_reference_day)/2091.5)+
                                         factor_w[pollutant[n],week[n]] + 
                                         factor_d[day[n],pollutant[n]]) ));
                                         
  }                                     

}





