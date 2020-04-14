#include <Rcpp.h>
#include <math.h> 
#include <R.h>
#include <Rmath.h> 
#include <stdlib.h> 
#include <stdio.h>

using namespace Rcpp;

double awert(double x){

   /* This function provides an absolute value of a double x */
 
   if(x < 0)
     return(-x);
   else
    return(x);
}


double epanc(double x){

   /* This function calculates the value of the Epanechnikov kernel at point x */

   if(x > 1 || x < -1)
     return(0.0);
   else
     return (3.0 / 4.0) * (1 - x * x);
}


double s_t_0(double u, double h, int T){
 
   /* This function computes the kernel constant S_T0 needed for the local linear weights 
      of the kernel average hat(psi)_T(u,h) */ 

   int t, t_min, t_max;
   double result;

   result = 0.0;

   t_min = (int)(floor((u - h) * (float)T));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(ceil((u + h) * (float)T));      
   if(t_max > T) t_max = T;

   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)T - u) / h);
   }
   return(result / (T * h));
}


double s_t_1(double u, double h, int T){
 
   /* This function computes the kernel constant S_T1 needed for the local linear weights 
      of the kernel average hat(psi)_T(u,h) */ 

   int t, t_min, t_max;
   double result;

   result = 0.0;

   t_min = (int)(floor((u - h) * (float)T));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(ceil((u + h) * (float)T));      
   if(t_max > T) t_max = T;

   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)T - u) / h) * ((t / (float)T - u) / h);
   }
   return(result / (T * h));
}

double s_t_2(double u, double h, int T){
   
   /* This function computes the kernel constant S_T2 needed for the local linear weights 
    of the kernel average hat(psi)_T(u,h) */ 
   
   int t, t_min, t_max;
   double result;
   
   result = 0.0;
   
   t_min = (int)(floor((u - h) * (float)T));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(ceil((u + h) * (float)T));      
   if(t_max > T) t_max = T;
   
   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)T - u) / h) * ((t / (float)T - u) / h) * ((t / (float)T - u) / h);
   }
   return(result / (T * h));
}

NumericVector kernel_averages(int T, NumericVector gset, NumericVector correct, NumericVector data, double sigmahat, int N){
   /* Inputs: 
      T        	length of time series
      gset        vector of location-bandwidth points (u_1,...,u_N,h_1,...h_N)
      correct     vector of corrections of lenghth N (lambda(h_1), ..., lambda(h_N))
      data        time series of length T (y_1, ..., y_T)
      N           number of location-bandwidth points in gset
       
      Output:
      result      vector of normalized kernel averages without corrections 
                  (\phi(u_1,h_1)/sigmahat,...,\phi(u_N,h_N)/sigmahat) and vector of absolute values with
                  corrections  (abs(\phi(u_1,h_1)/sigmahat) - \lambda(h_1),...,abs(\phi(u_N,h_N)/sigmahat) - \lambda(h_N))
                  stacked together
   */

   int i, t, t_min, t_max;
   double u, h, x, temp, result_temp, weight_norm;

   NumericVector result(2 * N);

   for(i = 0; i < N; i++){
      u = gset[i];
      h = gset[i+N];
      t_min = (int)(floor((u - h) * (float)T));      
      if(t_min < 1) t_min = 1;
      t_max = (int)(ceil((u + h) * (float)T));      
      if(t_max > T) t_max = T;

      result_temp = 0;
      weight_norm = 0;
      
      if(u > h && u < 1.0-h){
        for(t = t_min; t < (t_max+1); t++){
           x = ((t / (float)T - u) / h);           
           temp = epanc(x) * x;
           result_temp = result_temp + temp * data[t-1];
           weight_norm = weight_norm + temp * temp;
	      }
      } else {
        for(t = t_min; t < (t_max+1); t++){
           x = ((t / (float)T - u) / h);
           temp = epanc(x) * (s_t_2(u, h, T) - s_t_1(u, h, T) * x);
           result_temp = result_temp + temp * data[t-1];
           weight_norm = weight_norm + temp * temp;
         }
      }
      result[i] = result_temp / (sqrt(weight_norm) * sigmahat);
      result[i + N] = awert(result[i]) - correct[i];
   }
  return (result);
}

// [[Rcpp::export]]
NumericVector gaussian_statistics(int T, int N_ts, int SimRuns, Rcpp::NumericVector gset, int N,
                                        Rcpp::NumericVector sigma_vec){
  /* Inputs: 
    T        	 length of time series
    N_ts        number of time series
    SimRuns     number of simulations to produce the quantiles
    gset        vector of location-bandwidth points (u_1,...,u_N,h_1,...h_N)
    correct     vector of corrections of lenghth N (lambda(h_1), ..., lambda(h_N))
    data        time series of length T (y_1, ..., y_T)
    N           number of location-bandwidth points in gset
  
  Output:
    Phi_vec     vector of length SimRuns tht consist of the calculated Gaussian statistic for each simulation run
  */
   Rcpp::NumericVector Phi_vec(SimRuns);
   Rcpp::NumericVector correct(N);
   Rcpp::NumericMatrix Phi(N_ts, N_ts);
   NumericVector vals(2 * N);
   NumericVector vals_cor(N);
   NumericVector data(T);
   int pos, i, j, k;
   double sigmahat;
  
   for (k = 0; k < N; k++){
      correct(k) = sqrt(2 * log(1 / (2 * gset[k+N])));
   }
   
   if (N_ts == 1){
      for(pos = 0; pos < SimRuns; pos++){
         NumericVector Z(T);
         Z = rnorm(T);
         data = sigma_vec(0) * Z;
         sigmahat = sigma_vec(0);
         vals = kernel_averages(T, gset, correct, data, sigmahat, N);
         vals_cor = vals[Rcpp::Range(N, 2 * N - 1)];
         Phi_vec(pos) = max(vals_cor);
      }
   } else {
      for(pos = 0; pos < SimRuns; pos++){
         NumericMatrix Z(T, N_ts);
         for (i = 0; i < N_ts; i++){
            Z(_, i) = rnorm(T);
         }
         for (i = 0; i < N_ts - 1; i++){
            for (j = i + 1; j < N_ts; j++){
               data = sigma_vec(i) * (Z(_, i) - mean(Z(_, i)))- sigma_vec(j) * (Z(_, j) - mean(Z(_, j)));
               sigmahat = sqrt(sigma_vec(i) * sigma_vec(i) + sigma_vec(j) * sigma_vec(j));
               vals = kernel_averages(T, gset, correct, data, sigmahat, N);
               vals_cor = vals[Rcpp::Range(N, 2 * N - 1)];
               Phi(i, j) = max(vals_cor);
            }
         }
         Phi_vec[pos] = max(Phi);
      }
   }
   return(Phi_vec);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix multiscale_statistics_multiple(int T, int N_ts, Rcpp::NumericMatrix data, Rcpp::NumericVector gset,
                              int N, Rcpp::NumericVector sigma_vec){

   /* Inputs: 
    T           length of time series
    N_ts        number of time series
    data        matrix of N_ts time series of length T
    gset        vector of location-bandwidth points (u_1,...,u_N,h_1,...h_N)
    N           number of location-bandwidth points in gset
    sigma_vec   vector of estimated sqrt(long-run variance) with an entry for each time series
    
    Output:
    Psi         matrix of pairwise multiscale statistics Psi_ij
    */
   Rcpp::NumericVector correct(N);
   Rcpp::NumericMatrix Psi(N_ts, N_ts);
   NumericVector vals(2 * N);
   NumericVector vals_cor(N);
   NumericVector time_series(T);
   int i, j, k;
   double sigmahat;
   
   for (k = 0; k < N; k++){
      correct[k] = sqrt(2 * log(1 / (2 * gset[k+N])));
   }
   
   for (i = 0; i < N_ts - 1; i++){
      for (j = i + 1; j < N_ts; j++){
         time_series = data(_, i) - data(_, j);
         sigmahat = sqrt(sigma_vec(i) * sigma_vec(i) + sigma_vec(j) * sigma_vec(j));
         vals = kernel_averages(T, gset, correct, time_series, sigmahat, N);
         vals_cor = vals[Rcpp::Range(N, 2 * N - 1)];
         Psi(i, j) = max(vals_cor);
      }
   }
   return(Psi);
}

// [[Rcpp::export]]
Rcpp::List multiscale_statistics(int T, Rcpp::NumericVector data, Rcpp::NumericVector gset,
                                 int N, double sigma){
   
   /* Inputs: 
    T        	 length of time series
    gset        vector of location-bandwidth points (u_1,...,u_N,h_1,...h_N)
    correct     vector of corrections of lenghth N (lambda(h_1), ..., lambda(h_N))
    data        time series of length T
    N           number of location-bandwidth points in gset
    
    Output:
    list with the following elements
    vals     vector of kernel averages (sign included)
    vals_cor vector of absolute value of normalized kernel averages with correction 
                (abs(\phi(u_1,h_1)/sigmahat) - \lambda(h_1),...,abs(\phi(u_N,h_N)/sigmahat) - \lambda(h_N))
    stat     multiscale statistics calculated as max(vals_cor)
    */
   Rcpp::NumericVector correct(N);
   Rcpp::NumericVector vals(2 * N);
   Rcpp::NumericVector vals_cor(N);
   int k;
   double Psi;

   for (k = 0; k < N; k++){
      correct[k] = sqrt(2 * log(1 / (2 * gset[k+N])));
   }
   
   vals = kernel_averages(T, gset, correct, data, sigma, N);
   vals_cor = vals[Rcpp::Range(N, 2 * N - 1)];
   Psi = max(vals_cor);
   return Rcpp::List::create(Rcpp::Named("vals") = vals,
                             Rcpp::Named("vals_cor") = vals_cor,
                             Rcpp::Named("stat") = Psi);
}