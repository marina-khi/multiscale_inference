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


double s_t_0(double u, double h, int t_len){
 
   /* This function computes the kernel constant S_T0 needed for the local linear weights 
      of the kernel average hat(psi)_T(u,h) */ 

   int t, t_min, t_max;
   double result;

   result = 0.0;

   t_min = (int)(floor((u - h) * (float)t_len));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(ceil((u + h) * (float)t_len));      
   if(t_max > t_len) t_max = t_len;

   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)t_len - u) / h);
   }
   return(result / (t_len * h));
}


double s_t_1(double u, double h, int t_len){
 
   /* t_lenhis function computes the kernel constant S_T1 needed for the local linear weights 
      of the kernel average hat(psi)_T(u,h) */ 

   int t, t_min, t_max;
   double result;

   result = 0.0;

   t_min = (int)(floor((u - h) * (float)t_len));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(ceil((u + h) * (float)t_len));      
   if(t_max > t_len) t_max = t_len;

   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)t_len - u) / h) * ((t / (float)t_len - u) / h);
   }
   return(result / (t_len * h));
}

double s_t_2(double u, double h, int t_len){
   
   /* This function computes the kernel constant S_T2 needed for the local linear weights 
    of the kernel average hat(psi)_T(u,h) */ 
   
   int t, t_min, t_max;
   double result;
   
   result = 0.0;
   
   t_min = (int)(floor((u - h) * (float)t_len));      
   if(t_min < 1) t_min = 1;
   t_max = (int)(ceil((u + h) * (float)t_len));      
   if(t_max > t_len) t_max = t_len;
   
   for(t = t_min; t < (t_max+1); t++){
      result += epanc((t / (float)t_len - u) / h) * ((t / (float)t_len - u) / h) * ((t / (float)t_len - u) / h);
   }
   return(result / (t_len * h));
}

NumericVector kernel_averages(int t_len, NumericVector gset, NumericVector correct,
                              NumericVector data, double sigmahat, int n,
                              int deriv_order = 0){
   /* Inputs: 
      t_len       Integer, length of time series
      gset        vector of location-bandwidth points (u_1,...,u_n,h_1,...h_n)
      correct     vector of corrections of lenghth n (lambda(h_1), ..., lambda(h_n))
      data        time series of length t_len (y_1, ..., y_t_len)
      sigmahat    value by which we normalise our test statstics
      n           number of location-bandwidth points in gset
      deriv_order Order of the derivative of the trend that is investigated
    
           
      Output:
      result      vector of normalized kernel averages without corrections 
                  (\phi(u_1,h_1)/sigmahat,...,\phi(u_n,h_n)/sigmahat) and vector of absolute values with
                  corrections  (abs(\phi(u_1,h_1)/sigmahat) - \lambda(h_1),...,abs(\phi(u_n,h_n)/sigmahat) - \lambda(h_n))
                  stacked together
   */

   int i, t, t_min, t_max;
   double u, h, x, temp, result_temp, weight_norm;

   NumericVector result(2 * n);

   for(i = 0; i < n; i++){
      u = gset[i];
      h = gset[i + n];
      t_min = (int)(floor((u - h) * (float)t_len));      
      if(t_min < 1) t_min = 1;
      t_max = (int)(ceil((u + h) * (float)t_len));      
      if(t_max > t_len) t_max = t_len;

      result_temp = 0;
      weight_norm = 0;
      
      if(u > h && u < 1.0-h){
        for(t = t_min; t < (t_max+1); t++){
           x = ((t / (float)t_len - u) / h);           
           temp = epanc(x) * x;
           result_temp = result_temp + temp * data[t-1];
           weight_norm = weight_norm + temp * temp;
	      }
      } else {
        for(t = t_min; t < (t_max+1); t++){
           x = ((t / (float)t_len - u) / h);
           if(deriv_order == 1){
              temp = epanc(x) * (s_t_0(u, h, t_len) * x - s_t_1(u, h, t_len));
           } else if(deriv_order == 0){
              temp = epanc(x) * (s_t_2(u, h, t_len) - s_t_1(u, h, t_len) * x);
           } else {temp = 0;}
           result_temp = result_temp + temp * data[t-1];
           weight_norm = weight_norm + temp * temp;
         }
      }
      result[i] = result_temp / (sqrt(weight_norm) * sigmahat);
      result[i + n] = awert(result[i]) - correct[i];
   }
  return (result);
}

//' Simulates distribution of the gaussian statistics
//'
//' @param t_len           Integer. Length of time series for analysis.
//' @param n_ts        Integer. Number of time series.
//' @param sim_runs    Integer. Number of simulations needed to produce the quantiles.
//' @param gset        A double vector of location-bandwidth points (u_1,...,u_n,h_1,...h_n).
//' @param ijset       An integer vector of indices of the countries for the comparison
//'                    (i_1, i_2, ..., j_{n_comparisons}).
//' @param sigma_vec   A double vector of the sqrt of long-run variances.
//' @param deriv_order Integer. Order of the derivative of the trend that is investigated.
//'                    Default is 0 => we analyse whether the trend itself >< than 0.
//' @param epidem      Boolean. Equal to true in cases we are fitting epiemic model. Default is false.
//' 
//' @return Phi_vec A double vector of length sim_runs that consists of the calculated Gaussian statistic for each simulation run
// [[Rcpp::export]]
NumericVector simulate_gaussian(int t_len, int n_ts, int sim_runs,
                                Rcpp::NumericVector gset, Rcpp::IntegerVector ijset,
                                Rcpp::NumericVector sigma_vec,
                                int deriv_order = 0, bool epidem = false){
   int n = gset.length() / 2;
   Rcpp::NumericVector Phi_vec(sim_runs);
   Rcpp::NumericVector correct_a(n);
   Rcpp::NumericVector correct_b(n);
   Rcpp::NumericMatrix Phi(n_ts, n_ts);
   NumericVector vals(2 * n);
   NumericVector vals_cor(n);
   NumericVector data(t_len);
   int pos, i, j, k, t, l, n_comparisons;
   double sigmahat, u, h, len, result_temp;
   
   for (k = 0; k < n; k++){
      len = 2 * gset[k + n];
      /* Here we need the whole length of the interval! */
      correct_b[k] = sqrt(2 * log(1 / len));
      correct_a[k] = sqrt(log(exp(1) / len)) / log(log(exp(exp(1)) / len));
   }
   
   if (n_ts == 1){
      for (pos = 0; pos < sim_runs; pos++) {
         NumericVector Z(t_len);
         Z = rnorm(t_len);
         data = sigma_vec(0) * Z;
         sigmahat = sigma_vec(0);
         vals = kernel_averages(t_len, gset, correct_b, data, sigmahat, n, deriv_order);
         vals_cor = vals[Rcpp::Range(n, 2 * n - 1)];
         Phi_vec(pos) = max(vals_cor);
      }
   } else if (epidem) {
      
      n_comparisons = ijset.length() / 2;

      for (pos = 0; pos < sim_runs; pos++) {
         NumericMatrix Z(t_len, n_ts);
         for (i = 0; i < n_ts; i++) {
            Z(_, i) = rnorm(t_len);
         }
         for (l = 0; l < n_comparisons; l++) {
            i = ijset[l] - 1;
            j = ijset[l + n_comparisons] - 1;
            NumericVector result(n);
            for(k = 0; k < n; k++){
               u = gset[k];
               h = gset[k + n];

               result_temp = 0;
               for(t = 1; t < (t_len + 1); t++){
                  if (t / (float)t_len >= (u - h) && t / (float)t_len <= (u + h)){
                     result_temp += Z(t - 1, i) - Z(t - 1, j);
                  }
               }
               result[k] = correct_a[k] * (awert(result_temp) / (sqrt(2 * t_len * 2 * h)) - correct_b[k]);
            }
            Phi(i, j) = max(result);
         }
         Phi_vec[pos] = max(Phi);
      }
   } else {
      n_comparisons = ijset.length() / 2;
      
      for(pos = 0; pos < sim_runs; pos++) {
         NumericMatrix Z(t_len, n_ts);
         for (i = 0; i < n_ts; i++) {
            Z(_, i) = rnorm(t_len);
         }
         for (l = 0; l < n_comparisons; l++) {
            i = ijset[l] - 1;
            j = ijset[l + n_comparisons] - 1;
            data = sigma_vec[i] * (Z(_, i) - mean(Z(_, i)))- sigma_vec[j] * (Z(_, j) - mean(Z(_, j)));
            sigmahat = sqrt(sigma_vec[i] * sigma_vec[i] + sigma_vec[j] * sigma_vec[j]); 
            vals = kernel_averages(t_len, gset, correct_b, data, sigmahat, n);
            vals_cor = vals[Rcpp::Range(n, 2 * n - 1)];
            Phi(i, j) = max(vals_cor);
         }
         Phi_vec[pos] = max(Phi);
      }
   }
   return(Phi_vec);
}

//' Calculates the statistics (pairwise and overall) in case of multiple time series
//'
//' @param t_len       Integer. Length of time series for analysis.
//' @param n_ts        Integer. Number of time series.
//' @param data        Double matrix t_len x n_ts. Each column consists of one time series.
//' @param gset        A double vector of location-bandwidth points (u_1,..., u_n, h_1, ...h_n).
//' @param ijset       An integer vector of indices of the countries for the comparison
//'                    (i_1, i_2, ..., j_{n_comparisons}).
//' @param sigma_vec   A double vector of length n_ts of the sqrt of long-run variances.
//' @param epidem      Boolean. Equal to true in cases we are fitting epiemic model. Default is false.
//' 
//' @return stat       A double matrix of pairwise multiscale statistics Psi_ij.
//' @return vals_cor   Matrix with n rows and n_ts * (n_ts - 1) / 2 columns where each column
//'                    contains values of the normalised Kernel averages in order to be able to perform 
//'                    the test on every separate interval.
// [[Rcpp::export]]
Rcpp::List compute_multiple_statistics(int t_len, int n_ts, Rcpp::NumericMatrix data,
                                       Rcpp::NumericVector gset, Rcpp::IntegerVector ijset,
                                       Rcpp::NumericVector sigma_vec,
                                       bool epidem = false){
   int n = gset.length() / 2;
   Rcpp::NumericVector correct_b(n);
   Rcpp::NumericVector correct_a(n);
   Rcpp::NumericMatrix stat(n_ts, n_ts);
   NumericVector vals(2 * n);
   NumericVector vals_cor(n);
   NumericVector time_series(t_len);
   int i, j, k, l;
   int n_comparisons = ijset.length() / 2;
   Rcpp::NumericMatrix vals_cor_mat(n, n_comparisons);
   double sigmahat, len;

   for (k = 0; k < n; k++){
      len = 2 * gset[k + n];
      /* Here we need the whole length of the interval! */
      correct_b[k] = sqrt(2 * log(1 / len));
      correct_a[k] = sqrt(log(exp(1) / len)) / log(log(exp(exp(1)) / len));
   }
   
   if (epidem){
      
      sigmahat = sqrt(mean(sigma_vec * sigma_vec));

      int t;
      double u, h, result_temp, weight_norm;
      
      for (l = 0; l < n_comparisons; l++) {
         i = ijset[l] - 1;
         j = ijset[l + n_comparisons] - 1;
         NumericVector result(n);
         for(k = 0; k < n; k++){
            u = gset[k];
            h = gset[k + n];
               
            result_temp = 0;
            weight_norm = 0;
               
            for(t = 1; t < (t_len + 1); t++){
               if (t / (float)t_len >= (u - h) && t / (float)t_len <= (u + h)){
                  result_temp += data(t-1, i) - data(t-1, j);
                  weight_norm += data(t-1, i) + data(t-1, j);
               }
               vals_cor[k] = correct_a[k] * (awert(result_temp) / (sqrt(weight_norm) * sigmahat) - correct_b[k]);
            }
         }
         stat(i, j) = max(vals_cor);
         vals_cor_mat(_, l) = vals_cor;
      }
   }  else {
      
      for (l = 0; l < n_comparisons; l++) {
         i = ijset[l] - 1;
         j = ijset[l + n_comparisons] - 1;
         time_series = data(_, i) - data(_, j);
         sigmahat = sqrt(sigma_vec(i) * sigma_vec(i) + sigma_vec(j) * sigma_vec(j));
         vals = kernel_averages(t_len, gset, correct_b, time_series, sigmahat, n);
         vals_cor = vals[Rcpp::Range(n, 2 * n - 1)];
         stat(i, j) = max(vals_cor);
         vals_cor_mat(_, l) = vals_cor;
      }
   }
   return(Rcpp::List::create(Rcpp::Named("vals_cor_matrix") = vals_cor_mat,
                             Rcpp::Named("stat") = stat));
}


//' Calculates statistics in case of one time series
//'
//' @param t_len           Integer. Length of time series for analysis.
//' @param data        Double vector of length t_len that consists the time series for analysis.
//' @param gset        A double vector of location-bandwidth points (u_1,...,u_n,h_1,...h_n).
//' @param sigma       Double. Equal to the sqrt of the long-run variance.
//' @param deriv_order Integer. Order of the derivative of the trend that is investigated.
//'                    Default is 0 => we analyse whether the trend itself >< than 0.
//' 
//' @return stat       A double matrix of pairwise multiscale statistics Psi_ij.
//' @return vals       A double vector of length n of kernel averages (sign included).
//' @return vals_cor   A double vector of absolute value of normalized kernel averages with
//'                    correction \eqn{(abs(\phi(u_1,h_1)/sigmahat) - \lambda(h_1), ...,
//'                    abs(\phi(u_n,h_n)/sigmahat) - \lambda(h_n))}
//' @return stat       Double. Our multiscale statistics calculated as max(vals_cor).
// [[Rcpp::export]]
Rcpp::List compute_statistics(int t_len, Rcpp::NumericVector data, Rcpp::NumericVector gset,
                              double sigma, int deriv_order) {
   int n = gset.length() / 2;
   Rcpp::NumericVector correct(n);
   Rcpp::NumericVector vals(2 * n);
   int k;
   double stat;

   for (k = 0; k < n; k++){
      correct[k] = sqrt(2 * log(1 / (2 * gset[k + n])));
   }
   
   vals = kernel_averages(t_len, gset, correct, data, sigma, n, deriv_order);
   stat  = max(vals[Rcpp::Range(n, 2 * n - 1)]);
   return Rcpp::List::create(Rcpp::Named("vals") = vals[Rcpp::Range(0, n - 1)],
                             Rcpp::Named("vals_cor") = vals[Rcpp::Range(n, 2 * n - 1)],
                             Rcpp::Named("stat") = stat);
}