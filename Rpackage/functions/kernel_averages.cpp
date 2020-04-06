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


// [[Rcpp::export]]
Rcpp::List kernel_averages(int T, Rcpp::NumericVector gset, Rcpp::NumericVector correct, Rcpp::NumericVector data, double sigmahat, int N){

   /* Inputs: 
      T        	  length of time series
      gset        vector of location-bandwidth points (u_1,...,u_N,h_1,...h_N)
      correct     vector of corrections of lenghth N (lambda(h_1), ..., lambda(h_N))
      data        time series of length T (y_1, ..., y_T)
      N           number of location-bandwidth points in gset
       
      Output:
      list with the following components
      vals        vector of normalised kernel averages 
                  (\phi(u_1,h_1)/sigmahat, \phi(u_2,h_2)/sigmahat,...,\phi(u_N,h_N)/sigmahat)
      vals_cor    vector of absolute value of normalised kernel averages with correction 
                  (abs(\phi(u_1,h_1)/sigmahat) - \lambda(h_1),...,abs(\phi(u_N,h_N)/sigmahat) - \lambda(h_N))
      stat        multiscale statistic calculated as max(vals_cor)
   */

   int i, t, t_min, t_max;
   double u, h, x, temp, result_temp, weight_norm;

   NumericVector result(N);
   NumericVector result_cor(N);
   NumericVector stat_ms(1);
   
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
      result_cor[i] = awert(result[i]) - correct[i];
      if(i == 0){
        stat_ms[0] = result_cor[i];
      } else if (stat_ms[0] < result_cor[i]){
        stat_ms[0] = result_cor[i];
      }
   }
  return Rcpp::List::create(Rcpp::Named("vals") = result,
                            Rcpp::Named("vals_cor") = result_cor,
                            Rcpp::Named("stat") = stat_ms);
}