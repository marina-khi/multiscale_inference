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
Rcpp::NumericVector kernel_averages(int T, Rcpp::NumericVector gset, Rcpp::NumericVector data, double sigmahat, int N){

   /* Inputs: 
      T        	length of time series
      gset        vector of location-bandwidth points (u_1,...,u_N,h_1,...h_N)
      data        time series of length T (y_1, ..., y_T)
      N           number of location-bandwidth points in gset
       
      Outputs: 
      kernel averages
                  vector of kernel averages ( \phi(u_1,h_1), \phi(u_2,h_2),...,\phi(u_N,h_N) )  
   */

   int i, t, t_min, t_max;
   double u, h, x, temp, result_temp, weight_norm;

   Rcpp::NumericVector result(N);
   
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
   }

    return result;
}