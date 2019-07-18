#include <Rcpp.h>
#include <math.h> 
#include <R.h>
#include <Rmath.h> 
#include <stdlib.h> 
#include <stdio.h>


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
      result += epanc((t / (float)T - u) / h) * ((t/(float)T - u) / h);
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
      result += epanc((t / (float)T - u) / h) * ((t/(float)T - u) / h) * ((t/(float)T - u) / h);
   }
   return(result / (T * h));
}

// [[Rcpp::export]]
Rcpp::NumericVector sizer_weights(int T, Rcpp::NumericVector gset, int N){

   /* Inputs: 
      T        	  length of time series
      gset        vector of location-bandwidth points (u_1,...,u_N,h_1,...h_N)
      N           number of location-bandwidth points in gset
       
      Outputs: 
      weight      vector of local linear weights ( w_1(u_1,h_1),...,w_T(u_1,h_1),   
                  w_1(u_2,h_2),...,w_T(u_2,h_2),...,w_1(u_N,h_N),...,w_T(u_N,h_N) ) 
   */

  int i, t, pos, pos1, pos2, t_min, t_max;
  double u, h, x;
   
  Rcpp::NumericVector weight1(N * T);
  Rcpp::NumericVector weight2(N);
  Rcpp::NumericVector weight(N * T);
   
  pos1 = 0;
  pos2 = 0; 	

  for(i = 0; i < N; i++){
      u = gset[i];
      h = gset[i+N];
      t_min = (int)(floor((u - h) * (float)T));      
      if(t_min < 1) t_min = 1;
      t_max = (int)(ceil((u + h) * (float)T));      
      if(t_max > T) t_max = T;
      pos1 = i*T;
  
      if(u > h && u < 1.0-h){
        for(t = t_min; t < (t_max+1); t++){
          x = ((t  / (float)T - u) / h);           
          weight1[pos1 + (t_min-1)] = epanc(x) * x;
          weight2[pos2] = weight2[pos2] + epanc(x) * x * x;
	        pos1++;
	      }
      } else {
        for(t = t_min; t < (t_max+1); t++){
          x = ((t / (float)T - u) / h);
          weight1[pos1 + (t_min-1)] = epanc(x) * (s_t_0(u, h, T) * x - s_t_1(u, h, T));
          weight2[pos2] = weight2[pos2] + epanc(x) * (s_t_2(u, h, T) - s_t_1(u, h, T) * x);
	      pos1++;
	      }
      }
      pos2++;
  }
  
  pos = 0;
    
  for(i = 0; i < N; i++){
      for(t = 0; t < T; t++){
 	      weight[pos] =  weight1[pos] / weight2[i];     
        pos++;
      }
    }     
  return(weight);
  free(weight1);
  free(weight2); 
}

