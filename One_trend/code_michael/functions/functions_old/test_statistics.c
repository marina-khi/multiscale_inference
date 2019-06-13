
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

   int t;
   double result;
   result = 0;
   for(t = 0; t < T; t++){
      result += epanc(((t + 1) / (float)T - u) / h);
   }
   return(result / (T * h));
}


double s_t_1(double u, double h, int T){

   /* This function computes the kernel constant S_T1 needed for the local linear weights 
      of the kernel average hat(psi)_T(u,h) */ 

   int t;
   double result;
   result = 0;
   for(t = 0; t < T; t++){
      result += epanc(((t + 1) / (float)T - u) / h) * (((t + 1) / (float)T - u) / h);
   }
   return(result / (T * h));
}


double psi_statistic(double data[], int T, double u, double h){

   /* This function computes the kernel average hat(psi)_T(u, h) for a given point (u,h) */
  
   int t;
   double x, weight, nom, denom, result;
   nom = 0;
   denom = 0;

   if(u > h && u < 1-h) {
     for(t = 0; t < T; t++){
        x = (((t + 1) / (float)T - u) / h);
        weight = epanc(x) * x;
      nom += weight * data[t];
      denom += weight * weight;
     }
   } else {
     for(t = 0; t < T; t++){
        x = (((t + 1) / (float)T - u) / h);
        weight = epanc(x) * (s_t_0(u, h, T) * x - s_t_1(u, h, T));
        nom += weight * data[t];
        denom += weight * weight;
     }
   }
   result = nom/sqrt(denom);

   return(result);
}


void test_statistics(double *data, int *T, double *gset, int *N, double *h_vec, int *N_h, int *N_u, double *correct, double *sigmahat, double *values, double *Psi_multiscale, double *Psi_uncorrected, double *Psi_rowwise){

   /* Inputs: 
      data	  time series observations, vector of length T
      T        	  length of time series
      gset        vector of location-bandwidth points (u_1,...,u_N,h_1,...h_N)
      N           number of location-bandwidth points in gset
      h_vec       vector of bandwidth levels
      N_h         length of vector h_vec
      N_u         vector of length N_h, N_u[i] specifies the number of 
                  location-bandwidth points (u,h) with h=h_vec[i]  
      correct     vector of additive correction terms lambda of length N_h      
      sigmahat	  estimator for sqrt of long-run error variance 

      Outputs: 
      values            vector of individual test statistics (psi_statistic/sigmahat) 
                        for each (u,h) from gset 
      Psi_multiscale    value of the multiscale statistic
      Psi_uncorrected   value of the multiscale statistic (uncorrected version)
      Psi_rowwise       vector of length N_h, Psi_rowwise[i] is the maximum of the 
                        individual test statistics abs(psi_statistic/sigmahat)
                        corresponding to the bandwidth h=h_vec[i]
   */

   int i, j, pos;
   double temp;
   double *abs_values, *abs_values_corrected;
             
   abs_values = (double*) malloc(sizeof(double) * N[0]);
   abs_values_corrected = (double*) malloc(sizeof(double) * N[0]);

   pos = 0;
 	
   for(i=0; i < N_h[0]; i++){
      for(j=0; j < N_u[i]; j++){
 	 temp = psi_statistic(data, T[0], gset[pos], gset[pos+N[0]]) / sigmahat[0];
         values[pos] = temp;
         abs_values[pos] = awert(temp);
         abs_values_corrected[pos] = awert(temp) - correct[i];
         
         if(i == 0 && j == 0){
           Psi_multiscale[0] = abs_values_corrected[0]; 
	 } else if (abs_values_corrected[pos] > Psi_multiscale[0]) {
           Psi_multiscale[0] = abs_values_corrected[pos];
         } 

         if(i == 0 && j == 0){
           Psi_uncorrected[0] = abs_values[0]; 
	 } else if (abs_values[pos] > Psi_uncorrected[0]) {
           Psi_uncorrected[0] = abs_values[pos];
         } 

         if(j == 0){
           Psi_rowwise[i] = abs_values[pos];
	 } else if (abs_values[pos] > Psi_rowwise[i]) {
           Psi_rowwise[i] = abs_values[pos];
         }         

         pos++; 
      }
   }      
}
