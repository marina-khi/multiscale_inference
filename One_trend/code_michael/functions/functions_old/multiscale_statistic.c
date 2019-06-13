
/* The two main functions compute
(1) the values of the kernel averages hat(psi)_T(u,h) for each (u,h) in grid
(2) the value of the multiscale statistic hat(Psi).
The functions are described in more detail below. */


#include <math.h> 
#include <R.h>
#include <Rmath.h> 
#include <stdlib.h> 
#include <stdio.h>


/* Auxiliary functions */


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


/* Main functions */


double psi_average(double data[], int T, double u, double h){

   /* This function computes the kernel average hat(psi)_T(u, h) 

      Inputs: 
      data    time series of observations
      T       sample size
      u       location, point in [0,1]
      h       bandwidth

      Output: 
      result  value of hat(psi)_T(u, h) 
   */

   int t;
   double x, weight, nom, denom, result;
   nom = 0;
   denom = 0;

   for(t = 0; t < T; t++){
      x = (((t + 1) / (float)T - u) / h);
      weight = epanc(x) * (s_t_0(u, h, T) * x - s_t_1(u, h, T));
      nom += weight * data[t];
      denom += weight * weight;
   }
   
   result = nom/sqrt(denom);

   return(result);
}


void multiscale_statistic(double *data, int *T, double *grid, double *correct, int *N, double *sigmahat, double *psi_values, double *test_values, double *max_value){

   /* Inputs: 
      data	   time series observations, vector of length T
      T        	   length of time series
      grid         grid vector (u_1,...,u_N, h_1,...,h_N)
      correct      vector of additive correction terms (lambda_1,...lambda_N)      
      N        	   length of the grid
      sigmahat	   estimator for sqrt of long-run error variance 

      Outputs: 
      psi_values   vector of individual test statistics (psi_average/sigmahat) for each (u,h) from the grid 
      test_values  vector of individual test statistics (|psi_average/sigmahat| - correct) for each (u,h) from the grid
      max_value    overall value of the multiscale statistic
 
   */

   int i;
   double vals;
 	
   for(i=0; i < N[0]; i++){
      vals = psi_average(data, T[0], grid[i], grid[i + N[0]]) / sigmahat[0];
      psi_values[i] = vals;
      test_values[i] = awert(vals) - correct[i];
      if(i == 0){
     	max_value[0] = test_values[0];
      } else if (test_values[i] > max_value[0]) {
    	max_value[0] = test_values[i];
      }
   }
}
