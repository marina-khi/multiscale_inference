
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


void sizer_weights(int *T, double *gset, int *N, double *weight){

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
   double *weight1, *weight2; 

   weight1 = (double*) malloc(sizeof(double) * N[0] * T[0]);
   weight2 = (double*) malloc(sizeof(double) * N[0]);
 
   for(i = 0; i < N[0]*T[0]; i++){
      weight1[i] = 0.0;
   }
   for(i = 0; i < N[0]; i++){
      weight2[i] = 0.0;
   }

   pos1 = 0;
   pos2 = 0; 	

   for(i = 0; i < N[0]; i++){
      u = gset[i];
      h = gset[i+N[0]];
      t_min = (int)(floor((u - h) * (float)T[0]));      
      if(t_min < 1) t_min = 1;
      t_max = (int)(ceil((u + h) * (float)T[0]));      
      if(t_max > T[0]) t_max = T[0];
      pos1 = i*T[0];
  
      if(u > h && u < 1.0-h){
        for(t = t_min; t < (t_max+1); t++){
           x = ((t  / (float)T[0] - u) / h);           
           weight1[pos1 + (t_min-1)] = epanc(x) * x;
           weight2[pos2] += epanc(x) * x * x;
	   pos1++;
	}
      } else {
        for(t = t_min; t < (t_max+1); t++){
           x = ((t / (float)T[0] - u) / h);
           weight1[pos1 + (t_min-1)] = epanc(x) * (s_t_0(u, h, T[0]) * x - s_t_1(u, h, T[0]));
           weight2[pos2] += epanc(x) * (s_t_2(u, h, T[0]) - s_t_1(u, h, T[0]) * x);
	   pos1++;
	}
      }
      pos2++;
   }

   pos = 0;
    
   for(i = 0; i < N[0]; i++){
      for(t = 0; t < T[0]; t++){
 	 weight[pos] =  weight1[pos] / weight2[i];     
         pos++;
      }
    }     

   free(weight1);
   free(weight2); 
}
