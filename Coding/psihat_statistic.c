/*
Filename: "psihat_statistic.c"
Return a vectors of sequentially summed values
Arguments:
start -- value to start the sum at
size -- the number of elements to return
sumVect -- the vector of summed output values
*/

#include <math.h> 
#include <R.h>
#include <Rmath.h> 
#include <stdlib.h> 
#include <stdio.h>

double awert(double x){
	/*
	This function provides an absolute value of a double x
	*/

	if (x < 0)
		return(-x);
	else
		return(x);
}

double epanc(double x){
	/*
	This function calculates the value of a Epanechnikov kernel at point x
	*/
	if (x > 1 || x < -1)
		return(0.0);
	else
		return (3.0 / 4.0) * (1 - x * x);
}

double biweight(double x){
	/*
	This function calculates the value of a Biweight kernel at point x
	*/
	if (x > 1 || x < -1)
		return(0.0);
	else
		return (15.0 / 16.0) * (1 - x * x) * (1 - x * x);
}

double biweight_derivative(double x){
	/*
	This function calculates the value of a derivative of the Biweight kernel at point x
	*/
	if (x > 1 || x < -1)
		return(0.0);
	else
		return (-15.0 / 4.0) * x * (1 - x * x);
}


double psi_average(double data[], int T, double u, double h, int k_function){
	/*
	Kernel average function \psi_T(u, h) that takes u, h, data, T=length(data)
	and the type of kernel function as arguments. 
	The data can be y_data for \hat{ \Psi }_T or independent gaussian rv z_temp = sigma * z for \Psi^star_T.
	The output is one value for each u and h.
	*/
	int i;
	double result, result_temp, k, k_norm;
	result_temp = 0;
	k_norm = 0;
	k = 0;
	if (k_function == 1) {
		for (i = 0; i < T; i++) {
			k = epanc((u - i /(float)T) / h);
			result_temp += k * data[i];
			k_norm += k * k;
/*			Rprintf("We are here, %f, %f, %f\n", k, result_temp, k_norm);*/
		} 		    	
	}
	else {
		for (i = 0; i < T; i++) {
			k = biweight_derivative((u - i / (float)T) / h);
			result_temp += k * data[i];
			k_norm += k * k;
/*	   		Rprintf("We are here, %f, %f, %f\n", k, result_temp, k_norm); */   	
		}
	}
	result = result_temp / sqrt(k_norm);
/*	Rprintf("We are here, %f, %f, %f, %d\n", result, result_temp, k_norm, k_function);*/
	return(result);
}


void psihat_statistic(double *y_data, int *T, double *g_t_set, int *N, int *k_function, double *sigmahat, double *maximum, double *values){
             /* y_data		list of y_t values y= (y_1,...,y_n)
                T        	length of time series
                N        	length of the grid
                g_t_set   	grid vector grid= (u_1,...,u_N, g_1,...,g_N, lambda_1,...lambda_N)
                k_function	integer that indicates kernel function, 1 = epanechnikov kernel, 2 = derivative of biweight kernel
                maximum		return: value of the test statistic
				values		return: vector of values psi_average for each (u, h) from the grid
			*/
	int i, kernel_ind, ndim;
	double max;
    
    kernel_ind = k_function[0];
    ndim = N[0];
 	
 	for (i=0; i < ndim; i++) {
    	values[i] = awert(psi_average(y_data, T[0], g_t_set[i], g_t_set[i + ndim], kernel_ind) / sigmahat[0]) - g_t_set[2 *ndim +i];
/*    	Rprintf("%f, %f, %f, %f\n",tmp1, tmp2, tmp3, values[i]);*/
    	if (i == 0) {
    		max = values[i];
    	} else if (values[i] > max) {
    			max = values[i];
/*    			Rprintf("%f %d\n",max,i);*/
    	}
  	}
  	maximum[0] = max;	
}