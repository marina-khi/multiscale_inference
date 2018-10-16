/*
Filename: "psihat_statistic.c"
Returns value of \Psi^hat_statistic together with the values of
kernel averages \psi_T(u, h) for each u and h from g_t_set.
Arguments are described in the function itself.
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

double s_t_1(double u, double h, int T) {
	int i;
	double result;
	result = 0;
	for (i = 0; i < T; i++) {
		result += epanc(((i + 1) / (float)T - u) / h) * (((i + 1) / (float)T - u) / h);
	}
	return(result / (T * h));
}

double s_t_2(double u, double h, int T) {
	int i;
	double result;
	result = 0;
	for (i = 0; i < T; i++) {
		result += epanc(((i + 1) / (float)T - u) / h) * (((i + 1) / (float)T - u) / h) * (((i + 1) / (float)T - u) / h);
	}
	return(result / (T * h));
}

double s_t_0(double u, double h, int T) {
	int i;
	double result;
	result = 0;
	for (i = 0; i < T; i++) {
		result += epanc(((i + 1) / (float)T - u) / h);
	}
	return(result / (T * h));
}


double local_linear_smoothing(double data[], int T, double u, double h){
	/*
	Kernel average function \psi_T(u, h) that takes u, h, data, T=length(data)
	as arguments. The output is one value for each u and h.
	*/
	int i;
	double x, result, result_temp, k, k_norm;
	result_temp = 0;
	k_norm = 0;
	k = 0;
	for (i = 0; i < T; i++) {
		x = (((i + 1) / (float)T - u) / h);
		k = epanc(x) * (s_t_0(u, h, T) * x - s_t_1(u, h, T));
		result_temp += k * data[i];
		k_norm += k * x;
/*		Rprintf("We are here, %f, %f, %f\n", k, result_temp, k_norm); */   	
	}
	result = result_temp / (h * k_norm);
/*	Rprintf("We are here, %f, %f, %f, %d\n", result, result_temp, k_norm);*/
	return(result);
}


void SiZer_matrix(double *y_data, int *T, double *g_t_set, int *N, double *values_with_sign){
             /* y_data		list of y_t values y= (y_1,...,y_n)
                T        	length of time series
                N        	length of the grid
                g_t_set   	grid vector grid= (u_1,...,u_N, g_1,...,g_N, lambda_1,...lambda_N)
				sigmahat	appropriate estimator for sigma
                maximum		return: value of the test statistic
				values		return: vector of values psi_average for each (u, h) from the grid
			*/
	int i;
	double tmp1;
 	
 	for (i=0; i < N[0]; i++) {
		tmp1 = local_linear_smoothing(y_data, T[0], g_t_set[i], g_t_set[i + N[0]]);
		values_with_sign[i] = tmp1;
/*    	Rprintf("%f, %f, %f, %f\n",tmp1, tmp2, tmp3, values[i]);*/
  	}
}