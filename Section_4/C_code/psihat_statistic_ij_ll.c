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

double biweight_derivative(double x){
	/*
	This function calculates the value of a derivative of the Biweight kernel at point x
	*/
	if (x > 1 || x < -1)
		return(0.0);
	else
		return (-15.0 / 4.0) * x * (1 - x * x);
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


double psi_average_ij_ll(double data_i[], double data_j[], int T, double u, double h){
	/*
	Kernel average function \psi_T(u, h) that takes u, h, data, T=length(data) as arguments. 
	The data can be y_data for \hat{ \Psi }_T or independent gaussian rv z_temp = sigma * z for \Psi^star_T.
	The output is one value for each u and h.
	*/
	int i;
	double result, result_temp, k, k_norm;
	result_temp = 0;
	k_norm = 0;
	k = 0;
	for (i = 0; i < T; i++) {
		k = epanc(((i + 1) / (float)T - u) / h) * (s_t_2(u, h, T) - s_t_1(u, h, T) * (((i + 1) / (float)T - u) / h));
		result_temp += k * (data_i[i] - data_j[i]);
		k_norm += k * k;
/*		Rprintf("We are here, %f, %f, %f\n", k, result_temp, k_norm);*/
	} 		    	
	result = result_temp / sqrt(k_norm);
/*	Rprintf("We are here, %f, %f, %f, %d\n", result, result_temp, k_norm, k_function);*/
	return(result);
}


void psihat_statistic_ij_ll(double *y_data_i, double *y_data_j, int *T, double *g_t_set, int *N, double *sigmahat, double *maximum, double *values, double *values_with_sign){
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
		tmp1 = psi_average_ij_ll(y_data_i, y_data_j, T[0], g_t_set[i], g_t_set[i + N[0]]) / sigmahat[0];
    	values[i] = awert(tmp1) - g_t_set[2 * N[0] + i];
		values_with_sign[i] = tmp1;
/*    	Rprintf("%f, %f, %f, %f\n",tmp1, tmp2, tmp3, values[i]);*/
    	if (i == 0) {
    		maximum[0] = values[i];
    	} else if (values[i] > maximum[0]) {
    		maximum[0] = values[i];
/*    			Rprintf("%f %d\n",max,i);*/
    	}
  	}
}