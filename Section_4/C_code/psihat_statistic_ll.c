/*
Filename: "psihat_statistic_ij.c"
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
	double x, result, result_temp, k, k_norm;
	result_temp = 0;
	k_norm = 0;
	k = 0;
	for (i = 0; i < T; i++) {
		x = (((i + 1) / (float)T - u) / h);
		k = epanc(x) * (s_t_2(u, h, T) - s_t_1(u, h, T) * x);
		result_temp += k * (data_i[i] - data_j[i]);
		k_norm += k * k;
/*		Rprintf("We are here, %f, %f, %f\n", k, result_temp, k_norm);*/
	} 		    	
	result = result_temp / sqrt(k_norm);
/*	Rprintf("We are here, %f, %f, %f, %d\n", result, result_temp, k_norm, k_function);*/
	return(result);
}


double psihat_statistic_ij_ll(double *y_data_i, double *y_data_j, int T, double *g_t_set, int N, double sigmahat){
	int i;
	double tmp1, maximum, values[N];
 	
 	for (i=0; i < N; i++) {
		tmp1 = psi_average_ij_ll(y_data_i, y_data_j, T, g_t_set[i], g_t_set[i + N]) / sigmahat;
    	values[i] = awert(tmp1) - g_t_set[2 * N + i];
/*		values_with_sign[i] = tmp1;*/
/*    	Rprintf("%f, %f, %f, %f\n",tmp1, tmp2, tmp3, values[i]);*/
    	if (i == 0) {
    		maximum = values[i];
    	} else if (values[i] > maximum) {
    		maximum = values[i];
/*    			Rprintf("%f %d\n",max,i);*/
    	}
  	}
  	return(maximum);
}

void psihat_statistic_ll(double *y_data, int *T, double *g_t_set, int *N, int *N_ts, double *sigmahat, double *statistic, double *statistic_result){
             /* y_data		list of y_t values y= (y_1,...,y_n)
                T        	length of time series
                N        	length of the grid
                N_ts		numbet of time series
                g_t_set   	grid vector grid= (u_1,...,u_N, g_1,...,g_N, lambda_1,...lambda_N)
				sigmahat	appropriate estimator for sigma
                maximum		return: value of the test statistic
				values		return: vector of values psi_average for each (u, h) from the grid
			*/
	int i, j, k, n;
	double y_data_i[T[0]], y_data_j[T[0]];
 	
 	k = 0;
 	for (i = 0; i < N_ts[0] - 1; i++) {
 		for (j = i + 1; j < N_ts[0]; j++){
 			for (n = 0; n < T[0]; n++){
 				y_data_i[n] = y_data[i * T[0] + n];
 				y_data_j[n] = y_data[j * T[0] + n];
 			}
 			statistic[k] = psihat_statistic_ij_ll(y_data_i, y_data_j, T[0], g_t_set, N[0], sigmahat[0]);
 			if (k == 0) {
    			statistic_result[0] = statistic[k];
    		} else if (statistic[k] > statistic_result[0]) {
    			statistic_result[0] = statistic[k];
/*    			Rprintf("%f %d\n",max,i);*/
    		}
    		k ++;
 		}
  	}
}