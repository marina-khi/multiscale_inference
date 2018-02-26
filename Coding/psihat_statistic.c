/*
Filename: "psihat_statistic.c"
Return a vectors of sequentially summed values
Arguments:
start -- value to start the sum at
size -- the number of elements to return
sumVect -- the vector of summed output values
*/

#include <math.h> 
#include <Rinternals.h> 
#include <R.h>
#include <stdlib.h> 
#include <Rmath.h>

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
	double result, k, k_norm;
	result = 0;
	k_norm = 0;
	if (k_function == 1) {
		for (i = 0; i < T; i++) {
			k = epanc((u - i /(float)T) / h);
			result += k * data[i];
			k_norm += k * k;
		}
	}
	else {
		for (i = 0; i < T; i++) {
			k = biweight_derivative((u - i / (float)T) / h);
			result += k * data[i];
			k_norm += k * k;
		}
	}
	return(result / sqrt(k_norm));
}


void sumSeq(int *start, int *size, int *sumVect){
    /*
    This function provides a simple sequential sum
    where F[n] = F[n-1] + n
    */
    int i, j ;
    j = 0 ;
    for(i = *start; i < (*start + *size); i++){
        if(i == *start){
            sumVect[j] = i ;
        }
        else{
            sumVect[j] = sumVect[j-1] + i ;
        }
        j ++ ;
    }
}

void fiboSeq(int *size, int *sumVect){
    /*
    This function returns the Fibonacci sequence
    where F[n] = F[n-1] + F[n-2]
    */
    int i ;
    sumVect[0] = 0 ;
    sumVect[1] = 1 ;
    for(i = 2; i < *size; i++){
        sumVect[i] = sumVect[i-1] + sumVect[i-2] ;
    }
}