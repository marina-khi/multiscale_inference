/*
Filename: "estimating_sigma.c"
*/

#include <math.h> 
#include <R.h>
#include <Rmath.h> 
#include <stdlib.h> 
#include <stdio.h>

void estimating_sigma(double *y, int *T, int *L1, int *L2, double *sigmahat){
             /* y_data		list of y_t values y= (y_1,...,y_n)
                T        	length of time series
                N        	length of the grid
                g_t_set   	grid vector grid= (u_1,...,u_N, g_1,...,g_N, lambda_1,...lambda_N)
                k_function	integer that indicates kernel function, 1 = epanechnikov kernel, 2 = derivative of biweight kernel
				sigmahat	appropriate estimator for sigma
                maximum		return: value of the test statistic
				values		return: vector of values psi_average for each (u, h) from the grid
			*/
	int r, t, j;
	double gamma_hat_zero, gamma_hat_one, a_hat_1, sigma2;
	gamma_hat_zero = 0;
 
	/*Step 1 as in Section 5.2*/
	for (r = L1[0]; r < L2[0]; r++) {
		for (t = r + 1; t < T[0]; t++) {
			gamma_hat_zero += 1 / (2 * (T[0] - r) * (L2[0] - L1[0] + 1)) * (y[t - 1] - y[t - r - 1]) * (y[t - 1] - y[t - r - 1]);
			}
		}

	gamma_hat_one = gamma_hat_zero;
	for (j = 2; j < T[0]; j++) {
		gamma_hat_one -= 1 / (2 * (T[0] - 1)) * (y[j - 1] - y[j - 2]) * (y[j - 1] - y[j - 2]);
	}

	/*Step 2 as in Section 5.2 */
	a_hat_1 = gamma_hat_one / gamma_hat_zero;
	
	/*Step 3 as in Section 5.2*/
	sigma2 = (gamma_hat_zero * (1 - a_hat_1 * a_hat_1)) / ((1 - a_hat_1) * (1 - a_hat_1));
	sigmahat[0] = sqrt(sigma2);
}