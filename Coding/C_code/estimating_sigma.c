/*
Filename: "estimating_sigma.c"
*/

#include <math.h> 
#include <R.h>
#include <Rmath.h> 
#include <stdlib.h> 
#include <stdio.h>

void estimating_sigma(double *y, int *T, int *L1, int *L2, double *sigmahat){
             /* y			list of y_t values y= (y_1,...,y_T)
                T        	length of time series
                L1, L2     	Tuning parameters, integers
				sigmahat	return: estimator for square root of long-run error variance, sigma
			*/
	int r, t, j;
	double gamma_hat_zero, gamma_hat_one, a_hat_1, sigma2;
	gamma_hat_zero = 0;
 
	/*Step 1 as in Section 5.2*/
	for (r = L1[0]; r <= L2[0]; r++) {
		for (t = r + 1; t <= T[0]; t++) {
			gamma_hat_zero += 1 / (2 * (float)(T[0] - r) * (float)(L2[0] - L1[0] + 1)) * (y[t - 1] - y[t - r - 1]) * (y[t - 1] - y[t - r - 1]);
			}
		}

	gamma_hat_one = gamma_hat_zero;
	for (j = 2; j <= T[0]; j++) {
		gamma_hat_one -= 1 / (2 * (float)(T[0] - 1)) * (y[j - 1] - y[j - 2]) * (y[j - 1] - y[j - 2]);
	}

	/*Step 2 as in Section 5.2 */
	a_hat_1 = gamma_hat_one / gamma_hat_zero;
	
	/*Step 3 as in Section 5.2*/
	sigma2 = (gamma_hat_zero * (1 - a_hat_1 * a_hat_1)) / ((1 - a_hat_1) * (1 - a_hat_1));
	sigmahat[0] = sqrt(sigma2);
}