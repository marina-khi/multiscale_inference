#include <Rcpp.h>
#include <math.h> 
#include <R.h>
#include <Rmath.h> 
#include <stdlib.h> 
#include <stdio.h>

// [[Rcpp::export]]
Rcpp::List EstimateLrvHvK(Rcpp::NumericVector y, int T, int L1, int L2){
             /* Function that estimates AR(1) parameters of time series y based on the method from Hall and van Keilegom (2003).
              Args:
                y			    list of y_t values y= (y_1,...,y_T)
                T         length of time series
                L1, L2    Tuning parameters, integers
              Returns:
				        sigma_hat: estimator for square root of long-run error variance, sigma
				        a_hat_1:   estimator for a_1 coefficient of AR(1)
				        sigma_eta: estimator for square root of the variance of the innovation in AR(1) model
			*/
	int r, t, j;
	double gamma_hat_zero, gamma_hat_one;
	gamma_hat_zero = 0;
 
	/*Step 1 as in Section 5.2*/
	for (r = L1; r <= L2; r++) {
		for (t = r + 1; t <= T; t++) {
			gamma_hat_zero = gamma_hat_zero + 1 / (2 * (float)(T - r) * (float)(L2 - L1 + 1)) * (y[t - 1] - y[t - r - 1]) * (y[t - 1] - y[t - r - 1]);
		}
	}

	gamma_hat_one = gamma_hat_zero;
	for (j = 2; j <= T; j++) {
		gamma_hat_one = gamma_hat_one - 1 / (2 * (float)(T - 1)) * (y[j - 1] - y[j - 2]) * (y[j - 1] - y[j - 2]);
	}

	/*Step 2 as in Section 5.2 */
	double a_hat_1 = gamma_hat_one / gamma_hat_zero;
		
	/*Step 3 as in Section 5.2*/
	double sigma_eta = sqrt(gamma_hat_zero * (1 - a_hat_1 * a_hat_1));
	double sigma2 =  (sigma_eta * sigma_eta)  / ((1 - a_hat_1) * (1 - a_hat_1));
	double sigma_hat = sqrt(sigma2);
	return Rcpp::List::create(Rcpp::Named("a_hat") = a_hat_1,
                           Rcpp::Named("sigma_eta") = sigma_eta,
                           Rcpp::Named("sigma_hat") = sigma_hat);
}