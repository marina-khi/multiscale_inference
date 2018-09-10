estimating_sigma_new_method <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model    
  gamma_hat_temp <- sapply(0:(p+1), FUN = sample_autocovariance, y_data)
  if (p == 1){
    Gamma_matrix_hat <- gamma_hat_temp[2]
    gamma_hat        <- gamma_hat_temp[3]
    a_hat            <- gamma_hat / Gamma_matrix_hat     
  } else {
    Gamma_matrix_hat <- matrix(data = NA, nrow = p, ncol = p)
    for (i in 1:p){
      for (j in 1:p){
        if (i >= j) {
            Gamma_matrix_hat[i, j] = gamma_hat[i - j + 2]
        } else {
            Gamma_matrix_hat[i, j] = gamma_hat[j - i]
          }
      }
    }
    gamma_hat <-  gamma_hat_temp[3:(p+2)]
    a_hat     <- inv(Gamma_matrix_hat) %*% gamma_hat
  }
  sigma_eta_squared <- sum(a_hat * gamma_hat_temp[2:(p+1)]) - gamma_hat_temp[1]
  sigma_eta <- sqrt(sigma_eta_squared)
  sigma_hat <- sigma_eta / (abs(1 - sum(a_hat)))
  return(list(sigma_hat, a_hat, sigma_eta))
}