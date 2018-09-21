estimating_sigma_overidentified <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model
  
  y_data_diff <- c(0, diff(y_data))
  
  gamma_vector_hat <- sapply(2:(p+2), FUN = sample_autocovariance, y_data_diff)
  
  Gamma_matrix_hat <- matrix(data = NA, nrow = (p + 1), ncol = p)
  for (i in 1:(p+1)){
    for (j in 1:p){
      if (i  >= j){
        Gamma_matrix_hat[i, j] = sample_autocovariance(i - j + 1, y_data_diff)
      } else {
        Gamma_matrix_hat[i, j] = sample_autocovariance(j - i + 1, y_data_diff)
      }
    }
  }
  
  a_hat     <- solve(t(Gamma_matrix_hat) %*% Gamma_matrix_hat) %*% t(Gamma_matrix_hat) %*% gamma_vector_hat
  return(a_hat)
}