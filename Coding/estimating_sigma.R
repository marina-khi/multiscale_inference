estimating_sigma <- function(y, p, L1, L2){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  T = length(y)
  gamma_hat_zero = 0
  for (r in L1:L2) {
    for (t in (r+1):T) {
      gamma_hat_zero = gamma_hat_zero + 1/(2*(T-r)) * (y[t] - y[t-r]) * (y[t] - y[t-r])
#      cat("gamma_hat_zero:", gamma_hat_zero, "t = ", t, "r = ", r, "\n")
    }
  }
  gamma_hat_zero = gamma_hat_zero / (L2 - L1 + 1)
  
  gamma_hat <- vector(mode = "double",length = T)
  
  for (l in 1:(T-1)) {
    gamma_hat[l] = gamma_hat_zero
    for (j in (l+1):T) {
      gamma_hat[l] = gamma_hat[l] - 1/(T-l) * (y[j] - y[j-l]) * (y[j] - y[j-l])   
    }
  }
  
  
#  sigma2 = sigma_eta2 / ((1 - sum(a_hat))^2)
#  return(sqrt(sigma2))
  return(gamma_hat)
}