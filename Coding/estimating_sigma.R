estimating_sigma_for_AR1 <- function(y, L1, L2){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  T = length(y)
  
  #Step 1 as in Section 5.2
  gamma_hat_zero = 0
  gamma_hat <- vector(mode = "double",length = 1)
  for (r in L1:L2) {
    for (t in (r+1):T) {
      gamma_hat_zero = gamma_hat_zero + 1/(2 * (T - r) * (L2 - L1 + 1)) * (y[t] - y[t-r]) * (y[t] - y[t-r])
    }
  }
  
  gamma_hat[1] = gamma_hat_zero
  for (j in 2:T) {
    gamma_hat[1] = gamma_hat[1] - 1/(2 * (T - 1)) * (y[j] - y[j-1]) * (y[j] - y[j-1])   
  }

#  for (l in 1:L1) {
#    gamma_hat[l] = gamma_hat_zero
#    for (j in (l+1):T) {
#      gamma_hat[l] = gamma_hat[l] - 1/(2 * (T - l)) * (y[j] - y[j-l]) * (y[j] - y[j-l])   
#    }
#  }
  
  #Step 2 as in Section 5.2
  #cat("gamma_hat_zero:", gamma_hat_zero, "gamma_hat[1] = ", gamma_hat[1], "\n")
  a_hat_1 = gamma_hat[1] / gamma_hat_zero

  #Step 3 as in Section 5.2
  sigma_eta2 = gamma_hat_zero * (1 - a_hat_1 * a_hat_1)
  sigma2 = sigma_eta2 / ((1 - a_hat_1) * (1 - a_hat_1))
  return(sigma2)
}

estimating_sigma_for_AR1_with_param <- function(y, L1, L2){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  T = length(y)
  
  #Step 1 as in Section 5.2
  gamma_hat_zero = 0
  for (r in L1:L2) {
    for (t in (r+1):T) {
      gamma_hat_zero = gamma_hat_zero + 1/(2 * (T - r) * (L2 - L1 + 1)) * (y[t] - y[t-r]) * (y[t] - y[t-r])
    }
  }
  
  #Step 2
  fittedmodel = ar.ols(x = y, order.max = 1, demean = FALSE, intercept = FALSE)
  a_hat_1 = fittedmodel$ar[1]

  #Step 3 as in Section 5.2
  sigma_eta2 = gamma_hat_zero * (1 - a_hat_1 * a_hat_1)
  sigma2 = sigma_eta2 / ((1 - a_hat_1) * (1 - a_hat_1))
  return(sigma2)
}