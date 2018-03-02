#This is auxiliary file with lines that were used during writing programs.
library(microbenchmark)

#Auxiliary lines for checking function time
microbenchmark(psihat_statistic(y_data_1, g_t_set, kernel_ind, sigmahat))
microbenchmark(psihat_statistic_old(y_data_1, g_t_set, epanechnikov_kernel, sigmahat))

#Comparing C function with R function 
result <-psihat_statistic(y_data_1, g_t_set, kernel_ind, sigmahat)
g_t_set_with_values <- result[[1]]
psihat_statistic_value <- result[[2]]

#result_old <-psihat_statistic_old(y_data_1, g_t_set, epanechnikov_kernel, sigmahat)
g_t_set_with_values <- result_old[[1]]
psihat_statistic_value <- result_old[[2]]

#a_t_set_1[which.min(a_t_set_1$u),]


estimating_sigma <- function(y, p, L1, L2){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  T = length(y)
  
  #Step 1 as in Section 5.2
  gamma_hat_zero = 0
  gamma_hat <- vector(mode = "double",length = T)
  for (r in L1:L2) {
    for (t in (r+1):T) {
      gamma_hat_zero = gamma_hat_zero + 1/(2 * (T - r) * (L2 - L1 + 1)) * (y[t] - y[t-r]) * (y[t] - y[t-r])
      #      cat("gamma_hat_zero:", gamma_hat_zero, "t = ", t, "r = ", r, "\n")
    }
  }
  
  for (l in 1:(T-1)) {
    gamma_hat[l] = gamma_hat_zero
    for (j in (l+1):T) {
      gamma_hat[l] = gamma_hat[l] - 1/(2 * (T - l)) * (y[j] - y[j-l]) * (y[j] - y[j-l])   
    }
  }
  
  #Step 2 as in Section 5.2
  gamma_matrix <- matrix(, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p){
      if (i!=j) {
        gamma_matrix[i,j] = gamma_hat[abs(i-j)]
      } else {gamma_matrix[i, j] = gamma_hat_zero} 
    }
  }
  
  #  sigma2 = sigma_eta2 / ((1 - sum(a_hat))^2)
  #  return(sqrt(sigma2))
  return(gamma_hat)
}
