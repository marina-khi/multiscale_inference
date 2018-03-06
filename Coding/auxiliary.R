#This is auxiliary file with lines that were used during writing programs.
#library(microbenchmark)

#Auxiliary lines for checking function time
#microbenchmark(psihat_statistic(y_data_1, g_t_set, kernel_ind, sigmahat))
#microbenchmark(psihat_statistic_old(y_data_1, g_t_set, epanechnikov_kernel, sigmahat))

#Comparing C function with R function 
#result <-psihat_statistic(y_data_1, g_t_set, kernel_ind, sigmahat)
#g_t_set_with_values <- result[[1]]
#psihat_statistic_value <- result[[2]]

#result_old <-psihat_statistic_old(y_data_1, g_t_set, epanechnikov_kernel, sigmahat)
#g_t_set_with_values <- result_old[[1]]
#psihat_statistic_value <- result_old[[2]]

#a_t_set_1[which.min(a_t_set_1$u),]

#Adding a function that is const on the whole interval
#const <- sqrt(sigma*sigma/noise_to_signal)
#y_data_2 <- const + rnorm(T, 0, sigma)

estimating_sigma_for_AR1_old_version <- function(y, L1, L2){
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
  return(sqrt(sigma2))
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
  return(sqrt(sigma2))
}


# Function that calculates the multiscale statistic \hat{\Psi}_T.
# It takes the following entities as arguments.
#   y-data: the data
#   g_t_set: range of different locations u and bandwidths h
#   k_function: type of kernel function
#   sigmahat: the estimator of the square root of the long-run error variance \sigma^2
# It produces the value of the test statistic as an output
#psihat_statistic_old <- function(y_data, g_t_set, kernel_function = epanechnikov_kernel, sigmahat) {
#  g_t_set$values <- abs(mapply(psi_average, u = g_t_set$u, h = g_t_set$h,
#                                    MoreArgs = list(data = y_data, k_function = kernel_function))/sigmahat) - .subset2(g_t_set, 'lambda')
#  result = max(g_t_set$values)
#  return(list(g_t_set, result))
#}

# Function that calculates the auxiliary statistic \Psi^star_T.
# The only difference with the previous function is in the return values.
#psistar_statistic_old <- function(y_data, g_t_set, kernel_function = epanechnikov_kernel, sigmahat) {
#  g_t_set$values <- abs(mapply(psi_average, u = g_t_set$u, h = g_t_set$h,
#                               MoreArgs = list(data = y_data, k_function = kernel_function))/sigmahat) - .subset2(g_t_set, 'lambda')
#  result = max(g_t_set$values)
#  return(result)
#}

# Epanechnikov kernel function,
# which is defined as f(x) = 3/4(1-x^2)
# for |x|<=1 and 0 elsewhere
epanechnikov_kernel <- function(x)
{
  if (abs(x)<=1)
  {
    result = 3/4 * (1 - x*x)
  } else {
    result = 0
  }
  return(result)
}


# Additive correction tern \lambda(h) that depends only on the bandwidth h
lambda <- function(h)
{
  result = tryCatch(sqrt(2*log(1/(2*h))), warning = function(w) print("h is exceeding h_max"))
  return(result)
}


# Kernel weight function that calculates
# the vector of (\omega_{1,T}(u,h), \ldots, \omega_{T,T}(u,h)).
# Meanwhile, the function also calculates
# the norm ||K||_{u,h,T} of the kernel function.
omega <- function(T, u, h, k_function)
{
  result_temp = c()
  K_norm_temp = 0 
  for (i in 1:T) {
    x = (u - i/T)/h
    k = k_function(x)
    result_temp[i] = k
    K_norm_temp = K_norm_temp + k^2
    #    cat("k:", k, " omega_t_temp:", result_temp[i], "i = ", i, "\n")
  }
  K_norm = sqrt(K_norm_temp)
  result = result_temp / K_norm
  
  return(result)
}


# Kernel average function \psi_T(u,h) that takes u, h, data
# and the type of kernel function as arguments.
# The data can be y_data for \hat{\Psi}_T or independent gaussian rv
# z_temp = sigma*z for \Psi^star_T.
# The output is one value for each u and h.
psi_average <- function(data, u, h, k_function)
{
  T = length(data)
  result = sum(omega(T, u, h, k_function)*data)
  return(result)
}
