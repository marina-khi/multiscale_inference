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

# Additive correction tern \lambda(h)
# that depends only on the bandwidth h

lambda <- function(h)
{
  result = tryCatch(sqrt(2*log(1/(2*h))), warning = function(w) print("h is exceeding h_max"))
#  cat("lambda:", result)
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
  }
  K_norm = sqrt(K_norm_temp)
  result = result_temp / K_norm
#  cat("K_norm:", K_norm, " omega_t:", result)
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
#  cat("psihat_average:", result)
  return(result)
}

# Function that calculates the multiscale statistic \hat{\Psi}_T.
# It takes the following entities as arguments.
#   y-data: the data
#   g_t_set: range of different locations u and bandwidths h
#   k_function: type of kernel function
#   sigmahat: the estimator of the square root of the long-run error variance \sigma^2
# It produces the value of the test statistic as an output

psihat_statistic <- function(y_data, g_t_set, k_function = epanechnikov_kernel, sigmahat) {
  g_t_set_card = nrow(g_t_set)
  for (i in 1:g_t_set_card) {
    g_t_set[['values']][i] <- abs(psi_average(y_data, .subset2(g_t_set,'u')[i], .subset2(g_t_set,'h')[i], k_function)/sigmahat) - .subset2(g_t_set, 'lambda')[i]
  }
  result = max(g_t_set$values)
#  cat("Statistic:", result)
  return(result)
}

psihat_statistic_temp <- function(y_data, g_t_set, kernel_function = epanechnikov_kernel, sigmahat) {
  g_t_set$values <- abs(mapply(psi_average, u = g_t_set$u, h = g_t_set$h,
                                    MoreArgs = list(data = y_data, k_function = kernel_function))/sigmahat) - .subset2(g_t_set, 'lambda')
  result = max(g_t_set$values)
  return(result)
}

# Function that calculates the auxiliary statistic \Psi^star_T.
# It takes the following entities as arguments.
#   z_t: the independent standard normal random variables
#   g_t_set: range of different locations u and bandwidths h
#   k_function: type of kernel function
#   sigma: the square root of the long-run error variance \sigma^2
# It produces the value of the test statistic as an output
# which is used further to calculate the critical values

psistar_statistic <- function(z, g_t_set, k_function = epanechnikov_kernel, sigma) {
  g_t_set_card = nrow(g_t_set)
  z_temp = sigma*z
  for (i in 1:g_t_set_card) {
    g_t_set[['values']][i] <- abs(psi_average(z_temp, g_t_set[['u']][i], g_t_set[['h']][i], k_function)/sigma) - lambda(g_t_set[['h']][i])
  }
  result = max(g_t_set$values)
#  cat("Statistic with a star:", result)
  return(result)
}