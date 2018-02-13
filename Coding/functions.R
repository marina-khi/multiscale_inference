# function example - get measures of central tendency
# and spread for a numeric vector x. The user has a
# choice of measures and whether the results are printed.

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

lambda <- function(h)
{
  result = tryCatch(sqrt(2*log(1/(2*h))), warning = function(w) print("h is exceeding h_max"))
#  cat("lambda:", result)
  return(result)
}

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

psihat_average <- function(y_data, u, h, k_function)
{
  T = length(y_data)
  result = sum(omega(T, u, h, k_function)*y_data)
  cat("psihat_average:", result)
  return(result)
}

psihat_statistic <- function(y_data, g_t_set, k_function = epanechnikov_kernel, sigmahat) {
  g_t_set_card = nrow(g_t_set)
  for (i in 1:g_t_set_card) {
    g_t_set$values[i] <- abs(psihat_average(y_data, g_t_set$u[i], g_t_set$h[i], k_function)/sigmahat) - lambda(g_t_set$h[i])
  }
  result = max(g_t_set$values)
  cat("Statistic:", result)
#return(result)
}