# function example - get measures of central tendency
# and spread for a numeric vector x. The user has a
# choice of measures and whether the results are printed.

lambda <-function(h)
{
  result<- sqrt(2*log(1/(2*h)))
  return(result)
}

omega <- function(T, u, h, k_function)
{
  result=c()
  for (i in 1:T) { 
    result[i]=k_function((u-i/T)/h)/K_norm
  }
  return(result)
}

psihat_average <- function(y_data, u, h, k_function)
{
  T <- length(y_data)
  result<-sum(omega(T, u, h, k_function)*y_data)
  return(result)
}

psihat_statistic <- function(y_data, g_t_set, k_function = epanechnikov, sigmahat) {
  g_t_set$values = abs(psihat_average(y_data, g_t_set$u, g_t_set$h, k_function/sigmahat) - lambda(g_t_set$h)
  result <- max(g_t_set$values)
  cat("Statistic", result)
  return(result)
}