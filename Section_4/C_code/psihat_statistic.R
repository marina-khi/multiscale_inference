psihat_statistic <- function(y_data, N_ts, g_t_set, kernel_ind, sigmahat){
  #Wrapper of a C function from psihat_statistic.C
  #Function that calculates the multiscale statistic \hat{\Psi}_T.
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # kernel_ind      type of kernel function used - 1 for epanechnikov, 2 for biweight
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T        <- as.integer(nrow(y_data))
  N        <- as.integer(nrow(g_t_set))
  sigmahat <- as.double(sigmahat)
  kernel   <- as.integer(kernel_ind)
  N_ts     <- as.integer(N_ts)
  
  storage.mode(sigmahat)    <- "double"
  
  g_t_set_vec       <- unlist(g_t_set[c('u', 'h', 'lambda')])
  y_data_vec        <- unlist(y_data)
  values            <- vector(mode = "double", length = N)
  values_with_sign  <- vector(mode = "double", length = N)
  maximum           <- vector(mode = "double", length = 1)
  maximum_ij        <- vector(mode = 'double', length = N_ts * (N_ts - 1) / 2)
  
  result <- .C("psihat_statistic", y_data_vec, T, N_ts, g_t_set_vec, N, kernel, sigmahat, maximum_ij, maximum, values, values_with_sign)
  
  #  cat("kernel = ", result[[5]], ", T=", result[[2]], ", N=", result[[4]], ", sigmahat = ", result[[6]],
  #      "maximum = ", result[[7]], "len of values = ", length(result[[8]]), sep=" ")
  
  statistic_vector         <- result[[8]]
  statistic                <- result[[9]]
  
  return(list(statistic_vector, statistic)) 
}

psistar_statistic <- function(z_data, N_ts, g_t_set, kernel_ind = 1, sigmahat){
  #Wrapper of a C function from psihat_statistic.C
  #Function that calculates the auxiliary statistic \Psi^star_T.
  #The only difference with the previous function is in the return values.
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # kernel_ind      type of kernel function used - 1 for epanechnikov, 2 for biweight
  # sigmahat        appropriate estimate of sqrt(long-run variance)


  T        <- as.integer(nrow(z_data))
  N        <- as.integer(nrow(g_t_set))
  sigmahat <- as.double(sigmahat)
  kernel   <- as.integer(kernel_ind)
  N_ts     <- as.integer(N_ts)
  
  storage.mode(sigmahat)    <- "double"
  
  g_t_set_vec       <- unlist(g_t_set[c('u', 'h', 'lambda')])
  z_data_vec        <- unlist(z_data)
  values            <- vector(mode = "double", length = N)
  values_with_sign  <- vector(mode = "double", length = N)
  maximum           <- vector(mode = "double", length = 1)
  maximum_ij        <- vector(mode = 'double', length = N_ts * (N_ts - 1) / 2)
  
  result <- .C("psihat_statistic", z_data_vec, T, N_ts, g_t_set_vec, N, kernel, sigmahat, maximum_ij, maximum, values, values_with_sign)
  
  #  cat("kernel = ", result[[5]], ", T=", result[[2]], ", N=", result[[4]], ", sigmahat = ", result[[6]],
  #      "maximum = ", result[[7]], "len of values = ", length(result[[8]]), sep=" ")
  
  statistic                <- result[[9]]
  return(statistic) 
}