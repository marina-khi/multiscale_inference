psihat_statistic_ij <- function(y_data_i, y_data_j, g_t_set, sigmahat, kernel_method){
  #Wrapper of a C function from psihat_statistic.C
  #Function that calculates the multiscale statistic \hat{\Psi}_T.
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T        <- as.integer(length(y_data_i))
  N        <- as.integer(nrow(g_t_set))

  storage.mode(y_data_i)    <- "double"
  storage.mode(y_data_j)    <- "double"
  storage.mode(sigmahat)    <- "double"
  
  g_t_set_vec       <- unlist(g_t_set[c('u', 'h', 'lambda')])
  values            <- vector(mode = "double",length = N)
  values_with_sign  <- vector(mode = "double",length = N)
  maximum           <- vector(mode = "double", length = 1)
  
  if (kernel_method == "nw"){
    result <- .C("psihat_statistic_ij", y_data_i, y_data_j, T, g_t_set_vec, N, sigmahat, maximum, values, values_with_sign)
  } else if (kernel_method == "ll"){
    result <- .C("psihat_statistic_ij_ll", y_data_i, y_data_j, T, g_t_set_vec, N, sigmahat, maximum, values, values_with_sign)
  } else {
    print('Given method is currently not supported')
  }

  g_t_set$values           <- result[[8]]
  g_t_set$values_with_sign <- result[[9]]
  statistic                <- result[[7]]

  return(list(g_t_set, statistic)) 
}

psistar_statistic_ij <- function(z_data_i, z_data_j, g_t_set, sigmahat, kernel_method){
  #Wrapper of a C function from psihat_statistic.C
  #Function that calculates the auxiliary statistic \Psi^star_T.
  #The only difference with the previous function is in the return values.
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T                         <- as.integer(length(z_data_i))
  N                         <- as.integer(nrow(g_t_set))
  storage.mode(z_data_i)    <- "double"
  storage.mode(z_data_j)    <- "double"
  
  storage.mode(sigmahat)    <- "double"
  
  g_t_set_vec       <- unlist(g_t_set[c('u', 'h', 'lambda')])
  values            <- vector(mode = "double", length = N)
  values_with_sign  <- vector(mode = "double", length = N)
  maximum           <- vector(mode = "double", length = 1)

  if (kernel_method == "nw"){
    result <- .C("psihat_statistic_ij", z_data_i, z_data_j, T, g_t_set_vec, N, sigmahat, maximum, values, values_with_sign)
  } else if (kernel_method == "ll"){
    result <- .C("psihat_statistic_ij_ll", z_data_i, z_data_j, T, g_t_set_vec, N, sigmahat, maximum, values, values_with_sign)
  } else {
    print('Given method is currently not supported')
  }
  
  statistic <- result[[7]]
  
  return(statistic) 
}