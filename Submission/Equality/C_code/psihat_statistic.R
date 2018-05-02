psihat_statistic <- function(y_data, N_ts, g_t_set, sigmahat_vector_2, kernel_method){
  #Wrapper of a C function from psihat_statistic.C
  #Function that calculates the multiscale statistic \hat{\Psi}_T.
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T        <- as.integer(nrow(y_data))
  N        <- as.integer(nrow(g_t_set))
  
  storage.mode(sigmahat_vector_2)    <- "double"
  storage.mode(N_ts)                 <- "integer"
  
  y_data_vec        <- unlist(y_data)
  g_t_set_vec       <- unlist(g_t_set[c('u', 'h', 'lambda')])
  statistic_vector  <- vector(mode = "double", length = N_ts * (N_ts - 1) / 2)
  statistic_result  <- vector(mode = "double", length = 1)
  
  if (kernel_method == "nw"){
    result <- .C("psihat_statistic_nw", y_data_vec, T, g_t_set_vec, N, N_ts, sigmahat_vector_2, statistic_vector, statistic_result)
  } else if (kernel_method == "ll"){
    result <- .C("psihat_statistic_ll", y_data_vec, T, g_t_set_vec, N, N_ts, sigmahat_vector_2, statistic_vector, statistic_result)
  } else {
    print('Given method is currently not supported')
  }

  statistic       <- result[[7]]
  statistic_value <- result[[8]]
  return(list(statistic, statistic_value)) 
}

psistar_statistic <- function(z_data, N_ts, g_t_set, sigmahat_vector_2, kernel_method){
  #Wrapper of a C function from psihat_statistic.C
  #Function that calculates the multiscale statistic \Psi*star_T.
  #The only difference with the previous function is in the return values.
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T        <- as.integer(nrow(z_data))
  N        <- as.integer(nrow(g_t_set))
  
  storage.mode(sigmahat_vector_2)    <- "double"
  storage.mode(N_ts)                 <- "integer"
  
  z_data_vec        <- unlist(z_data)
  g_t_set_vec       <- unlist(g_t_set[c('u', 'h', 'lambda')])
  statistic_vector  <- vector(mode = "double", length = N_ts * (N_ts - 1) / 2)
  statistic_result  <- vector(mode = "double", length = 1)

  if (kernel_method == "nw"){
    result <- .C("psihat_statistic_nw", z_data_vec, T, g_t_set_vec, N, N_ts, sigmahat_vector_2, statistic_vector, statistic_result)
  } else if (kernel_method == "ll"){
    result <- .C("psihat_statistic_ll", z_data_vec, T, g_t_set_vec, N, N_ts, sigmahat_vector_2, statistic_vector, statistic_result)
  } else {
    print('Given method is currently not supported')
  }
  
  statistic_value <- result[[8]]
  return(statistic_value) 
}





