psihat_statistic <- function(y_data, g_t_set, kernel_ind, sigmahat){
  #Wrapper of a C function from psihat_statistic.C
  #Function that calculates the multiscale statistic \hat{\Psi}_T.
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # kernel_ind      type of kernel function used - 1 for epanechnikov, 2 for biweight
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T <- as.integer(length(y_data))
  N <- as.integer(nrow(g_t_set))
  sigmahat <- as.double(sigmahat)
  kernel <- as.integer(kernel_ind)
  
  #storage.mode(kernel_ind)  <- "integer"
  storage.mode(sigmahat)    <- "double"
  
  g_t_set_vec       <- unlist(g_t_set[c('u', 'h', 'lambda')])
  values            <- vector(mode = "double",length = N)
  values_with_sign  <- vector(mode = "double",length = N)
  maximum           <- vector(mode = "double", length = 1)
  
  result <- .C("psihat_statistic", y_data, T, g_t_set_vec, N, kernel, sigmahat, maximum, values, values_with_sign)

#  cat("kernel = ", result[[5]], ", T=", result[[2]], ", N=", result[[4]], ", sigmahat = ", result[[6]],
#      "maximum = ", result[[7]], "len of values = ", length(result[[8]]), sep=" ")
  
  g_t_set$values  <- result[[8]]
  g_t_set$values_with_sign <- result[[9]]
  statistic       <- result[[7]]

  return(list(g_t_set, statistic)) 
}

psistar_statistic <- function(y_data, g_t_set, kernel_ind = 1, sigmahat){
  #Wrapper of a C function from psihat_statistic.C
  #Function that calculates the auxiliary statistic \Psi^star_T.
  #The only difference with the previous function is in the return values.
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # kernel_ind      type of kernel function used - 1 for epanechnikov, 2 for biweight
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T <- as.integer(length(y_data))
  N <- as.integer(nrow(g_t_set))
  
  storage.mode(kernel_ind)  <- "integer"
  storage.mode(sigmahat)    <- "double"
  
  g_t_set_vec <- unlist(g_t_set[c('u', 'h', 'lambda')])
  values      <- vector(mode = "double",length = N)
  values_with_sign  <- vector(mode = "double",length = N)
  maximum     <- vector(mode = "double", length = 1)
  
  result    <- .C("psihat_statistic", y_data, T, g_t_set_vec, N, kernel_ind, sigmahat, maximum, values, values_with_sign)
  statistic <- result[[7]]
  
  return(statistic) 
}