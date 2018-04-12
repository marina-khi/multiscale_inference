testing_different_time_trends <- function(N, y_data, alpha, kernel_method, sigmahat_vector_2){
  
  T_data <- nrow(y_data)
  
  ###########################################
  #Calculating gaussian quantile for T_tempr#
  ###########################################
  
  g_t_set                    <- creating_g_set(T_data, kernel_method)
  gaussian_quantile          <- calculating_gaussian_quantile(T_data, N, g_t_set, kernel_method, alpha)
  pairwise_gaussian_quantile <- calculating_gaussian_quantile(T_data, 2, g_t_set, kernel_method, alpha)

  #########################################
  #Calculating the statistic for real data#
  #########################################
  
  #Calculating sqrt root of the average long-run variance
  sigmahat_tempr <- sqrt(sum(sigmahat_vector_2)/N)
  
  matrix_of_statistic <- matrix(, nrow = N, ncol = N)
  matrix_of_statistic_new <- matrix(, nrow = N, ncol = N)
  for (i in (1 : (N - 1))){
    for (j in ((i + 1):N)){
      sigmahat_new = sqrt(sigmahat_vector_2[i] + sigmahat_vector_2[j])
      result = psihat_statistic_ij(y_data[[i]], y_data[[j]], g_t_set, sqrt(2) * sigmahat_tempr, kernel_method)
      result_new = psihat_statistic_ij(y_data[[i]], y_data[[j]], g_t_set, sigmahat_new, kernel_method)
      matrix_of_statistic[i, j] = result[[2]]
      matrix_of_statistic_new[i, j] = result_new[[2]]
      cat(matrix_of_statistic[i, j], "with i =", i, "and j = ",j, "\n")
      #g_t_set = result[[1]]
    }
  }
  statistic = max(as.vector(matrix_of_statistic), na.rm=TRUE)
  statistic_new = max(as.vector(matrix_of_statistic_new), na.rm=TRUE)
  
  #And now the testing itself
  if (statistic_new > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", statistic_new,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", statistic_new,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  }
  
  if (statistic > gaussian_quantile) {
    cat("For old sigma estimate: We reject H_0 with probability", alpha, "Psihat_statistic = ", statistic,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  } else {
    cat("For old sigma estimate: We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", statistic,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  }
  
  return(list(matrix_of_statistic, matrix_of_statistic_new, gaussian_quantile, pairwise_gaussian_quantile))
}