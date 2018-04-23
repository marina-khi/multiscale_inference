testing_different_time_trends <- function(N_ts, y_data, alpha, kernel_method, sigmahat_vector_2){
  
  T_data <- nrow(y_data)

  ###########################################
  #Calculating gaussian quantile for T_tempr#
  ###########################################
  
  g_t_set                    <- creating_g_set(T_data, kernel_method)
  
  filename = paste("distribution/distr_for_application_T_", T, "_and_N_ts_", N_ts, "_and_method_", kernel_method, ".RData", sep = "")
  if(!file.exists(filename)) {
    gaussian_statistic_distribution <- replicate(1000, {
      z_matrix <- matrix(0, nrow = T_data, ncol = N_ts)
      for (i in 1:N_ts){
        z_matrix[, i] <- rnorm(T_data, 0, sqrt(sigmahat_vector_2[i]))
      }
      psistar_statistic(z_matrix, N_ts, g_t_set, sigmahat_vector_2, kernel_method)
    })
    save(gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)  
  cat("Gaussian quantile is", gaussian_quantile, "\n")
  
  filename = paste("distribution/pairwise_distr_for_application_T_", T, "_and_method_", kernel_method, ".RData", sep = "")
  if(!file.exists(filename)) {
    pairwise_gaussian_statistic_distribution <- replicate(1000, {
      z_matrix <- matrix(0, nrow = T_data, ncol = 2)
      for (i in 1:2){
        z_matrix[, i] <- rnorm(T_data, 0, sqrt(sigmahat_vector_2[i]))
      }
      psistar_statistic(z_matrix, N_ts, g_t_set, sigmahat_vector_2, kernel_method)
    })
    save(pairwise_gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }
  
  pairwise_gaussian_quantile <- quantile(pairwise_gaussian_statistic_distribution, probs = (1 - alpha), type = 1) 
  cat("Pairwise gaussian quantile is", pairwise_gaussian_quantile, "\n")
  
  #########################################
  #Calculating the statistic for real data#
  #########################################
  
  statistic_matrix <- matrix(0, nrow = N_ts, ncol = N_ts)
  results          <- psihat_statistic(y_data, N_ts, g_t_set, sigmahat_vector_2, kernel_method)
  statistic        <- result[[2]]
  statistic_matrix[lower.tri(statistic_matrix, diag=FALSE)] <- result[[1]]
  
  #And now the testing itself
  if (statistic > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", statistic_new,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", statistic_new,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  }
  
  return(list(matrix_of_statistic, gaussian_quantile, pairwise_gaussian_quantile))
}