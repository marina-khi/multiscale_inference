testing_different_time_trends <- function(N_ts, y_data, month_column, alpha, kernel_method, sigmahat_vector_2){
  
  T_data <- nrow(y_data)

  ###########################################
  #Calculating gaussian quantile for T_tempr#
  ###########################################
  
  g_t_set <- creating_g_set(T_data, kernel_method)
  
  filename = paste("distribution/distr_for_application_T_", T_data, "_and_N_ts_", N_ts, "_and_method_", kernel_method, ".RData", sep = "")
  if(!file.exists(filename)) {
    gaussian_statistic_distribution <- replicate(1000, {
      z_matrix <- matrix(, nrow = T_data, ncol = N_ts)
      for (i in 1:N_ts){
        z_matrix[, i] <- rnorm(T_data, 0, sqrt(sigmahat_vector_2[i]))
      }
      z_matrix <- cbind(z_matrix, month_column) #Adding auxiliary month column to "deseasonalize" simulated data
      z_matrix[1:N_ts] <- lapply(z_matrix[1:N_ts], function(x) x - ave(x, z_matrix[N_ts + 1], FUN=mean))
      z_matrix <- z_matrix[-(N_ts + 1)]
      psistar_statistic(z_matrix, N_ts, g_t_set, sigmahat_vector_2, kernel_method)
    })
    save(gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)  
  cat("Gaussian quantile is", gaussian_quantile, "\n")
  

  #########################################
  #Calculating the statistic for real data#
  #########################################
  
  statistic <- psihat_statistic(y_data, N_ts, g_t_set, sigmahat_vector_2, kernel_method)

  #And now the testing itself
  if (statistic[[2]] > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", statistic[[2]],
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", statistic[[2]],
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  }
  
  return(list(statistic[[2]], gaussian_quantile))
}