simulations_size <- function(a_hat, sigma, N_ts, N_rep, different_T, different_alpha, kernel_method){
  size_vec  <- c()

  for (T_size in different_T){
    simulated_data           <- matrix(NA, nrow = T_size, ncol = N_ts)
    colnames(simulated_data) <- 1:N_ts
    L1                       <- floor(sqrt(T_size))
    L2                       <- floor(2 * sqrt(T_size))
    g_t_set_sim              <- creating_g_set(T_size, kernel_method)

    simulated_statistic = replicate(N_rep, {
        sigmahat_vector_2 <- c()
        for (i in 1:N_ts){
          simulated_data[, i] <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
          sigma_i = estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_i * sigma_i)
        }
        sigmahat_averaged_2 <- rep(sum(sigmahat_vector_2)/N_ts, N_ts)
        simulated_data <- simulated_data - colMeans(simulated_data)[col(simulated_data)]
        result = psihat_statistic(simulated_data, N_ts, g_t_set_sim, sigmahat_averaged_2, kernel_method)[[2]]
        result
      })
    for (alpha in different_alpha){
      filename_distr = paste("Equality/distribution/distr_T_", T_size, "_and_N_", N_ts, "_and_method_", kernel_method, ".RData", sep = "")
      gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha, filename_distr)
      size <- sum(simulated_statistic > gaussian_quantile)/N_rep
      size_vec = c(size_vec, size)
      cat("Ratio of rejection under H0 is ", size, "with alpha = ", alpha, "and T = ", T_size, "\n")
    }
  }
  filename = paste0("Plots/", N_ts, "_stations_sizetable_method_", kernel_method, ".tex")
  creating_matrix_and_texing(size_vec, different_T, different_alpha, filename)
  return(size_vec)
}


simulations_power <- function(a_hat, sigma, N_ts, N_rep, different_T, different_alpha, kernel_method){
  power_vec <- c()
  
  for (b in c(0.75, 1.00, 1.25)){
    power_b_vec <- c()
    for (T_size in different_T){
      simulated_data           <- matrix(NA, nrow = T_size, ncol = N_ts)
      colnames(simulated_data) <- 1:N_ts
      L1                       <- floor(sqrt(T_size))
      L2                       <- floor(2 * sqrt(T_size))
      g_t_set_sim              <- creating_g_set(T_size, kernel_method)
      
      m <- numeric(T_size)
      for (j in 1:T_size){m[j] = (j - 0.5*T_size) * (b/T_size)}
      
      simulated_statistic = replicate(N_rep, {
        sigmahat_vector_2   <- c()
        simulated_data[, 1] <- m + arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
        for (i in 2:N_ts){
          simulated_data[, i] <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
          sigma_i             <- estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2   <- c(sigmahat_vector_2, sigma_i * sigma_i)
        }
        sigmahat_averaged_2 <- rep(sum(sigmahat_vector_2)/N_ts, N_ts)
        simulated_data      <- simulated_data - colMeans(simulated_data)[col(simulated_data)]
        result              <- psihat_statistic(simulated_data, N_ts, g_t_set_sim, sigmahat_averaged_2, kernel_method)[[2]]
        result
      })
      for (alpha in different_alpha){
        filename_distr = paste("Equality/distribution/distr_T_", T_size, "_and_N_", N_ts, "_and_method_", kernel_method, ".RData", sep = "")
        gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha, filename_distr)
        power <- sum(simulated_statistic > gaussian_quantile)/N_rep
        power_b_vec = c(power_b_vec, power)
        cat("Ratio of rejection under H1 is ", power, "with alpha = ", alpha, "T = ", T_size, "and b =", b, "\n")
      }
    }
    filename = paste0("Plots/", N_ts, "_stations_powertable_method_", kernel_method, "_with_b_", b*100, ".tex")
    creating_matrix_and_texing(power_b_vec, different_T, different_alpha, filename)
  }
  power_vec <- c(power_vec, power_b_vec)
  return(power)
}