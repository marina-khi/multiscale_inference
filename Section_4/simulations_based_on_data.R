simulations_based_on_data <- function(N_ts, N_rep, y_data, different_T, different_alpha, kernel_method){

  ##############################
  #Defining necessary constants#
  ##############################
  
  a_hat <- 0.267 #This I need to redo!!!!
  
  #####################
  #Simulating the size#
  #####################
  
  size <- c()
  
  for (T_size in different_T){
    simulated_data           <- matrix(NA, nrow = T_size, ncol = N_ts)
    colnames(simulated_data) <- 1:N_ts
    L1                       <- floor(sqrt(T_size))
    L2                       <- floor(2 * sqrt(T_size))
    g_t_set_sim              <- creating_g_set(T_size, kernel_method)
  
    for (alpha in different_alpha){
      #Calculating gaussian quantiles for given T and alpha
      gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha)
      
      size_of_the_test = replicate(N_rep, {
        sigmahat_vector_2_simulated <- c()
        for (i in 1:N_ts){
          simulated_data[, i] <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, 0.59), n = T_size)
          sigma_i = estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2_simulated <- c(sigmahat_vector_2_simulated, sigma_i * sigma_i)
        }
        sigmahat <- sqrt(sum(sigmahat_vector_2_simulated)/N_ts)
        result = psihat_statistic(simulated_data, N_ts, g_t_set_sim, sqrt(2) * sigmahat, kernel_method)[[2]]
        if (result > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      size = c(size, sum(size_of_the_test)/N_rep)
      cat("Ratio of rejection under H0 is ", sum(size_of_the_test)/N_rep, "with alpha = ", alpha, "and T = ", T_size, "\n")
    }
  }
  
  ######################
  #Simulating the power#
  ######################
  
  #Simulating the power
  
  power <- c()
  
  for (b in c(1.25, 1.875, 2.5)){
    for (T_size in different_T){
      simulated_data           <- matrix(NA, nrow = T_size, ncol = N_ts)
      colnames(simulated_data) <- 1:N_ts
      L1                       <- floor(sqrt(T_size))
      L2                       <- floor(2 * sqrt(T_size))
      g_t_set_sim              <- creating_g_set(T_size, kernel_method)
      
      m <- numeric(T_size)
      for (j in 1:T_size){if (j/T < 0.6) {m[i] = 0} else {m[j] = (j - 0.6*T)*b}}
      
      for (alpha in different_alpha){
        gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha)
        
        power_of_the_test = replicate(N_rep, {
          sigmahat_vector_2_simulated <- c()
          simulated_data[, 1]         <- m + arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, 0.59), n = T_size)
          for (i in 2:N_ts){
            simulated_data[, i]         <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, 0.59), n = T_size)
            sigma_i                     <- estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
            sigmahat_vector_2_simulated <- c(sigmahat_vector_2_simulated, sigma_i * sigma_i)
          }
          sigmahat <- sqrt(sum(sigmahat_vector_2_simulated)/N_ts)
          result   <- psihat_statistic(simulated_data, N_ts, g_t_set_sim, sqrt(2) * sigmahat, kernel_method)[[2]]
          if (result > gaussian_quantile) {d = 1} else {d = 0}
          d
        })
        power = c(power, sum(power_of_the_test)/N_rep)
        cat("Ratio of rejection under H1 is ", sum(power_of_the_test)/N_rep, "with alpha = ", alpha, "T = ", T_size, "and b = ", b, "\n")
      }
    }
  }
  return(list(size, power))
}