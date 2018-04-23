simulations_size <- function(N_ts, N_rep, different_T, different_alpha, kernel_method){
  a_hat <- 0.267
  sigma <- 0.59
  size  <- c()
  for (T_size in different_T){
    simulated_data           <- matrix(NA, nrow = T_size, ncol = N_ts)
    colnames(simulated_data) <- 1:N_ts
    L1                       <- floor(sqrt(T_size))
    L2                       <- floor(2 * sqrt(T_size))
    g_t_set_sim              <- creating_g_set(T_size, kernel_method)

    for (alpha in different_alpha){
      gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha)
      size_of_the_test = replicate(N_rep, {
        sigmahat_vector_2 <- c()
        for (i in 1:N_ts){
          simulated_data[, i] <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
          sigma_i = estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_i * sigma_i)
        }
        result = psihat_statistic(simulated_data, N_ts, g_t_set_sim, sigmahat_vector_2, kernel_method)[[2]]
        if (result > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      size = c(size, sum(size_of_the_test)/N_rep)
      cat("Ratio of rejection under H0 is ", sum(size_of_the_test)/N_rep, "with alpha = ", alpha, "and T = ", T_size, "\n")
    }
  }
  filename = paste0("../Plots/", N_ts, "stations_sizetable_method_", kernel_method, ".tex")
  creating_matrix_and_texing(size, different_T, different_alpha, filename)
  return(size)
}


simulations_power <- function(N_ts, N_rep, different_T, different_alpha, kernel_method){
  a_hat <- 0.267
  sigma <- 0.59
  power <- c()
  
  for (b in c(1.00, 1.25)){
    power_b <- c()
    for (T_size in different_T){
      simulated_data           <- matrix(NA, nrow = T_size, ncol = N_ts)
      colnames(simulated_data) <- 1:N_ts
      L1                       <- floor(sqrt(T_size))
      L2                       <- floor(2 * sqrt(T_size))
      g_t_set_sim              <- creating_g_set(T_size, kernel_method)
      
      m <- numeric(T_size)
      for (j in 1:T_size){m[j] = (j - 0.5*T_size) * (b/T_size)}
      
      for (alpha in different_alpha){
        gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha)
        
        power_of_the_test = replicate(N_rep, {
          sigmahat_vector_2_simulated <- c()
          simulated_data[, 1]         <- m + arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
          for (i in 2:N_ts){
            simulated_data[, i]         <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
            sigma_i                     <- estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
            sigmahat_vector_2_simulated <- c(sigmahat_vector_2_simulated, sigma_i * sigma_i)
          }
          result   <- psihat_statistic(simulated_data, N_ts, g_t_set_sim, sigmahat_vector_2_simulated, kernel_method)[[2]]
          if (result > gaussian_quantile) {d = 1} else {d = 0}
          d
        })
        power_b = c(power_b, sum(power_of_the_test)/N_rep)
        cat("Ratio of rejection under H1 is ", sum(power_of_the_test)/N_rep, "with alpha = ", alpha, "T = ", T_size, "and b = ", b, "\n")
      }
    }
    filename = paste0("../Plots/", N_ts, "_stations_powertable_method_", kernel_method, "_with_b_", b*100, ".tex")
    creating_matrix_and_texing(power_b, different_T, different_alpha, filename)
  }
  power <- c(power, power_b)
  return(power)
}

simulations_clustering <- function(N_ts, N_rep, different_T, different_alpha, kernel_method){
  a_hat <- 0.267
  sigma <- 0.59
  clusters <- c()
 
  for (T_size in different_T){
    simulated_data           <- matrix(NA, nrow = T_size, ncol = N_ts)
    colnames(simulated_data) <- 1:N_ts
    L1                       <- floor(sqrt(T_size))
    L2                       <- floor(2 * sqrt(T_size))
    g_t_set_sim              <- creating_g_set(T_size, kernel_method)
    
    m1 <- numeric(T_size)
    m2 <- numeric(T_size)
    for (j in 1:T_size){
      m1[j] = (j - 0.5*T_size) * (1/T_size)
      m2[j] = (j - 0.5*T_size) * (-1/T_size)
    }
    
    for (alpha in different_alpha){
      gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha)
      clustering_results   = replicate(N_rep, {
        sigmahat_vector_2_simulated <- c()
        for (i in 1:(N_ts/3)){
          simulated_data[, i]         <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
          sigma_i                     <- estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2_simulated <- c(sigmahat_vector_2_simulated, sigma_i * sigma_i)
        }
        for (i in (N_ts/3 + 1):(2 * N_ts/3)){
          simulated_data[, i]         <- m1 + arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
          sigma_i                     <- estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2_simulated <- c(sigmahat_vector_2_simulated, sigma_i * sigma_i)
        }
        for (i in (2 * N_ts/3 + 1):N_ts){
          simulated_data[, i]         <- m2 + arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
          sigma_i                     <- estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2_simulated <- c(sigmahat_vector_2_simulated, sigma_i * sigma_i)
        }
        results          <- psihat_statistic(simulated_data, N_ts, g_t_set_sim, sigmahat_vector_2_simulated, kernel_method)
        statistic_vector <- results[[1]]
        statistic_value  <- results[[2]]
        if (statistic_value > gaussian_quantile) {
          statistic_matrix                                          <- matrix(0, N_ts, N_ts)
          statistic_matrix[lower.tri(statistic_matrix, diag=FALSE)] <- statistic_vector
          statistic_matrix                                          <- as.dist(statistic_matrix)
          clustering_result <- hclust(statistic_matrix, method = "complete")
          groups            <- cutree(clustering_result, h = gaussian_quantile)
          number_of_groups  <- max(groups)
          if (number_of_groups == 3){
            
          }
        } else {
          number_of_groups <- 1
          groups           <- rep(1, N_ts)
        }
        list(number_of_groups, groups)
      })
      filename2 = paste0("distribution/clustering_results_for_T_", T_size, "_and_alpha_", alpha*100, ".RData")
      save(clustering_results, file = filename2)      
    }
  }
  return(NULL)
}

clustering_analysis <- function(N_ts, different_T, different_alpha){
  correct_number_of_groups_vec   <- c()
  correctly_specified_groups_vec <- c()
  for (T_size in different_T){
    for (alpha in different_alpha){
      
      filename2 = paste0("distribution/clustering_results_for_T_", T_size, "_and_alpha_", alpha*100, ".RData")
      load(file = filename2)
      N_rep <- length(clustering_results)/2 #Because for each repetition we are storing 2 values in clustering_results
      
      numbers_of_groups <- c()
      groups <- matrix(NA, ncol = N_rep, nrow = N_ts)
      
      for (i in 1:N_rep){
        numbers_of_groups[i] <- clustering_results[[2*i - 1]][1]
        groups[, i] <- clustering_results[[2*i]]
      }
      
      correct_specification      <- c(rep(1, N_ts/3) , rep(2, N_ts/3), rep(3, N_ts/3))
      correct_number_of_groups   <- 0
      correctly_specified_groups <- 0
      
      for (i in 1:N_rep){
        if (numbers_of_groups[i] == 3) {
          correct_number_of_groups = correct_number_of_groups + 1
          groups132 <- recode(groups[, i], "2=3;3=2")
          groups213 <- recode(groups[, i], "1=2;2=1")
          groups231 <- recode(groups[, i], "1=2;2=3;3=1")
          groups312 <- recode(groups[, i], "1=3;2=1;3=2")
          groups321 <- recode(groups[, i], "1=3;3=1")
          difference <- min(sum(!correct_specification == groups132), sum(!correct_specification == groups213), 
                            sum(!correct_specification == groups231), sum(!correct_specification == groups312),
                            sum(!correct_specification == groups321), sum(!correct_specification == groups[, i]))
          if (difference == 0){
            correctly_specified_groups = correctly_specified_groups + 1  
          }
        }
      }
      correct_number_of_groups_vec   <- c(correct_number_of_groups_vec, correct_number_of_groups/N_rep)
      correctly_specified_groups_vec <- c(correctly_specified_groups_vec, correctly_specified_groups/N_rep)
      cat("Percentage of detecting true number of clusters", correct_number_of_groups/N_rep, "with alpha = ", alpha, "T = ", T_size, "\n")
      cat("Percentage of detecting true clustering", correctly_specified_groups/N_rep, "with alpha = ", alpha, "T = ", T_size, "\n")
    }
  }
  filename = paste0("../Plots/", N_ts, "_stations_number_of_groups_method_", kernel_method, ".tex")
  creating_matrix_and_texing(correct_number_of_groups_vec, different_T, different_alpha, filename)
  filename2 = paste0("../Plots/", N_ts, "_stations_groups_method_", kernel_method, ".tex")
  creating_matrix_and_texing(correctly_specified_groups_vec, different_T, different_alpha, filename2)
}