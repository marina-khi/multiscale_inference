simulations_size <- function(N_ts, N_rep, different_T, different_alpha, kernel_method){
  a_hat <- 0.267
  sigma <- 0.59
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
        simulated_data <- simulated_data - colMeans(simulated_data)[col(simulated_data)]
        result = psihat_statistic(simulated_data, N_ts, g_t_set_sim, sigmahat_vector_2, kernel_method)[[2]]
        result
      })
    for (alpha in different_alpha){
      gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha)
      size <- sum(simulated_statistic > gaussian_quantile)/N_rep
      size_vec = c(size_vec, size)
      cat("Ratio of rejection under H0 is ", size, "with alpha = ", alpha, "and T = ", T_size, "\n")
    }
  }
  filename = paste0("../Plots/", N_ts, "_stations_sizetable_method_", kernel_method, ".tex")
  creating_matrix_and_texing(size, different_T, different_alpha, filename)
  return(size)
}


simulations_power <- function(N_ts, N_rep, different_T, different_alpha, kernel_method){
  a_hat <- 0.267
  sigma <- 0.59
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
        sigmahat_vector_2_simulated <- c()
        simulated_data[, 1]         <- m + arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
        for (i in 2:N_ts){
          simulated_data[, i]         <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, sigma), n = T_size)
          sigma_i                     <- estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2_simulated <- c(sigmahat_vector_2_simulated, sigma_i * sigma_i)
        }
        simulated_data <- simulated_data - colMeans(simulated_data)[col(simulated_data)]
        result   <- psihat_statistic(simulated_data, N_ts, g_t_set_sim, sigmahat_vector_2_simulated, kernel_method)[[2]]
        result
      })
      for (alpha in different_alpha){
        gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha)
        power <- sum(simulated_statistic > gaussian_quantile)/N_rep
        power_b_vec = c(power_b_vec, power)
        cat("Ratio of rejection under H1 is ", power, "with alpha = ", alpha, "T = ", T_size, "and b =", b, "\n")
      }
    }
    filename = paste0("../Plots/", N_ts, "_stations_powertable_method_", kernel_method, "_with_b_", b*100, ".tex")
    creating_matrix_and_texing(power_b_vec, different_T, different_alpha, filename)
  }
  power_vec <- c(power_vec, power_b_vec)
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
    
    simulated_statistic = replicate(N_rep, {
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
      simulated_data <- simulated_data - colMeans(simulated_data)[col(simulated_data)]
      results        <- psihat_statistic(simulated_data, N_ts, g_t_set_sim, sigmahat_vector_2_simulated, kernel_method)[[1]]
      results
    })
    for (alpha in different_alpha){    
      gaussian_quantile = calculating_gaussian_quantile(T_size, N_ts, g_t_set_sim, kernel_method, alpha)
      groups_vec <- matrix(NA, ncol = N_rep, nrow = N_ts)
      number_of_groups_vec <- c()
      for (i in 1:N_rep){
        statistic_vector <- simulated_statistic[, i]
        statistic_value  <- max(statistic_vector)
        if (statistic_value > gaussian_quantile) {
          statistic_matrix                                          <- matrix(0, N_ts, N_ts)
          statistic_matrix[lower.tri(statistic_matrix, diag=FALSE)] <- statistic_vector
          statistic_matrix                                          <- as.dist(statistic_matrix)
          clustering_result <- hclust(statistic_matrix, method = "complete")
          groups            <- cutree(clustering_result, h = gaussian_quantile)
          number_of_groups  <- max(groups)
        } else {
          number_of_groups <- 1
          groups           <- rep(1, N_ts)
        }
      groups_vec[, i]      <- groups
      number_of_groups_vec <- c(number_of_groups_vec, number_of_groups)
      }
      cat("Number of groups:", number_of_groups_vec, "for T = ", T_size, "and alpha = ", alpha, "\n")
    }
    clustering_results <- rbind(number_of_groups_vec, groups_vec)
    filename2 = paste0("distribution/clustering_results_for_T_", T_size, "_and_alpha_", alpha*100, ".RData")
    save(clustering_results, file = filename2)      
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
      N_rep <- ncol(clustering_results)
      
      correct_specification      <- c(rep(1, N_ts/3) , rep(2, N_ts/3), rep(3, N_ts/3))
      correct_number_of_groups   <- 0
      correctly_specified_groups <- 0
      
      for (i in 1:N_rep){
        if (clustering_results[1, i] == 3) {
          correct_number_of_groups = correct_number_of_groups + 1
          groups132  <- recode(clustering_results[2:(N_ts + 1), i], "2=3;3=2")
          groups213  <- recode(clustering_results[2:(N_ts + 1), i], "1=2;2=1")
          groups231  <- recode(clustering_results[2:(N_ts + 1), i], "1=2;2=3;3=1")
          groups312  <- recode(clustering_results[2:(N_ts + 1), i], "1=3;2=1;3=2")
          groups321  <- recode(clustering_results[2:(N_ts + 1), i], "1=3;3=1")
          difference <- min(sum(!correct_specification == groups132), sum(!correct_specification == groups213), 
                            sum(!correct_specification == groups231), sum(!correct_specification == groups312),
                            sum(!correct_specification == groups321), sum(!correct_specification == clustering_results[2:(N_ts + 1), i]))
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