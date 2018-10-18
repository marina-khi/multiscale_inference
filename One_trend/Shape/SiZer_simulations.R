SiZer_simulations_power <- function(a_1, sigma_eta, N_rep, slopes, different_alpha, different_T){
  set.seed(1) #For reproducibility  
  
  kernel_ind <- 2
  sigmahat   <- sqrt(sigma_eta^2/((1 - a_1)^2))  
  
  power_vect <- c()

  for (slope in slopes){  
    for (T_size in different_T){
      different_i <- seq(from = 5/T_size, to = 1, by = 5/T_size)
      different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
      gamma = c()
      for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
        gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
      }
    
      #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
      true_var <- gamma[1] / T_size
      for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}
    
      T_star   <- gamma[1]/true_var

      for (alpha in different_alpha){

        #THIS PART IS ONLY FOR FAST CALCULATIONS OF THE MATRIX. IF YOU CHANGE different_i OR different_h, YOU CAN'T USE IT!!!
        filename_aux = paste0("Shape/distribution/SiZer_matrix_T_", T_size, "_a_1_", a_1, "_sigma_eta_", sigma_eta, "_alpha_", alpha*100, ".RData")
        if(!file.exists(filename_aux)) {
          SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)
          save(SiZer_matrix, file = filename_aux)
        } else {
          load(filename_aux)
        }
        #USE THIS INSTEAD
        #SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)  
        
        gaussian_quantile <- calculating_gaussian_quantile_ll(T_size, SiZer_matrix, "comparison", kernel_ind, alpha)

        power_of_the_test <- replicate(N_rep, {
          line_trend  <- numeric(T_size)
          for (i in 1:T_size) {line_trend[i] = (i - 0.5*T_size) * slope/T_size}
          y_data_ar_1_with_trend <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta)) + line_trend

          g_t_set     <- psihat_statistic_ll(y_data_ar_1_with_trend, SiZer_matrix, kernel_ind, sigmahat)[[1]]

          results_our   <- c()
          results_their <- c()

          for (row in 1:nrow(g_t_set)){
            i              = g_t_set[row, 'u']
            h              = g_t_set[row, 'h']
            q_h            = g_t_set[row, 'q_h']
            sd_m_hat_prime = g_t_set[row, 'sd']
            
            XtWX_inverse_XtW   = g_t_set$XtWX_inv_XtW[[row]]
            m_hat_prime <- (XtWX_inverse_XtW %*% y_data_ar_1_with_trend)[2]
            
            if (m_hat_prime - q_h * sd_m_hat_prime > 0){
              results_their = c(results_their, 1)
            } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
              results_their = c(results_their, -1)
            } else {
              results_their = c(results_their, 0)
            }

            if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
              results_our = c(results_our, 1)
            } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
              results_our = c(results_our, -1)
            } else {
              results_our = c(results_our, 0)
            }
          }
          c(sum(abs(results_our)), sum(abs(results_their)))
        })

        cat("a_1 = ", a_1, ", slope = ", slope, ", T = ", T_size, ", alpha = ", alpha, "\n",
            "SiZer results. Percentage of total rejection: ", (rowSums(power_of_the_test != 0)/N_rep)[2],
            "Our results. Percentage of total rejection: ", (rowSums(power_of_the_test != 0)/N_rep)[1], "\n")
        power_vect <- c(power_vect, (rowSums(power_of_the_test != 0)/N_rep)[1], (rowSums(power_of_the_test != 0)/N_rep)[2])
      }
    }
  }
  return(power_vect)
}

SiZer_simulations_size <- function(a_1, sigma_eta, N_rep, different_alpha, different_T){
  set.seed(1)
  
  kernel_ind <- 2
  sigmahat   <- sqrt(sigma_eta^2/((1 - a_1)^2))  
  
  size  <- c()

  for (T_size in different_T){
    different_i <- seq(from = 5/T_size, to = 1, by = 5/T_size)
    different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
    gamma = c()
    for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
      gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
    }
    
    #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
    true_var <- gamma[1] / T_size
    for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}
    
    T_star   <- gamma[1]/true_var
    
    for (alpha in different_alpha){
      #THIS PART IS ONLY FOR FAST CALCULATIONS OF THE MATRIX. IF YOU CHANGE different_i OR different_h, YOU CAN'T USE IT!!!
      filename_aux = paste0("Shape/distribution/SiZer_matrix_T_", T_size, "_a_1_", a_1, "_sigma_eta_", sigma_eta, "_alpha_", alpha*100, ".RData")
      if(!file.exists(filename_aux)) {
        SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)
        save(SiZer_matrix, file = filename_aux)
      } else {
        load(filename_aux)
      }
      #USE THIS INSTEAD
      #SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)  
      
      gaussian_quantile <- calculating_gaussian_quantile_ll(T_size, SiZer_matrix, "comparison", kernel_ind, alpha)

      #Simulating size
      size_of_the_test <- replicate(N_rep, {
        y_data_ar_1 <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
        g_t_set     <- psihat_statistic_ll(y_data_ar_1, SiZer_matrix, kernel_ind, sigmahat)[[1]]

        results_our   <- c()
        results_their <- c()

        for (row in 1:nrow(g_t_set)){
          i                = g_t_set[row, 'u']
          h                = g_t_set[row, 'h']
          q_h              = g_t_set[row, 'q_h']
          sd_m_hat_prime   = g_t_set[row, 'sd']
          XtWX_inverse_XtW = g_t_set$XtWX_inv_XtW[[row]]

          m_hat_prime <- (XtWX_inverse_XtW %*% y_data_ar_1)[2]

          if (m_hat_prime - q_h * sd_m_hat_prime > 0){
            results_their = c(results_their, 1)
          } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
            results_their = c(results_their, -1)
          } else {
            results_their = c(results_their, 0)
          }
          if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
            results_our = c(results_our, 1)
          } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
            results_our = c(results_our, -1)
          } else {
            results_our = c(results_our, 0)
          }
        }
        c(sum(abs(results_our)), sum(abs(results_their)))
      })

      cat("a_1 = ", a_1, ", T = ", T_size, ", alpha = ", alpha, "\n",
          "Size of SiZer: ", (rowSums(size_of_the_test != 0)/N_rep)[2],
          ", size of our method: ", (rowSums(size_of_the_test != 0)/N_rep)[1], "\n")
      size <- c(size, (rowSums(size_of_the_test != 0)/N_rep)[1], (rowSums(size_of_the_test != 0)/N_rep)[2])    }
  }
  return(size)
}
