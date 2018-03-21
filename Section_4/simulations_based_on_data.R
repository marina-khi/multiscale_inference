simulations_based_on_data <- function(N_ts, N_rep, different_T, different_alpha, kernel_method){

  ##############################
  #Defining necessary constants#
  ##############################
  
  a_hat <- 0.5
  
  if (kernel_method == "nw"){
    quantile_function = calculating_gaussian_quantile
    statistic_function = psihat_statistic_ij
    defining_set = creating_g_set
  } else if (kernel_method == "ll"){
    quantile_function = calculating_gaussian_quantile_ll
    statistic_function = psihat_statistic_ij_ll
    defining_set = creating_g_set_ll
  } else {
    print('Given method is currently not supported')
  }
  
  ######################################
  #Simulating the data for the AR1 case#
  ######################################
  
  size_ar1 <- c()
  
  for (T_size in different_T){
    simulated_data <- matrix(NA, nrow = T_size, ncol = N_ts)
    colnames(simulated_data) <- 1:N_ts
    L1 <- floor(sqrt(T_size))
    L2 <- floor(2 * sqrt(T_size))
    g_t_set_sim = defining_set(T_size)
  
    for (alpha in different_alpha){
      #Calculating gaussian quantiles for given T and alpha
      gaussian_quantile = quantile_function(T_size, N_ts, g_t_set_sim, alpha)
      
      size_of_the_test = replicate(N_rep, {
        sigmahat_vector_2 <- c()
        for (i in 1:N_ts){
          simulated_data[, i] <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, 1), n = T_size)
          sigma_i = estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
          sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_i * sigma_i)
        }
        #sigmahat <- sqrt(sum(sigmahat_vector_2)/N)
        
        matrix_of_statistic <- matrix(, nrow = N_ts, ncol = N_ts)
        for (i in (1 : (N - 1))){
          for (j in ((i + 1):N)){
            result = statistic_function(simulated_data[, i], simulated_data[, j], g_t_set_sim, sqrt(sigmahat_vector_2[i] + sigmahat_vector_2[j]))
            matrix_of_statistic[i, j] = result[[2]]
            #      cat(matrix_of_statistic[i, j], "with i =", i, "and j = ",j, "\n")
          }
        }
        result_with_trend = max(as.vector(matrix_of_statistic), na.rm=TRUE)
        if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      size_ar1 = c(size_ar1, sum(size_of_the_test)/N_rep)
      cat(sum(size_of_the_test)/N_rep, "with alpha = ", alpha, "and T = ", T_size, "\n")
    }
  }
  return(size_ar1)
}