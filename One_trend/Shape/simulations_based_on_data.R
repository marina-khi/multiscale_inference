simulations_general <- function(N, different_T, different_alpha, a_hat, sigma_eta, order, test_problem, kernel_t, filename, L1, L2, K1, K2){
  #Defining necessary parameters
  num_T     <- length(different_T)
  num_alpha <- length(different_alpha)
  
  #Recoding kernel functions and type of kernel estimator 
  if (test_problem == "zero"){
    kernel_ind = 1
  } else if (test_problem == "constant"){
    kernel_ind = 2
  } else {
    print('Given testing problem is currently not supported')
  }
 
  if (kernel_t == "nw"){
    quantile_function = calculating_gaussian_quantile
    statistic_function = psihat_statistic
  } else if (kernel_t == "ll"){
    quantile_function = calculating_gaussian_quantile_ll
    statistic_function = psihat_statistic_ll
  } else {
    print('Given method is currently not supported')
  }

  #Calculating the size of the test

  #We don't have the trend function, the null hypothesis is true and the errors are from AR(1)
  size_ar1 <- c()
  for (T in different_T){

    g_t_set =creating_g_set(T, kernel_t)
    
    for (alpha in different_alpha){
      #Calculating gaussian quantiles for given T and alpha
      gaussian_quantile = quantile_function(T, g_t_set, test_problem, kernel_ind, alpha)
      cat("Gaussian quantile = ", gaussian_quantile, "with T = ", T, "and alpha = ", alpha, "\n")
      
      #Replicating test procedure N times
      size_of_the_test_ar_1 = replicate(N, {
        y_data_ar_1 = arima.sim(model = list(ar = a_hat), n = T, innov = rnorm(T, 0, sigma_eta))
        sigmahat <- estimating_variance_new(y_data_ar_1, L1, L2, order, K1, K2)[[1]]
        result_notrend_ar1 = statistic_function(y_data_ar_1, g_t_set, kernel_ind, sigmahat)[[2]]
        if (result_notrend_ar1 > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      size_ar1 <- c(size_ar1, sum(size_of_the_test_ar_1)/N)
      cat("Ratio of rejection under H0 is ", sum(size_of_the_test_ar_1)/N, "with T = ", T, "and alpha = ", alpha, "\n")
    }
  }
  #Creating a nice-looking matrix for size of the test and writing it to a tex file
  filename0 = paste0(filename, "_sizetable_", kernel_t, "_testing_", test_problem, ".tex")
  creating_matrix_and_texing(size_ar1, different_T, different_alpha, filename0)

  ###################################
  #Calculating the power of the test#
  ###################################
  

  #This is for partly linear function + AR(1) case
  for (slope in c(4.0, 3.5, 3.0)){
    power_ar1 = c()
    for (T in different_T){
      g_t_set = creating_g_set(T, kernel_t)
      grid_points <- seq(from = 1/T, to = 1, length.out = T)
      
      m = numeric(T)
      for (i in 1:T) {m[i] = (i - 0.5*T) * slope/T}
      
      #b = a/(0.4 * T) #Constant needed just for calculation
      #
      #for (i in 1:T) {if (i/T < 0.6) {m[i] = 0} else {m[i] = (i - 0.6*T)*b}}
      #lines(grid_points, m, lty = i)
      
      for (alpha in different_alpha){
        #Calculating gaussian quantiles for given T and alpha
        gaussian_quantile = quantile_function(T, g_t_set, test_problem, kernel_ind, alpha)
        cat("Gaussian quantile = ", gaussian_quantile, "with T = ", T, "and alpha = ", alpha, "\n")
        
        #Replicating test procedure N times
        size_of_the_test_with_trend = replicate(N, {
          #Adding a function that is 0 on the first half and linear on the second half
          y_data_ar_1_with_trend = m + arima.sim(model = list(ar = a_hat), innov = rnorm(T, 0, sigma_eta), n = T)
          sigmahat <- estimating_variance_new(y_data_ar_1_with_trend, L1, L2, order, K1, K2)[[1]]
          result_with_trend = statistic_function(y_data_ar_1_with_trend, g_t_set, kernel_ind, sigmahat)[[2]]
          if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
          d
        })
        power_ar1 = c(power_ar1, sum(size_of_the_test_with_trend)/N)
        cat("Ratio of rejection under H1 with partially linear trend is", sum(size_of_the_test_with_trend)/N, "with slope =", slope, ", T =", T, "and alpha =", alpha, "\n")
      }
    }
    filename0 = paste0(filename, "_powertable_slope_", slope*100, "_", kernel_t, "_testing_", test_problem, ".tex")
    creating_matrix_and_texing(power_ar1, different_T, different_alpha, filename0)
  }
}

simulations_based_on_data <- function(N_rep, y_data, different_alpha, different_slopes, order, test_problem, L1, L2, K1, K2){
  #Defining necessary parameters

  num_alpha <- length(different_alpha)
  T_data    <- length(y_data)

  result             <- estimating_variance_new(y_data, L1, L2, order, K1, K2)
  sigma_hat_data     <- result[[1]]
  a_hat_data         <- result[[2]]
  sigma_eta_hat_data <- result[[3]]
  
  #Recoding kernel functions and type of kernel estimator 
  if (test_problem == "zero"){
    kernel_ind = 1
  } else if (test_problem == "constant"){
    kernel_ind = 2
  } else {
    print('Given testing problem is currently not supported')
  }
  
  g_t_set = creating_g_set(T_data, "ll")
  
  #Calculating size and power of the test#
  
  size_and_power <- c()
  for (slope in different_slopes){

    m = numeric(T_data)
    for (i in 1:T_data) {if (i/T_data < 0.6) {m[i] = 0} else {m[i] = (i - 0.6*T_data)*slope/T_data}}

    for (alpha in different_alpha){
      #Calculating gaussian quantiles for given T and alpha
      gaussian_quantile = calculating_gaussian_quantile_ll(T_data, g_t_set, test_problem, kernel_ind, alpha)
      cat("Gaussian quantile = ", gaussian_quantile, "with T = ", T_data, "and alpha = ", alpha, "\n")
        
      #Replicating test procedure N_rep times
      size_of_the_test_with_trend = replicate(N_rep, {
        #Adding a function that is 0 on the first half and linear on the second half
        y_data_with_trend = m + arima.sim(model = list(ar = a_hat_data), innov = rnorm(T_data, 0, sigma_eta_hat_data), n = T_data)
        sigmahat <- estimating_variance_new(y_data_with_trend, L1, L2, order, K1, K2)[[1]]
        result_with_trend = psihat_statistic_ll(y_data_with_trend, g_t_set, kernel_ind, sigmahat)[[2]]
        if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      size_and_power = c(size_and_power, sum(size_of_the_test_with_trend)/N_rep)
      cat("Ratio of rejection under H1 with partially linear trend is", sum(size_of_the_test_with_trend)/N_rep, "with slope =", slope, ", T =", T_data, "and alpha =", alpha, "\n")
    }
  }
  creating_matrix_and_texing(size_and_power, different_slopes, different_alpha, "Paper/Plots/finite_sample_properties_based_on_application.tex")
}