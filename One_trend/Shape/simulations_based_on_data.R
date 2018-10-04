simulations_general <- function(N, different_T, different_alpha, different_slopes,
                                a_hat, sigma_eta, order, test_problem, kernel_t, L1, K1, K2){

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
  size_ar1 <- c()
  for (alpha in different_alpha){
    for (T in different_T){
      g_t_set =creating_g_set(T, kernel_t)

      #Calculating gaussian quantiles for given T and alpha
      gaussian_quantile = quantile_function(T, g_t_set, test_problem, kernel_ind, alpha)

      #Replicating test procedure N times
      size_of_the_test_ar_1 = replicate(N, {
        y_data_ar_1 <- arima.sim(model = list(ar = a_hat), n = T, innov = rnorm(T, 0, sigma_eta))
        sigmahat    <- estimating_variance_new(y_data_ar_1, L1, L1, order, K1, K2)[[1]]
        result_notrend_ar1 = statistic_function(y_data_ar_1, g_t_set, kernel_ind, sigmahat)[[2]]
        if (result_notrend_ar1 > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      size_ar1 <- c(size_ar1, sum(size_of_the_test_ar_1)/N)
      cat("a_1 = ", a_hat, ". Ratio of rejection under H0 is ", sum(size_of_the_test_ar_1)/N, "with T = ", T, "and alpha = ", alpha, "\n")
    }
  }

  #Calculating the power of the test#
  power_ar1 = c()  
  
  #This is for partly linear function + AR(1) case
  for (alpha in different_alpha){
    for (slope in different_slopes){
      for (T in different_T){
        g_t_set = creating_g_set(T, kernel_t)
        m = numeric(T)
        for (i in 1:T) {if (i/T < 0.6) {m[i] = 0} else {m[i] = (i - 0.6*T)*slope/T}}
        
        #Calculating gaussian quantiles for given T and alpha
        gaussian_quantile = quantile_function(T, g_t_set, test_problem, kernel_ind, alpha)

        #Replicating test procedure N times
        size_of_the_test_with_trend = replicate(N, {
          #Adding a function that is 0 on the first half and linear on the second half
          y_data_ar_1_with_trend = m + arima.sim(model = list(ar = a_hat), innov = rnorm(T, 0, sigma_eta), n = T)
          sigmahat <- estimating_variance_new(y_data_ar_1_with_trend, L1, L1, order, K1, K2)[[1]]
          result_with_trend = statistic_function(y_data_ar_1_with_trend, g_t_set, kernel_ind, sigmahat)[[2]]
          if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
          d
        })
        power_ar1 = c(power_ar1, sum(size_of_the_test_with_trend)/N)
        cat("a_1 = ", a_hat, ". Ratio of rejection under H1 is", sum(size_of_the_test_with_trend)/N, "with slope =", slope, ", T =", T, "and alpha =", alpha, "\n")
      }
    }
  }
  return(list(size_ar1, power_ar1))
}