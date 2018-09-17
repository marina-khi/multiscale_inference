SiZer_simulations <- function(T_size, a_1, sigma_eta, alpha, N_rep, PDFpath){
  kernel_ind  <- 2

  different_i <- seq(from = 1/T_size, to = 1, by = 1/T_size)
  different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
  
  #L1 = floor(sqrt(T_size))
  #L2 = floor(2 * sqrt(T_size))
  
  
  ############################################################
  #Calculating everything that does not depend on Y for SiZer#
  ############################################################

  gamma = c()
  for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
    gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
  }
  
  #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
  true_var <- gamma[1] / T_size
  for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}
  
  T_star   <- gamma[1]/true_var
  
  SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)  

  #################################################################
  #Calculating everything that does not depend on Y for our method#
  #################################################################
  
  #Gaussian statistic for our own method
  gaussian_statistic_distribution <- replicate(1000, {
    z = rnorm(T_size, 0, 1)
    psistar_statistic_ll(z, SiZer_matrix, kernel_ind, 1)
  })
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)
  
  sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))
  
  #################
  #Simulating size#
  #################
  
  pdf(file=PDFPath)  
  
  size_of_the_test <- replicate(N_rep, {
    y_data_ar_1 <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
    #sigmahat    <- estimating_sigma_for_AR1(y_data_ar_1, L1, L2)[[1]]
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
    c(sum(abs(results_their)), sum(abs(results_our)))
  })
  
  cat("Size of SiZer: ", (rowSums(size_of_the_test != 0)/N_rep)[1], ", size of our method: ", (rowSums(size_of_the_test != 0)/N_rep)[2], "\n")
  
  hist(size_of_the_test[1, size_of_the_test[1, ] != 0], main = paste0("Size of SiZer: ", (rowSums(size_of_the_test != 0)/N_rep)[1]), breaks = seq(0, max(size_of_the_test[1, ])+5, by = 5))
  hist(size_of_the_test[2, size_of_the_test[2, ] != 0], breaks = seq(0, max(size_of_the_test[2, ])+5, by = 5), main = paste0("Size of our method: ", (rowSums(size_of_the_test != 0)/N_rep)[2]))
  
  ##################
  #Simulating power#
  ##################
  
  for (a in c(3.5, 4, 4.5, 5)){
    power_of_the_test <- replicate(N_rep, {
      line_trend  <- numeric(T_size)
      for (i in 1:T_size) {line_trend[i] = (i - 0.5*T_size) * a/T_size}
      y_data_ar_1_with_trend <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta)) + line_trend
      #sigmahat    <- estimating_sigma_for_AR1(y_data_ar_1_with_trend, L1, L2)[[1]]
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
      c(sum(abs(results_their)), sum(abs(results_our)), sum(results_their == 1), sum(results_our == 1))
    })
    
    cat("a = ", a, "\n",
        "SiZer results. Percentage of total rejection: ", (rowSums(power_of_the_test != 0)/N_rep)[1],
        ", percentage of positive rejections: ", (rowSums(power_of_the_test != 0)/N_rep)[3], "\n",
        "Our results. Percentage of total rejection: ", (rowSums(power_of_the_test != 0)/N_rep)[2],
        ", percentage of positive rejections: ", (rowSums(power_of_the_test != 0)/N_rep)[4], "\n")
    hist(power_of_the_test[1, power_of_the_test[1, ] != 0], breaks = seq(0, max(power_of_the_test[1, ])+10, by = 10), main = paste0("Linear trend. Power of SiZer: ", (rowSums(power_of_the_test != 0)/N_rep)[1], ", a = ", a))
    hist(power_of_the_test[2, power_of_the_test[2, ] != 0], breaks = seq(0, max(power_of_the_test[2, ])+10, by = 10), main = paste0("Linear trend. Power of our method: ", (rowSums(power_of_the_test != 0)/N_rep)[2], ", a = ", a))
  }
  
  for (a in c(3.5, 4, 4.5, 5)){
    power_of_the_test <- replicate(N_rep, {
      kink_function = numeric(T_size)
      for (i in 1:T_size) {if (i/T_size < 0.5) {kink_function[i] = 0} else {kink_function[i] = (i - 0.5*T_size)*a/T_size}}
      y_data_ar_1_with_trend <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta)) + kink_function
      #sigmahat    <- estimating_sigma_for_AR1(y_data_ar_1_with_trend, L1, L2)[[1]]
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
      c(sum(abs(results_their)), sum(abs(results_our)), sum(results_their == 1), sum(results_our == 1))
    })
    
    cat("a = ", a, "\n",
        "SiZer results. Percentage of total rejection: ", (rowSums(power_of_the_test != 0)/N_rep)[1],
        ", percentage of positive rejections: ", (rowSums(power_of_the_test != 0)/N_rep)[3], "\n",
        "Our results. Percentage of total rejection: ", (rowSums(power_of_the_test != 0)/N_rep)[2],
        ", percentage of positive rejections: ", (rowSums(power_of_the_test != 0)/N_rep)[4], "\n")
    hist(power_of_the_test[1, power_of_the_test[1, ] != 0], breaks = seq(0, max(power_of_the_test[1, ])+10, by = 10), main = paste0("Kink function. Power of SiZer: ", (rowSums(power_of_the_test != 0)/N_rep)[1], ", a = ", a))
    hist(power_of_the_test[2, power_of_the_test[2, ] != 0], breaks = seq(0, max(power_of_the_test[2, ])+10, by = 10), main = paste0("Kink function. Power of our method: ", (rowSums(power_of_the_test != 0)/N_rep)[2], ", a = ", a))
  }
  
  dev.off()
}