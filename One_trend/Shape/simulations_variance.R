histograms_for_variance_estimators <- function(a_1_star, sigma_eta_star, T_size, p, slope, N_rep,
                                               pdfname_a_hat, pdfname_lrv, q, K1, K2, M1, M2, produce_plots){
  lrv_star   <- sigma_eta_star^2/((1 - a_1_star)^2)
  line_trend <- (1:T_size) * slope/T_size
    
  a_hat_hall_vect    <- c() 
  a_hat_method2_vect <- c()
  a_hat_oracle_vect  <- c()
  
  sigma_eta_hat_hall_vect    <- c() 
  sigma_eta_hat_method2_vect <- c()
  sigma_eta_hat_oracle_vect  <- c()
  
  lrv_hat_hall_vect    <- c() 
  lrv_hat_method2_vect <- c()
  lrv_hat_oracle_vect  <- c()
  
  for (i in 1:N_rep){
    eps               <- arima.sim(model = list(ar = a_1_star), n = T_size, innov = rnorm(T_size, 0, sigma_eta_star))
    y_data_with_trend <- eps + line_trend
    
    results_our_method         <- estimating_variance_new(y_data_with_trend, q, q, p, K1, K2)
    a_hat_method2_vect         <- c(a_hat_method2_vect, results_our_method[[2]])
    sigma_eta_hat_method2_vect <- c(sigma_eta_hat_method2_vect, results_our_method[[3]])
    lrv_hat_method2_vect       <- c(lrv_hat_method2_vect, results_our_method[[1]]^2)
    
    result_HvK         <- estimating_sigma_for_AR1(y_data_with_trend, M1, M2)
    a_hat_hall_vect    <- c(a_hat_hall_vect, result_HvK[[2]])
    sigma_eta_hat_hall <- calculating_sigma_eta(y_data_with_trend, result_HvK[[2]], p)
    lrv_hat_hall       <- sigma_eta_hat_hall^2/((1 - sum(result_HvK[[2]]))^2)
    sigma_eta_hat_hall_vect <- c(sigma_eta_hat_hall_vect, sigma_eta_hat_hall)
    lrv_hat_hall_vect       <- c(lrv_hat_hall_vect, lrv_hat_hall)

    res_oracle                <- lm(eps[2:T_size] ~ eps[1:(T_size-1)] - 1)
    a_hat_oracle_vect         <- c(a_hat_oracle_vect, as.vector(res_oracle$coef))
    sigma_eta_hat_oracle_vect <- c(sigma_eta_hat_oracle_vect, sqrt(mean(res_oracle$residuals^2)))
    lrv_hat_oracle_vect       <- c(lrv_hat_oracle_vect, mean(res_oracle$residuals^2)/(1-sum(as.vector(res_oracle$coef)))^2)
  }

  if (produce_plots == "yes"){
    plotting_variance_histograms(a_hat_method2_vect, a_hat_hall_vect, a_hat_oracle_vect, a_1_star, pdfname_a_hat,
                                 text1 = expression(widehat(a)), text2 = expression(widehat(a)[HvK]),
                                 text3 = expression(widehat(a)[oracle]),
                                 cutting_at_end = 'yes')
        
    plotting_variance_histograms(lrv_hat_method2_vect, lrv_hat_hall_vect, lrv_hat_oracle_vect, lrv_star, pdfname_lrv,
                                 text1 = expression(widehat(sigma)[]^2), text2 = expression(widehat(sigma)[HvK]^2),
                                 text3 = expression(widehat(sigma)[oracle]^2),
                                 cutting_at_end = 'no')
  } else if ((produce_plots == 'selected') && ((a_1_star == -0.95) || (a_1_star == 0.25)) && (q == 25) && (K2 == 10)){
      plotting_variance_histograms(a_hat_method2_vect, a_hat_hall_vect, a_hat_oracle_vect, a_1_star, pdfname_a_hat,
                                   text1 = expression(widehat(a)), text2 = expression(widehat(a)[HvK]),
                                   text3 = expression(widehat(a)[oracle]),
                                   cutting_at_end = 'yes')
      
      plotting_variance_histograms(lrv_hat_method2_vect, lrv_hat_hall_vect, lrv_hat_oracle_vect, lrv_star, pdfname_lrv,
                                   text1 = expression(widehat(sigma)[]^2), text2 = expression(widehat(sigma)[HvK]^2),
                                   text3 = expression(widehat(sigma)[oracle]^2),
                                   cutting_at_end = 'no')
  }
  
  # compute MSE values
  a_mse2       <- mean((a_hat_method2_vect-a_1_star)^2)
  a_mse_HvK    <- mean((a_hat_hall_vect-a_1_star)^2)
  a_mse_oracle <- mean((a_hat_oracle_vect-a_1_star)^2)
  
  lrv_mse2       <- mean((lrv_hat_method2_vect-lrv_star)^2)
  lrv_mse_HvK    <- mean((lrv_hat_hall_vect-lrv_star)^2)
  lrv_mse_oracle <- mean((lrv_hat_oracle_vect-lrv_star)^2)
  return(list(a_mse2, a_mse_HvK, a_mse_oracle, lrv_mse2, lrv_mse_HvK, lrv_mse_oracle))
}