estimating_variance <- function(a_1_star, sigma_eta_star, T_size, p, slope, sim.design, N_rep, 
                                pdfname_a_hat, pdfname_lrv, q, K2, M1, M2, produce_plots){
  lrv_star   <- sigma_eta_star^2/((1 - a_1_star)^2)

  a_hat_hall_vect    <- c() 
  a_hat_vect <- c()
  a_hat_oracle_vect  <- c()
  
  sigma_eta_hat_hall_vect    <- c() 
  sigma_eta_hat_vect <- c()
  sigma_eta_hat_oracle_vect  <- c()
  
  lrv_hat_hall_vect    <- c() 
  lrv_hat_vect <- c()
  lrv_hat_oracle_vect  <- c()
  
  for (i in 1:N_rep){
    data.simulated <- simulating_data(T_size, a_1_star, sigma_eta_star, sim.design = sim.design, slope.fac = slope)
    
    eps               <- data.simulated$eps
    y_data_with_trend <- data.simulated$data
    
    results_our_method <- AR_lrv(data=y_data_with_trend, q=q, r.bar=K2, p=p)
    a_hat_vect         <- c(a_hat_vect, results_our_method$ahat)
    sigma_eta_hat_vect <- c(sigma_eta_hat_vect, sqrt(results_our_method$vareta))
    lrv_hat_vect       <- c(lrv_hat_vect, results_our_method$lrv)
    
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
    plotting_variance_histograms(a_hat_vect, a_hat_hall_vect, a_hat_oracle_vect, a_1_star, pdfname_a_hat,
                                 text1 = expression(widehat(a)), text2 = expression(widehat(a)[HvK]),
                                 text3 = expression(widehat(a)[oracle]),
                                 cutting_at_end = 'yes')
        
    plotting_variance_histograms(lrv_hat_vect, lrv_hat_hall_vect, lrv_hat_oracle_vect, lrv_star, pdfname_lrv,
                                 text1 = expression(widehat(sigma)[]^2), text2 = expression(widehat(sigma)[HvK]^2),
                                 text3 = expression(widehat(sigma)[oracle]^2),
                                 cutting_at_end = 'no')
  } else if ((produce_plots == 'selected') && ((a_1_star == -0.95) || (a_1_star == 0.25)) && (q == 25) && (K2 == 10)){
      plotting_variance_histograms(a_hat_vect, a_hat_hall_vect, a_hat_oracle_vect, a_1_star, pdfname_a_hat,
                                   text1 = expression(widehat(a)), text2 = expression(widehat(a)[HvK]),
                                   text3 = expression(widehat(a)[oracle]),
                                   cutting_at_end = 'yes')
      
      plotting_variance_histograms(lrv_hat_vect, lrv_hat_hall_vect, lrv_hat_oracle_vect, lrv_star, pdfname_lrv,
                                   text1 = expression(widehat(sigma)[]^2), text2 = expression(widehat(sigma)[HvK]^2),
                                   text3 = expression(widehat(sigma)[oracle]^2),
                                   cutting_at_end = 'no')
  }
  
  # compute MSE values
  a_mse       <- mean((a_hat_vect-a_1_star)^2)
  a_mse_HvK    <- mean((a_hat_hall_vect-a_1_star)^2)
  a_mse_oracle <- mean((a_hat_oracle_vect-a_1_star)^2)
  
  lrv_mse       <- mean((lrv_hat_vect-lrv_star)^2)
  lrv_mse_HvK    <- mean((lrv_hat_hall_vect-lrv_star)^2)
  lrv_mse_oracle <- mean((lrv_hat_oracle_vect-lrv_star)^2)
  return(list(a_mse, a_mse_HvK, a_mse_oracle, lrv_mse, lrv_mse_HvK, lrv_mse_oracle))
}