SimulateVariance <- function(a_1_star, sigma_eta_star, T_size, Nsim,
                             pdfname_a_hat, pdfname_lrv,
                             sim.design = 'constant', slope = 1,
                             q = 25, r.bar = 10, M1 = 20, M2= 30,
                             produce_plots ='no'){
  # Compares three different methods of estimating the AR(1) parameters, the coefficient a_1 and the long-run variance in particular.
  #   First method is proposed in our, second method is proposed in Hall and Van Keilegom (2003) and the third method is the oracle estimators
  #   the rowwise version of the multiscale test (T_rw) compared to the power of SiZer test (T_SiZer). The 
  #   taht are computed under the assumption that the error process is observed.
  #   More detailed description of the methods can be found in Section 5.2
  #
  # Args:
  #   a_1_star: true AR(1) coefficient for the error distribution.
  #   sigma_eta_star: true standard deviation of the innovation term in the AR(1) process.
  #   T_size:  Sample size of the time series simulated.
  #   Nsim: number of simulations for power calculations. Default is 1000.
  #   sim.design: Type of the trend function. Can be "constant", "linear", "brokenline", "bump", "blocks" or "sine".
  #     Default is "constant".
  #   slope: Height of the bump signal or the slope of "linear" sim.design. Default is 1.
  #   q: tuning parameter for estimating the long-run variance from Section 4. Default is 25. 
  #
  # Returns:
  #   power: One vector of length 4 with each entry corresponding to the power of each testing procedure. The order is as follows.
  #     T_ms, T_uc, T_rw, T_SiZer.
  #   power_ms:    A vector of length equal to the number of bandwidths analysed, each entry corresponding to the rowwise power 
  #     of the multiscale testing procedure (T_ms) for this bandwidth.
  #   power_uncor: A vector of length equal to the number of bandwidths analysed, each entry corresponding to the rowwise power 
  #     of the uncorrected version of the multiscale testing procedure (T_uc) for this bandwidth.
  #   power_rows:  A vector of length equal to the number of bandwidths analysed, each entry corresponding to the rowwise power 
  #     of the the rowwise version of the multiscale testing procedure (T_rw) for this bandwidth.
  #   power_SiZer: A vector of length equal to the number of bandwidths analysed, each entry corresponding to the rowwise power 
  #     of the SiZer testing procedure (T_SiZer) for this bandwidth.
  #   h.grid.new:  A vector of the bandwidths analysed.
  
  
  lrv_star   <- sigma_eta_star^2/((1 - a_1_star)^2)

  a_hat_hall_vect    <- c() 
  a_hat_vect <- c()
  a_hat_oracle_vect  <- c()

  lrv_hat_hall_vect    <- c() 
  lrv_hat_vect <- c()
  lrv_hat_oracle_vect  <- c()
  
  for (i in 1:Nsim){
    #Simulating the data
    data.simulated <- simulating_data(T_size, a_1_star, sigma_eta_star, sim.design = sim.design, slope.fac = slope)
    
    eps               <- data.simulated$eps
    y_data_with_trend <- data.simulated$data
    
    results_our_method <- AR_lrv(data=y_data_with_trend, q=q, r.bar=r.bar, p=1)
    a_hat_vect         <- c(a_hat_vect, results_our_method$ahat)
    lrv_hat_vect       <- c(lrv_hat_vect, results_our_method$lrv)
    
    T_size                          <- as.integer(T_size)
    M1                              <- as.integer(M1)
    M2                              <- as.integer(M2)
    storage.mode(y_data_with_trend) <- "double"
    
    result_HvK         <- EstimateLrvHvK(y_data_with_trend, T_size, M1, M2)
    a_hat_hall_vect    <- c(a_hat_hall_vect, result_HvK$a_hat)
    lrv_hat_hall       <- (result_HvK$sigma_hat)^2
    lrv_hat_hall_vect  <- c(lrv_hat_hall_vect, lrv_hat_hall)

    res_oracle          <- lm(eps[2:T_size] ~ eps[1:(T_size-1)] - 1)
    a_hat_oracle_vect   <- c(a_hat_oracle_vect, as.vector(res_oracle$coef))
    lrv_hat_oracle_vect <- c(lrv_hat_oracle_vect, mean(res_oracle$residuals^2)/(1-sum(as.vector(res_oracle$coef)))^2)
  }

  if (produce_plots == "yes"){
    PlotHistograms(a_hat_vect, a_hat_hall_vect, a_hat_oracle_vect, a_1_star, pdfname_a_hat,
                                 text1 = expression(widehat(a)), text2 = expression(widehat(a)[HvK]),
                                 text3 = expression(widehat(a)[oracle]),
                                 cut.at.end = TRUE)
        
    PlotHistograms(lrv_hat_vect, lrv_hat_hall_vect, lrv_hat_oracle_vect, lrv_star, pdfname_lrv,
                                 text1 = expression(widehat(sigma)[]^2), text2 = expression(widehat(sigma)[HvK]^2),
                                 text3 = expression(widehat(sigma)[oracle]^2),
                                 cut.at.end = FALSE)
  } else if ((produce_plots == 'selected') && ((a_1_star == -0.95) || (a_1_star == 0.25)) && (q == 25) && (r.bar == 10)){
    PlotHistograms(a_hat_vect, a_hat_hall_vect, a_hat_oracle_vect, a_1_star, pdfname_a_hat,
                                   text1 = expression(widehat(a)), text2 = expression(widehat(a)[HvK]),
                                   text3 = expression(widehat(a)[oracle]),
                                   cut.at.end = TRUE)
      
    PlotHistograms(lrv_hat_vect, lrv_hat_hall_vect, lrv_hat_oracle_vect, lrv_star, pdfname_lrv,
                                   text1 = expression(widehat(sigma)[]^2), text2 = expression(widehat(sigma)[HvK]^2),
                                   text3 = expression(widehat(sigma)[oracle]^2),
                                   cut.at.end = FALSE)
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