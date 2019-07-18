SimulateVariance <- function(a_1_star, sigma_eta_star, T_size, Nsim,
                             pdfname_a_hat, pdfname_lrv,
                             sim.design = 'constant', slope = 1,
                             q = 25, r.bar = 10, M1 = 20, M2= 30,
                             produce.plots ='no'){
  # Compares three different methods of estimating the AR(1) parameters, in particular the coefficient a_1 and the long-run variance.
  #   First method is proposed in our paper, second method is proposed in Hall and Van Keilegom (2003) and the third method is the oracle
  #   one that is calculated as if the errors are directly observed. More detailed description of the methods can be found in Section 5.2
  #
  # Args:
  #   a_1_star: True AR(1) coefficient for the error distribution.
  #   sigma_eta_star: True standard deviation of the innovation term in the AR(1) process.
  #   T_size:  Sample size of the time series simulated.
  #   Nsim: number of simulations for power calculations. Default is 1000.
  #   pdfname_a_hat: Name (and path) of the pdf file which will contain the histograms for estimators of a_1.
  #   pdfname_lrv: Name (and path) of the pdf file which will contain the histograms for estimators of the long-run variance.
  #   sim.design: Type of the trend function. Can be "constant", "linear", "brokenline", "bump", "blocks" or "sine".
  #     Default is "constant".
  #   slope: Height of the bump signal or the slope of "linear" sim.design. Default is 1.
  #   q, r.bar: Tuning parameters for estimating the long-run variance from Section 4. Default are 25 and 10. 
  #   M1, M2: Tuning parameters for estimating the long-run variance from Hall and Van Keilegom (2003). Default are 20 and 30.
  #   produce.plots: If "yes", then the function produces histograms of the estimators for all the specifications used and 
  #     stores them in the corresponding files. If 'selected', then the function produces histograms of the estimators only for
  #     the following specifications: (a_1 = -0.95 or a_1 = 0.25) and q = 25 and r.bar = 10. If 'no', no histograms are produced.
  #     Default is "no".
  #   
  # Returns:
  #   a_mse: Value of the mean squared error for the estimator \widehat{a} of a_1 for our method.
  #   a_mse_HvK: Value of the mean squared error for the estimator \widehat{a}_{HvK} of a_1
  #     for the method from Hall and Van Keilegom (2003).
  #   a_mse_oracle: Value of the mean squared error for the oracle estimator \widehat{a}_{oracle} of a_1.
  #   lrv_mse: Value of the mean squared error for the estimator \widehat{\sigma}^2 of the long-run variance for our method.
  #   lrv_mse_HvK: Value of the mean squared error for the estimator \widehat{\sigma}^2_{HvK} of the long-run variance
  #     for the method from Hall and Van Keilegom (2003).
  #   lrv_mse_oracle: Value of the mean squared error for the oracle estimator \widehat{\sigma}^2_{oracle} of the long-run variance.
 
  #Load necessary functions  
  source("functions/long_run_variance.r")
  source("functions/sim.r")
  sourceCpp("functions/EstimateLrvHvK.cpp")
  
  lrv_star   <- sigma_eta_star^2/((1 - a_1_star)^2)

  a_hat_hall_vect    <- c()
  a_hat_vect <- c()
  a_hat_oracle_vect  <- c()

  lrv_hat_hall_vect    <- c() 
  lrv_hat_vect <- c()
  lrv_hat_oracle_vect  <- c()
  
  for (i in 1:Nsim){
    #Simulate the data
    data.simulated <- simulating_data(T_size, a_1_star, sigma_eta_star, sim.design = sim.design, slope.fac = slope)
    
    eps               <- data.simulated$eps
    y_data_with_trend <- data.simulated$data

    #Estimate with our method
    results_our_method <- AR_lrv(data=y_data_with_trend, q=q, r.bar=r.bar, p=1)
    a_hat_vect         <- c(a_hat_vect, results_our_method$ahat)
    lrv_hat_vect       <- c(lrv_hat_vect, results_our_method$lrv)
    
    #Estimate with the method proposed by Hall and Van Keilegom (2003)
    T_size                          <- as.integer(T_size)
    M1                              <- as.integer(M1)
    M2                              <- as.integer(M2)
    storage.mode(y_data_with_trend) <- "double"
    
    result_HvK         <- EstimateLrvHvK(y_data_with_trend, T_size, M1, M2)
    a_hat_hall_vect    <- c(a_hat_hall_vect, result_HvK$a_hat)
    sigma_eta_hat_hall <- calculating_sigma_eta(y_data_with_trend, result_HvK$a_hat, p = 1)
    lrv_hat_hall       <- sigma_eta_hat_hall^2/((1 - sum(result_HvK$a_hat))^2)
    lrv_hat_hall_vect  <- c(lrv_hat_hall_vect, lrv_hat_hall)

    #Estimate as if the errors are observed - oracle estimator

    res_oracle          <- lm(eps[2:T_size] ~ eps[1:(T_size-1)] - 1)
    a_hat_oracle_vect   <- c(a_hat_oracle_vect, as.vector(res_oracle$coef))
    lrv_hat_oracle_vect <- c(lrv_hat_oracle_vect, mean(res_oracle$residuals^2)/(1-sum(as.vector(res_oracle$coef)))^2)
  }

  if (produce.plots == "yes"){
    PlotHistograms(a_hat_vect, a_hat_hall_vect, a_hat_oracle_vect, a_1_star, pdfname_a_hat,
                                 text1 = expression(widehat(a)), text2 = expression(widehat(a)[HvK]),
                                 text3 = expression(widehat(a)[oracle]),
                                 cut.at.end = TRUE)
        
    PlotHistograms(lrv_hat_vect, lrv_hat_hall_vect, lrv_hat_oracle_vect, lrv_star, pdfname_lrv,
                                 text1 = expression(widehat(sigma)[]^2), text2 = expression(widehat(sigma)[HvK]^2),
                                 text3 = expression(widehat(sigma)[oracle]^2),
                                 cut.at.end = FALSE)
  } else if ((produce.plots == 'selected') && ((a_1_star == -0.95) || (a_1_star == 0.25)) && (q == 25) && (r.bar == 10)){
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