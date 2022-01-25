#####################
#Estimating variance#
#####################

sigma <- function(gdp_mat, n_ts, q = 15, r = 10){
  #Order selection
  q_vec <- 10:20
  r_vec <- 10:15
  order_results <- c()
  
  for (j in 1:n_ts){
    criterion_matrix <- expand.grid(q = q_vec, r = r_vec)
    
    criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
    
    for (i in 1:nrow(criterion_matrix)){
      FPE <- c()
      AIC <- c()
      AICC <- c()
      SIC <- c()
      HQ <- c()
      
      different_orders <- (1:9)
      
      for (order in different_orders){
        AR.struc      <- estimate_lrv(data = gdp_mat_augm[, j], q = criterion_matrix$q[[i]],
                                      r_bar = criterion_matrix$r[[i]], p = order)
        sigma_eta_hat <- sqrt(AR.struc$vareta)
        FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
        AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
        AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order +2)/t_len))
        SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
        HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
      }
      criterion_matrix$FPE[[i]]  <- which.min(FPE)
      criterion_matrix$AIC[[i]]  <- which.min(AIC)
      criterion_matrix$AICC[[i]] <- which.min(AICC)
      criterion_matrix$SIC[[i]]  <- which.min(SIC)
      criterion_matrix$HQ[[i]]   <- which.min(HQ)
    }
    maxim <- max(criterion_matrix[, 3:7])
    order_results <- c(order_results, maxim)
    #cat("For the country ", colnames(gdp_mat_augm)[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
  }
  
  #Calculating each sigma_i separately
  sigmahat_vector <- c()
  for (i in 1:n_ts){
    AR.struc        <- estimate_lrv(data = gdp_mat_augm[, i], q = q, r_bar = r,
                                    p = order_results[i])
                                   #p = 1)
    sigma_hat_i     <- sqrt(AR.struc$lrv)
    sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
  }
  return(sigmahat_vector)
}