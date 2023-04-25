#####################
#Estimating variance#
#####################

sigma <- function(matrix_, n_ts, q_ = 15, r_ = 10, procedure = "all_one"){
  #Order selection
  q_vec <- 15:25
  r_vec <- 10:15
  order_results <- c()
  order_results_BIC <- c()
  
  for (j in 1:n_ts){
    criterion_matrix <- expand.grid(q = q_vec, r = r_vec)
    
    criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$BIC  <- numeric(length = nrow(criterion_matrix))
    criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
    
    for (i in 1:nrow(criterion_matrix)){
      FPE <- c()
      AIC <- c()
      AICC <- c()
      BIC <- c()
      HQ <- c()
      
      different_orders <- (1:9)
      
      for (order in different_orders){
        AR.struc      <- estimate_lrv(data = matrix_[, j], q = criterion_matrix$q[[i]],
                                      r_bar = criterion_matrix$r[[i]], p = order)
        sigma_eta_hat <- sqrt(AR.struc$vareta)
        FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
        AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
        AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order +2)/t_len))
        BIC <- c(BIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
        HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
      }
      criterion_matrix$FPE[[i]]  <- which.min(FPE)
      criterion_matrix$AIC[[i]]  <- which.min(AIC)
      criterion_matrix$AICC[[i]] <- which.min(AICC)
      criterion_matrix$BIC[[i]]  <- which.min(BIC)
      criterion_matrix$HQ[[i]]   <- which.min(HQ)
    }
    maxim <- max(criterion_matrix[, 3:7])
    order_results <- c(order_results, maxim)
    order_results_BIC <- c(order_results_BIC, max(criterion_matrix$BIC))
    cat("For the country ", colnames(matrix_)[j], " the results are as follows: ",
        max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ",
        max(criterion_matrix$AICC), " ", max(criterion_matrix$BIC), " ",
        max(criterion_matrix$HQ), " \n")
  }
  
  #Calculating each sigma_hat_i separately
  sigmahat_vector <- c()
  ahat_vector     <- c()
  
  if (procedure == "all_one") {
    for (i in 1:n_ts){
      AR.struc        <- estimate_lrv(data = matrix_[, i], q = q_, r_bar = r_,
                                      p = 1)
      sigma_hat_i     <- sqrt(AR.struc$lrv)
      sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
      ahat_vector     <- c(ahat_vector, AR.struc$ahat)
    }
  } else {
    for (i in 1:n_ts){
      AR.struc        <- estimate_lrv(data = matrix_[, i], q = q_, r_bar = r_,
                                      p = order_results_BIC[i])
      sigma_hat_i     <- sqrt(AR.struc$lrv)
      sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
      ahat_vector     <- c(ahat_vector, AR.struc$ahat)
    }    
  }
  return(sigmahat_vector)
}