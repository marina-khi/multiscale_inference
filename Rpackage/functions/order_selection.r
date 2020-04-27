order_selection <- function(data, q = NULL, r = 5:15){
  # Function that calculates different information criterions for the number of time series
  # based on our variance estimator for a range of tuning parameters. 
  #
  # Arguments:
  # data        number of time series in a matrix. Column names of the matrix should be reasonable
  # q           vector of different tuning parameters to analyse. If not supplied,
  #             q is calculated as [2log(T)], [2log(T)] + 1, ..., [2sqrt(T)] + 1.
  # r           vector of different tuning parameters r.bar to analyse. If not supplied,
  #             r is taken to be 5, 6, ..., 15
  #
  # Outputs: list with the elements  
  # values      vector of individual test statistics (hat(psi)/sigmahat) 
  #             for each point (u,h) in the location-bandwidth grid
  # stat_ms     value of the multiscale statistic
  
  Tlen <- nrow(data)
  N_ts <- ncol(data)

  if (is.null(q)){
    q = seq(floor(2 * log(Tlen)), ceiling(2 * sqrt(Tlen)), by = 1) 
  }
  
  result_list <- list()
  order_results <- c()
  if (N_ts == 1){
    list_names <- c('the only time series')
  } else {
    list_names <- colnames(data)
  }
    
  for (j in 1:N_ts){
    criterion_matrix <- expand.grid(q = q, r = r)
    
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
        AR.struc      <- AR_lrv(data=data[, j], q=criterion_matrix$q[[i]], r.bar=criterion_matrix$r[[i]], p=order)
        sigma_eta_hat <- sqrt(AR.struc$vareta)
        FPE <- c(FPE, (sigma_eta_hat^2 * (Tlen + order)) / (Tlen - order))
        AIC <- c(AIC, Tlen * log(sigma_eta_hat^2) + 2 * order)
        AICC <- c(AICC, Tlen * log(sigma_eta_hat^2) + Tlen * (1 + order / Tlen)/(1 - (order +2)/Tlen))
        SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(Tlen) / Tlen)
        HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(Tlen)) / Tlen)
      }
      criterion_matrix$FPE[[i]]  <- which.min(FPE)
      criterion_matrix$AIC[[i]]  <- which.min(AIC)
      criterion_matrix$AICC[[i]] <- which.min(AICC)
      criterion_matrix$SIC[[i]]  <- which.min(SIC)
      criterion_matrix$HQ[[i]]   <- which.min(HQ)
    }
    maxim <- max(criterion_matrix[, 3:7])
    order_results <- c(order_results, maxim)
    cat("For ", list_names[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
    result_list[[list_names[j]]] <- criterion_matrix
  }
  result_list[['orders']] <- order_results
  return(result_list)
}