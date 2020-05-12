#' Calculates estimators of the standard deviations for a number of iid samples. 
#' It uses the standard first differences method.
#' @param data      Matrix of the samples. Each sample is located in one of the columns.
#' @param method    Method of estimation, either 'first', then first differences are used, or 'second',
#'                  then the two neighbours method is used, as described in Brown and Levine (2007)
#' @return lrv      Estimator of the long run variance of the error terms.
#' @return ahat     Vector of length p of estimated AR coefficients.
#' @return vareta   Estimator of the variance of the innovation term

estimate_iid_sd <- function(data, method = 'first'){
  N    = ncol(data)
  Tlen = nrow(data)
  sigmahat_vec <- c()
  if (method == 'first') {
    for (i in 1:N){
      variance_i   <- sum((data[2:Tlen, i] - data[1:(Tlen - 1), i])^2)/(2 * Tlen - 2)
      sigma_hat_i  <- sqrt(variance_i)
      sigmahat_vec <- c(sigmahat_vec, sigma_hat_i)
    }
  } else if (method == 'second'){
    for (i in 1:N){
      variance_i   <- 2 * sum((data[3:Tlen, i]/2 - data[2:(Tlen - 1), i] + data[1:(Tlen - 2), i]/2)^2)/(3 * Tlen - 6)
      sigma_hat_i  <- sqrt(variance_i)
      sigmahat_vec <- c(sigmahat_vec, sigma_hat_i)
    }
  } else {
    sigmahat_vec <- NULL
  }
  return(sigmahat_vec)
}