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

#Create a matrix (for size and power table for example) and write them in the tex file
creating_matrix_and_texing <- function(vect, vect_t, vect_alpha, filename){
  matrix_ <- matrix(vect, nrow = length(vect_t), ncol = length(vect_alpha), byrow = TRUE)
  rownames(matrix_) <- vect_t
  colnames(matrix_) <- vect_alpha
  
  addtorow     <- list()
  addtorow$pos <- list(0, 0)
  addtorow$command <- c("& \\multicolumn{3}{c}{nominal size $\\alpha$} \\\\\n",
                        "$T$ & 0.01 & 0.05 & 0.1 \\\\\n") 
  print.xtable(xtable(matrix_, digits = c(3), align = "cccc"), type = "latex",
               file = filename, add.to.row = addtorow, include.colnames = FALSE)
}

#Truncated functions for speeding up the process
statistics <- function(data, sigma_vec = 1, n_ts = 2, grid = NULL,
                       ijset = NULL, alpha = 0.05, sim_runs = 1000) {
  
  t_len <- nrow(data)
  
  #If grid is not supplied, we construct it by default
  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }
  
  #If ijset is not supplied, we compare all
  #possible pairs of time series.
  if (is.null(ijset)) {
    ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
    ijset <- ijset[ijset$i < ijset$j, ]
  }
  
  psi   <- compute_statistics(data = data, sigma = 1, sigma_vec = sigma_vec,
                              n_ts = n_ts, grid = grid, deriv_order = 0,
                              epidem = FALSE)
  stat  <- psi$stat
  gset_with_values <- psi$gset_with_values
  
  return(list(stat = stat, stat_pairwise = psi$stat_pairwise,
              ijset = ijset, gset_with_values = gset_with_values))
}

statistics_full <- function(data, sigma_vec = 1, n_ts = 2, grid = NULL,
                            ijset = NULL, alpha = c(0.05), sim_runs = 1000) {
  
  t_len <- nrow(data)
  
  #If grid is not supplied, we construct it by default
  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }
  
  #If ijset is not supplied, we compare all
  #possible pairs of time series.
  if (is.null(ijset)) {
    ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
    ijset <- ijset[ijset$i < ijset$j, ]
  }
  
  # Select (1-alpha) quantile of the multiscale statistic under the null
  quantiles <- compute_quantiles(t_len = t_len, grid = grid, n_ts = n_ts,
                                 ijset = ijset, sigma = 1,
                                 sim_runs = sim_runs,
                                 deriv_order = 0,
                                 correction = TRUE, epidem = FALSE)
  
  probs <- as.vector(quantiles$quant[1, ])
  quant <- as.vector(quantiles$quant[2, ])
  
  quant_vec <- c()
  for (alpha_ind in alpha){
    if (sum(probs == (1 - alpha_ind)) == 0)
      pos <- which.min(abs(probs - (1 - alpha_ind)))
    if (sum(probs == (1 - alpha_ind)) != 0)
      pos <- which.max(probs == (1 - alpha_ind))    
    quant_vec <- c(quant_vec, quant[pos])
  }
  
  psi   <- compute_statistics(data = data, sigma = 1, 
                              sigma_vec = sigma_vec, n_ts = n_ts,
                              grid = grid, deriv_order = 0,
                              epidem = FALSE)
  stat  <- psi$stat
  
  return(list(quant = quant_vec, stat = stat, stat_pairwise = psi$stat_pairwise,
              ijset = ijset, sim_runs = sim_runs, grid = grid))
}
