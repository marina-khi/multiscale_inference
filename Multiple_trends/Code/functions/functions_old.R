# Unnecessary pieces of code
#
#
#
#
# ################################################
# #Estimating the parameters from the application#
# ################################################
# 
# source("functions/data_loading_hp.R") #Now matrix hp_data contains all data
# 
# dates     <- unique(hp_data$year)
# n_ts      <- length(unique(hp_data$iso))
# t_len     <- nrow(hp_data) / n_ts
# countries <- unique(hp_data$iso)
# 
# q         <- 15 #Parameters for the estimation of sigma
# r         <- 10
# 
# #Estimating alpha and beta (Needed only for getting rid of the covariates' effect)
# source("functions/parameters_estimation.R")
# estimated   <- parameters_hp(hp_data, n_ts, t_len, countries)
# hp_log      <- estimated$hp_log #Original time series
# hp_log_augm <- estimated$hp_log_augm   #Augmented time series
# 
# gdp_log           <- matrix(NA, ncol = n_ts, nrow = t_len)
# colnames(gdp_log) <- countries
# 
# pop_log           <- matrix(NA, ncol = n_ts, nrow = t_len)
# colnames(pop_log) <- countries
# 
# ltrate            <- matrix(NA, ncol = n_ts, nrow = t_len)
# colnames(ltrate)  <- countries
# 
# infl              <- matrix(NA, ncol = n_ts, nrow = t_len)
# colnames(infl)    <- countries
# 
# i <- 1 
# for (country in countries){
#   tmp <- hp_data[hp_data$iso == country, ]
#   tmp <- tmp[order(tmp$year), ]
#   gdp_log[, i] <- tmp$log_gdp
#   pop_log[, i] <- tmp$log_pop
#   ltrate[, i]  <- tmp$ltrate
#   infl[, i]    <- tmp$infl
#   i = i + 1
# }
# 
# a_coef_matrix <- matrix(NA, ncol = 5, nrow = n_ts)
# colnames(a_coef_matrix) <- c("gdp_log", "pop_log", "ltrate", "infl", "error_process")
# rownames(a_coef_matrix) <- countries
# 
# innovation_var_matrix <- matrix(NA, ncol = 5, nrow = n_ts)
# colnames(innovation_var_matrix) <- c("gdp_log", "pop_log", "ltrate", "infl", "error_process")
# rownames(innovation_var_matrix) <- countries
# 
# sigmahat_vector <- c()
# 
# #Estimating coefficients for all covariates
# for (i in 1:n_ts){
#   AR.struc_gdp_log <- estimate_lrv(data = gdp_log[, i], q = q, r_bar = r, p = 1)
#   AR.struc_pop_log <- estimate_lrv(data = pop_log[, i], q = q, r_bar = r, p = 1)
#   AR.struc_ltrate  <- estimate_lrv(data = ltrate[, i], q = q, r_bar = r, p = 1)
#   AR.struc_infl    <- estimate_lrv(data = infl[, i], q = q, r_bar = r, p = 1)
#   AR.struc_errors  <- estimate_lrv(data = hp_log_augm[, i], q = q, r_bar = r, p = 1)
#   
#   a_coef_matrix[i, 1] <- AR.struc_gdp_log$ahat
#   a_coef_matrix[i, 2] <- AR.struc_pop_log$ahat
#   a_coef_matrix[i, 3] <- AR.struc_ltrate$ahat
#   a_coef_matrix[i, 4] <- AR.struc_infl$ahat
#   a_coef_matrix[i, 5] <- AR.struc_errors$ahat
#   
#   innovation_var_matrix[i, 1] <- AR.struc_gdp_log$vareta
#   innovation_var_matrix[i, 2] <- AR.struc_pop_log$vareta
#   innovation_var_matrix[i, 3] <- AR.struc_ltrate$vareta
#   innovation_var_matrix[i, 4] <- AR.struc_infl$vareta
#   innovation_var_matrix[i, 5] <- AR.struc_errors$vareta
#   
#   sigma_hat_i     <- sqrt(AR.struc_errors$lrv)
#   sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
# }
# 
# a_coef_vector         <- colSums(a_coef_matrix)/length(countries)
# innovation_var_vector <- colSums(innovation_var_matrix)/length(countries)
# 
# a_coef_matrix         <- rbind(a_coef_matrix, mean = a_coef_vector)
# innovation_var_matrix <- rbind(innovation_var_matrix, mean = innovation_var_vector)
# 
# 
# addtorow     <- list()
# addtorow$pos <- list(0)
# addtorow$command <- c("Country & GDP & Population & Long-term rate & Inflation & Error process\\\\\n") 
# print.xtable(xtable(a_coef_matrix, digits = c(3), align = "cccccc"),
#              type = "latex", file = "output/revision/a_coefficients.tex",
#              add.to.row = addtorow, include.colnames = FALSE)
# 
# addtorow     <- list()
# addtorow$pos <- list(0)
# addtorow$command <- c("Country & GDP & Population & Long-term rate & Inflation & Error process\\\\\n") 
# print.xtable(xtable(sqrt(innovation_var_matrix), digits = c(3), align = "cccccc"),
#              type = "latex", file = "output/revision/innovation_sds.tex",
#              add.to.row = addtorow, include.colnames = FALSE)
#
# #Function that simulates 4 covariates as AR(1) with the given
# #coefficients (a_x_vec_ and sigma_x_vec_),
# #the error terms as AR(1) also with the given coefficient a_ and sigma_,
# #the time series as y = beta_ %*% covariates + m_matrix_ + errors,
# #estimates the parameters, and then computes the test statistics
# repl_revision <- function(rep, t_len_, n_ts_, grid_, gaussian_sim = FALSE,
#                           a_ = 0, sigma_ = 1, beta_ = NULL,
#                           a_x_vec_ = c(0, 0, 0, 0),
#                           sigma_x_vec_ = c(1, 1, 1, 1),
#                           rho_ = 0, 
#                           m_matrix_ = NULL, q_ = 25, r_ = 10){
# 
#   library(MSinference)
#   library(dplyr)
#   
#   if (gaussian_sim){
#     z_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#     z_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#     sigma_vector  <- rep(1, n_ts_)
#     
#     for (i in 1:n_ts_){
#       z_matrix[, i]      <- rnorm(t_len_, 0, sigma_)
#       z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
#     }
#     
#     psi <- compute_statistics(data = z_augm_matrix,
#                               sigma_vec = sigma_vector,
#                               n_ts = n_ts_, grid = grid_)
#   } else {
#     y_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#     y_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#     error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#     
#     beta_hat_matrix <- matrix(NA, nrow = 4, ncol = n_ts_)
#     
#     library(mvtnorm)
#     big_sigma_matrix       <- matrix(rho_, nrow = n_ts_, ncol = n_ts_)    
#     diag(big_sigma_matrix) <- 1
#     alpha_vec <- rmvnorm(1, mean = rep(0, n_ts_), sigma = big_sigma_matrix)
#     
#     sigmahat_vector <- c()
#     
#     if (!is.null(beta_)){
#       x_matrix_1 <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#       x_matrix_2 <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#       x_matrix_3 <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#       x_matrix_4 <- matrix(NA, nrow = t_len_, ncol = n_ts_)
#       for (i in 1:n_ts_){
#         error_matrix[, i] <- arima.sim(model = list(ar = a_),
#                                        innov = rnorm(t_len_, 0, sigma_),
#                                        n = t_len_)
#         x_matrix_1[, i]   <- arima.sim(model = list(ar = a_x_vec_[1]),
#                                       innov = rnorm(t_len_, 0, sigma_x_vec_[1]),
#                                       n = t_len_)
#         x_matrix_2[, i]   <- arima.sim(model = list(ar = a_x_vec_[2]),
#                                        innov = rnorm(t_len_, 0, sigma_x_vec_[2]),
#                                        n = t_len_)
#         x_matrix_3[, i]   <- arima.sim(model = list(ar = a_x_vec_[3]),
#                                        innov = rnorm(t_len_, 0, sigma_x_vec_[3]),
#                                        n = t_len_)
#         x_matrix_4[, i]   <- arima.sim(model = list(ar = a_x_vec_[4]),
#                                        innov = rnorm(t_len_, 0, sigma_x_vec_[4]),
#                                        n = t_len_)
#         x_matrix          <- cbind(x_matrix_1[, i], x_matrix_2[, i], x_matrix_3[, i], x_matrix_4[, i])
#         
#         y_matrix[, i]     <- alpha_vec[i] + m_matrix_[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
#         
#         #First differences
#         x_diff_1  <- x_matrix_1[, i] - dplyr::lag(x_matrix_1[, i], n = 1, default = NA)
#         x_diff_2  <- x_matrix_2[, i] - dplyr::lag(x_matrix_2[, i], n = 1, default = NA)
#         x_diff_3  <- x_matrix_3[, i] - dplyr::lag(x_matrix_3[, i], n = 1, default = NA)
#         x_diff_4  <- x_matrix_4[, i] - dplyr::lag(x_matrix_4[, i], n = 1, default = NA)
#         y_diff    <- y_matrix[, i] - dplyr::lag(y_matrix[, i], n = 1, default = NA)
#         
#         #Estimating beta
#         x_diff_tmp <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3, x_diff_4))[-1, ]
#         y_diff_tmp <- as.matrix(y_diff)[-1, ]
#         
#         beta_hat_tmp         <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
#         beta_hat_matrix[, i] <- beta_hat_tmp
#         alpha_hat_tmp        <- mean(y_matrix[, i] - x_matrix %*% as.vector(beta_hat_tmp))
#         
#         y_augm_matrix[, i] <- y_matrix[, i] - x_matrix %*% as.vector(beta_hat_tmp) - alpha_hat_tmp
#         AR.struc           <- estimate_lrv(data = y_augm_matrix[, i], q = q_,
#                                            r_bar = r_, p = 1)
#         sigma_hat_i        <- sqrt(AR.struc$lrv)
#         sigmahat_vector    <- c(sigmahat_vector, sigma_hat_i) 
#       }
#     } else {
#       for (i in 1:n_ts_){
#         error_matrix[, i]  <- arima.sim(model = list(ar = a_),
#                                         innov = rnorm(t_len_, 0, sigma_),
#                                         n = t_len_)
#         y_matrix[, i]      <- alpha_vec[i] + m_matrix_[, i] + error_matrix[, i]
#         
#         #Estimating alpha_i
#         alpha_hat_tmp      <- mean(y_matrix[, i])
#         
#         y_augm_matrix[, i] <- y_matrix[, i] - alpha_hat_tmp
#         AR.struc           <- estimate_lrv(data = y_augm_matrix[, i], q = q_,
#                                            r_bar = r_, p = 1)
#         sigma_hat_i        <- sqrt(AR.struc$lrv)
#         sigmahat_vector    <- c(sigmahat_vector, sigma_hat_i)   
#       }     
#     }
#     psi <- compute_statistics(data = y_augm_matrix,
#                               sigma_vec = sigmahat_vector,
#                               n_ts = n_ts_, grid = grid_)    
#   }
#   results <- c(as.vector(psi$stat_pairwise), as.vector(beta_hat_matrix))
#   return(results)
# }
# #Function that plots the graphs with histograms
# plot_histogram <- function(pdfname, data_matrix, n_ts, names, star_value){
#   smallest_value <- floor(min(data_matrix)*5)/5
#   biggest_value  <- ceiling(max(data_matrix)*5)/5
#   
#   if (biggest_value - smallest_value <= 0.6){
#     breaks_grid <- seq(smallest_value, biggest_value, by = 0.05)    
#   } else {
#     breaks_grid <- seq(smallest_value, biggest_value, by = 0.2)
#   }
#   breaks_grid[length(breaks_grid)] <- biggest_value
#   
#   pdf(pdfname, width = 6, height = 9, paper="special")
#   par(mfrow = c(n_ts/2, 2))
#   par(mar = c(3, 3, 0.5, 0)) #Margins for each plot
#   par(oma = c(0.5, 0.5, 0.5, 0.2)) #Outer margins
#   for (i in 1:n_ts){
#     hist_ <- hist(data_matrix[, i], breaks = breaks_grid, plot = FALSE)
#     highestCount <- max(hist_$counts)
#     hist(data_matrix[, i], main = NULL, breaks = breaks_grid, freq=TRUE,
#          xlim = c(smallest_value, biggest_value), ylim = c(0,highestCount),
#          xlab = "", mgp = c(2, 0.5, 0), cex.lab = 1.1)
#     mtext(side = 1, text = names[i], line = 1.5, cex = 0.6)
#     segments(x0 = star_value, y0 = 0, x1 = star_value, y1 = highestCount,
#              col = "red", lwd = 1.5)
#   }
#   dev.off()
# }

#Function that simulates the covariates as AR(1), the error terms as AR(1)
#the time series as y = beta_ %*% covariates + m_matrix_ + errors,
#estimates the parameters, and then computes the test statistics
repl2 <- function(rep, t_len_, n_ts_, grid_, gaussian_sim = FALSE,
                  a_ = 0, sigma_ = 1, beta_ = 0, a_x_ = 0, sigma_x_ = 1,
                  m_matrix_ = NULL, q_ = 25, r_ = 10){
  library(MSinference)
  library(dplyr)
  
  if (gaussian_sim){
    z_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    z_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    sigma_vector  <- rep(1, n_ts_)
    
    for (i in 1:n_ts_){
      z_matrix[, i]      <- rnorm(t_len_, 0, sigma_)
      z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    }
    
    psi <- compute_statistics(data = z_augm_matrix,
                              sigma_vec = sigma_vector,
                              n_ts = n_ts_, grid = grid_)
  } else {
    y_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    y_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    
    sigmahat_vector <- c()
    
    if (beta_ != 0){
      x_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      for (i in 1:n_ts_){
        error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                       innov = rnorm(t_len_, 0, sigma_),
                                       n = t_len_)
        x_matrix[, i]     <- arima.sim(model = list(ar = a_x_),
                                       innov = rnorm(t_len_, 0, sigma_x_),
                                       n = t_len_)
        y_matrix[, i]     <- m_matrix_[, i] + beta_ * x_matrix[, i] + error_matrix[, i]
        
        #First differences
        x_diff    <- x_matrix[, i]- dplyr::lag(x_matrix[, i], n = 1, default = NA)
        y_diff    <- y_matrix[, i]- dplyr::lag(y_matrix[, i], n = 1, default = NA)
        
        #Estimating beta
        x_diff_tmp <- as.matrix(x_diff)[-1, ]
        y_diff_tmp <- as.matrix(y_diff)[-1, ]
        
        beta_hat_tmp  <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
        alpha_hat_tmp <- mean(y_matrix[, i] - x_matrix[, i] * as.vector(beta_hat_tmp))
        
        y_augm_matrix[, i] <- y_matrix[, i] - x_matrix[, i] * as.vector(beta_hat_tmp) - alpha_hat_tmp
        AR.struc           <- estimate_lrv(data = y_augm_matrix[, i], q = q_,
                                           r_bar = r_, p = 1)
        sigma_hat_i        <- sqrt(AR.struc$lrv)
        sigmahat_vector    <- c(sigmahat_vector, sigma_hat_i) 
      }
    } else {
      for (i in 1:n_ts_){
        error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                       innov = rnorm(t_len_, 0, sigma_),
                                       n = t_len_)
        y_matrix[, i]     <- m_matrix_[, i] + error_matrix[, i]
        
        #Estimating alpha_i
        alpha_hat_tmp     <- mean(y_matrix[, i])
        
        y_augm_matrix[, i]  <- y_matrix[, i] - alpha_hat_tmp
        AR.struc            <- estimate_lrv(data = y_augm_matrix[, i], q = q_,
                                            r_bar = r_, p = 1)
        sigma_hat_i         <- sqrt(AR.struc$lrv)
        sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)   
      }     
    }
    psi <- compute_statistics(data = y_augm_matrix,
                              sigma_vec = sigmahat_vector,
                              n_ts = n_ts_, grid = grid_)    
  }
  results <- as.vector(psi$stat_pairwise)
  return(results)
}

CalculateSize <- function(timeseries1, timeseries2, T, different_alpha, Nsim = 1000,
                          SimRuns = 1000, q_ = 25, r_ = 15,
                          remove.small.ess =FALSE){
  # Computes the statistical size of the multiscale test (T_ms), the uncorrected version of the multiscale test (T_uc),
  #   the rowwise version of the multiscale test (T_rw) compared to the size of SiZer test (T_SiZer). The 
  #   The size is calculated only for the case where the errors follow AR(1) process.
  #   IMPORTANT: For the computational purposes, the size of SiZer is calculated only for sample sies <= 1000.
  #
  # Args:
  #   T:  Sample size of the time series simulated.
  #   a1: AR(1) coefficient for simulating the error distribution.
  #   sigma_eta: standard deviation of the innovation term in the AR(1) process.
  #   different_alpha: Different confidence levels. Default is 5%.
  #   Nsim: number of simulations for power calculations. Default is 1000.
  #   SimRuns: number of simulations for calculating gaussian quantile for T_ms. Default is 5000.
  #   sigma.type: If 'estimated', then the long-run variance used in the test statistic is first estimated by
  #     AR_lrv() function. Otherwise, if 'true' the true theoretical value of the long-run variance is used.
  #   q_: tuning parameter for estimating the long-run variance from Section 4. Default is 25. 
  #   remove.small.ess: If TRUE, then we restrict attention only to those points (u, h) for which 
  #     the effective sample size for correlated data ESS^*(u, h) is not smaller than 5. Default is FALSE.
  #
  # Returns:
  #   size.ms:    A vector of length equal to the number of different significant levels alpha, each entry corresponding to the size 
  #     of the multiscale testing procedure (T_ms) for this alpha.
  #   size.SiZer: A vector of length equal to the number of different significant levels alpha, each entry corresponding to the size 
  #     of the SiZer testing procedure (T_SiZer) for this alpha.
  #   size.rw.SiZer: A list of length equal to the number of different significant levels alpha, each entry being a vector of length
  #     equal to the number of bandwidths analysed. The vector contains rowwise size of the SiZer testing procedure (T_SiZer).
  #   h.grid.new:  A vector of the bandwidths analysed.
  
  #Load necessary functions  
  # source("functions/ConstructGrid.r")
  # source("functions/multiscale_statistics.r")
  # source("functions/multiscale_quantiles.r")
  # source("functions/multiscale_testing.r")
  # source("functions/long_run_variance.r")
  # source("functions/sim.r")
  
  sourceCpp("functions/kernel_weights.cpp")
  sourceCpp("functions/SiZer_functions.cpp")
  
  #Construct grid
  grid      <- grid_construction(T)
  gset      <- grid$gset
  u.grid    <- sort(unique(gset[,1]))
  h.grid    <- sort(unique(gset[,2]))
  
  AR.struc1   <- estimate_lrv(data = timeseries1, q = q_, r_bar = r_, p = 1)
  sigma_hat_1 <- sqrt(AR.struc1$lrv)
  a_1         <- AR.struc1$ahat
  autocov1    <- (sigma_hat_1^2/(1-a_1^2)) * (a_1^seq(0,T-1,by=1))  
  
  AR.struc2   <- estimate_lrv(data = timeseries2, q = q_, r_bar = r_, p = 1)
  sigma_hat_2 <- sqrt(AR.struc2$lrv)
  a_2         <- AR.struc2$ahat
  autocov2    <- (sigma_hat_2^2/(1-a_2^2)) * (a_2^seq(0,T-1,by=1))  
  
  if (remove.small.ess){
    ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
    deletions <- ess$del
    grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
  }
  
  gset        <- grid$gset
  N           <- as.integer(dim(gset)[1])
  h.grid.new  <- sort(unique(grid$gset[,2]))
  
  # Compute kernel weights and critical value for multiscale test
  T                      <- as.integer(T) 
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp) 
  storage.mode(gset_cpp) <- "double"
  
  #wghts <- matrix(kernel_weights(T, gset_cpp, N), ncol = T, byrow = TRUE)
  
  sizer.wghts  <- SiZer_weights(T=T, grid=grid)
  sizer.std    <- SiZer_std(weights = sizer.wghts, autocov1 = autocov1,
                            autocov2 = autocov2, t_len = T)
  
  sizer.quants <- vector("list", length(different_alpha))
  for (k in 1:length(different_alpha)){
    sizer.quants[[k]] <- SiZer_quantiles(alpha=different_alpha[k], T=T,
                                         grid=grid, autocov1=autocov1,
                                         autocov2=autocov2)
  }
  
  size_matrix_temp_ms    <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
  size_matrix_temp_SiZer <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
  
  cat("","\n")
  cat("Carrying out the size simulations for the following specification: a_1 = ", a1, ", T = ", T,"\n")
  progbar <- txtProgressBar(min = 1, max = Nsim, style = 3, char = ".")
  
  for (i in 1:Nsim){
    #Simulating the time series
    data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
    data           <- data.simulated$data
    trend          <- data.simulated$trend
    
    #Estimating the coefficients for the ts and the long-run variance
    AR.struc  <- AR_lrv(data = data, q=q_, r.bar = 10, p = 1)
    a.hat     <- AR.struc$ahat
    vareta    <- AR.struc$vareta
    sigma_hat <- sqrt(AR.struc$lrv)
    stats     <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid)
    
    #Values of the multiscale statistic  
    vals     <- stats$values
    
    #Values of SiZer
    values     <- sizer.wghts %*% data
    sizer.vals <- as.vector(values)
    
    #Based on the same values of the test statistic, perform the test at different significance levels
    for (j in 1:length(different_alpha)){
      alpha         <- different_alpha[j]
      test.res      <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
      
      if (sum(abs(test.res$test_ms) == 1, na.rm = TRUE) > 0) {size_matrix_temp_ms[i, j] <- 1} else {size_matrix_temp_ms[i, j] <- 0}
      if (sum(abs(test.res$test_rows) == 1, na.rm = TRUE) > 0) {size_matrix_temp_rows[i, j] <- 1} else {size_matrix_temp_rows[i, j] <- 0}
      
      SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants[[j]], grid=grid)
      if (sum(abs(SiZer_results$test) == 1, na.rm = TRUE) >0) {size_matrix_temp_SiZer[i, j] <- 1} else {size_matrix_temp_SiZer[i, j] <- 0}
      
      for (k in 1:length(h.grid.new)){
        bw      <- h.grid.new[k]
        h_index <- match(bw, h.grid)
        
        if (sum(abs(test.res$test_ms[h_index, ]) == 1, na.rm = TRUE) > 0) {size_rowwise_temp_ms[[j]][i, k] <- 1} else {size_rowwise_temp_ms[[j]][i, k] <- 0}
        if (sum(abs(SiZer_results$test[h_index, ])== 1, na.rm = TRUE) > 0) {size_rowwise_temp_SiZer[[j]][i, k] <- 1} else {size_rowwise_temp_SiZer[[j]][i, k] <- 0}
      }
    }
    setTxtProgressBar(progbar, i)
  }
  close(progbar)
  
  size_ms    <- colSums(size_matrix_temp_ms)/Nsim
  size_SiZer <- colSums(size_matrix_temp_SiZer)/Nsim
  
  rm(size_matrix_temp_ms, size_matrix_temp_SiZer)
  
  return(list(size.ms = size_ms, size.SiZer = size_SiZer,
              h.grid = h.grid.new))
}

##################
#SIZER COMPARISON#
##################
SiZermap <- function(u.grid, h.grid, test.results, plot.title = NA){
  # computes SiZer map from the test results 
  
  col.vec <- c("red", "purple", "blue", "gray") 
  #col.vec <- c("#F7F7F7", "#969696", "#525252", "#636363") 
  temp    <- sort(unique(as.vector(test.results))) + 2
  temp    <- seq(min(temp),max(temp),by=1)
  col.vec <- col.vec[temp]
  
  image(x = u.grid, y = log(h.grid,10), z = t(test.results), col = col.vec,
        xlab = '', ylab = expression(log[10](h)), main = plot.title, xaxt = 'n',
        mgp = c(1,0.5,0))
}

# #For smoothing
# grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
# y_aggregated_vec <- rowMeans(y_matrix)
# smoothed         <- mapply(local_linear_smoothing, grid_points,
#                            MoreArgs = list(y_aggregated_vec, grid_points, bw = 7/t_len))
# 
# y_sizer_matrix[, 1] <- rep(y_aggregated_vec - smoothed, n_ts)

# #Local linear estimates for each of the ts
# smoothed_i  <- mapply(local_linear_smoothing, grid_points,
#                       MoreArgs = list(y_matrix[, i], grid_points, bw = 7/t_len))
# y_sizer_matrix[((i - 1) * t_len + 1):(i * t_len), 2] <- y_augm_matrix[, i] - smoothed_i