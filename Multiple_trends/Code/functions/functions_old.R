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
