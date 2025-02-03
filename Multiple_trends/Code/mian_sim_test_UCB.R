rm(list=ls())

library(MSinference)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
library(Rcpp)
library(tictoc)
library(foreach)
library(parallel)
library(doParallel)
library(mvtnorm)

#Load necessary functions  
source("functions/functions.r")
source("functions/functions_other.r")


##############################
#Defining necessary constants#
##############################

n_ts <- 2 #Number of time series

n_rep    <- 1000 #number of simulations for calculating size and power
sim_runs <- 1000 #number of simulations to calculate the Gaussian quantiles for MS test

different_T <- c(100, 250, 500) #Different lengths of time series
alpha       <- 0.05 #Confidence levels
different_b <- c(0, 0.5, 1, 2) #Zero is for calculating the size

#For the error process
a     <- 0.25
sigma <- 0.25

#For the covariate process
beta    <- c(1, 1, 1)
a_x_vec <- c(0.25, 0.25, 0.25) #VAR(1) coefficients
phi     <- 0.25                 #dependence between the innovations

#For the fixed effects
rho <- 0.25 #covariance between the fixed effects

#Parameters for the estimation of long-run-variance
q <- 25 
r <- 10

seed <- 246802468


######################
#Derivative constants#
######################
sigma_vector <- rep(sigma, n_ts)

phi_matrix       <- matrix(phi, nrow = 3, ncol = 3)
diag(phi_matrix) <- 1
a_matrix         <- diag(a_x_vec) 

big_sigma_matrix       <- matrix(rho, nrow = n_ts, ncol = n_ts)
diag(big_sigma_matrix) <- 1


################################
#Calculating the size and power#
################################

size_and_power_array <- array(NA, dim = c(length(different_T),
                                          length(different_b),
                                          1),
                              dimnames = list(t = different_T,
                                              b = different_b,
                                              alpha = alpha))
size_and_power_UCB_array <- array(NA, dim = c(length(different_T),
                                                length(different_b),
                                                1),
                                  dimnames = list(t = different_T,
                                                  b = different_b,
                                                  alpha = alpha))

for (t_len in different_T){
  set.seed(seed)
  k <- match(t_len, different_T)
  
  k_n <- floor(t_len^(1/3))
  m   <- floor(t_len / k_n)
  
  #Constructing the grid
  u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
  
  grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
  
  ####################################
  #Calculating the Gaussian quantiles#
  ####################################
  
  cat("Calculating the Gaussian quantiles\n")
  simulated_pairwise_gaussian <- matrix(NA, nrow = n_ts * n_ts, ncol = sim_runs)
  simulated_gaussian_UCB      <- c()
  
  for (val in 1:sim_runs){
    z_matrix      <- matrix(NA, nrow = t_len, ncol = n_ts)
    z_augm_matrix <- matrix(NA, nrow = t_len, ncol = n_ts)

    for (i in 1:n_ts){
      z_matrix[, i]      <- rnorm(t_len, 0, 1)
      z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    }

    psi <- compute_statistics(data = z_augm_matrix,
                              sigma_vec = rep(1, n_ts),
                              n_ts = n_ts, grid = grid)
    simulated_pairwise_gaussian[, val] <- as.vector(psi$stat_pairwise)
    
    gaussian_UCB <- mapply(UCB_estimation, grid_points,
                           MoreArgs = list(data_p = z_matrix[, 1],
                                           grid_p = grid_points,
                                           bw = 5/t_len))
    simulated_gaussian_UCB <- c(simulated_gaussian_UCB, max(gaussian_UCB))
  }
  
  simulated_gaussian <- apply(simulated_pairwise_gaussian, 2, max)
  
  probs     <- seq(0.5, 0.995, by = 0.005)
  quantiles <- as.vector(quantile(simulated_gaussian, probs = probs))
  quantiles <- rbind(probs, quantiles)
  
  colnames(quantiles) <- NULL
  rownames(quantiles) <- NULL
  
  quants <- as.vector(quantiles[2, ])

  quantiles_UCB <- as.vector(quantile(simulated_gaussian_UCB, probs = probs))
  quantiles_UCB <- rbind(probs, quantiles_UCB)
  
  colnames(quantiles_UCB) <- NULL
  rownames(quantiles_UCB) <- NULL
  
  quants_UCB <- as.vector(quantiles_UCB[2, ])
  
  if (sum(probs == (1 - alpha)) == 0) {
    pos <- which.min(abs(probs - (1 - alpha)))
  } else {
    pos   <- which.max(probs == (1 - alpha))  
  }

  quant     <- quants[pos]
  quant_UCB <- quants_UCB[pos]
  
    
  #################################
  #Testing for different scenarios#
  #################################
  
  for (b in different_b){
    simulated_pairwise_statistics <- matrix(NA, nrow = n_ts * n_ts, ncol = n_rep)
    result_UCB                    <- c()
    
    m_matrix <- matrix(0, nrow = t_len, ncol = n_ts)
    if (b == 0) {
      cat("SIZE SIMULATIONS\n")
    } else {
      cat("POWER SIMULATIONS WITH b = ", b, "\n")
      #Only the first trend function is non-zero:
      m_matrix[, 1] <- bump((1:t_len)/t_len) * b
    }

    for (val in 1:n_rep){
      y_matrix      <- matrix(NA, nrow = t_len, ncol = n_ts)
      y_augm_matrix <- matrix(NA, nrow = t_len, ncol = n_ts)
      error_matrix  <- matrix(NA, nrow = t_len, ncol = n_ts)
      
      alpha_vec     <- rmvnorm(1, mean = rep(0, n_ts), sigma = big_sigma_matrix)
      
      sigmahat_vec     <- rep(NA, n_ts)
      sigmahat_vec_UCB <- rep(NA, n_ts) 
    
      #UNIFORM CONFIDENCE BOUNDS
      estimated_trend_UCB <- matrix(NA, nrow = t_len, ncol = n_ts)
      upper_UCB           <- matrix(NA, nrow = t_len, ncol = n_ts)
      lower_UCB           <- matrix(NA, nrow = t_len, ncol = n_ts)
      
      for (i in 1:n_ts){
        error_matrix[, i] <- arima.sim(model = list(ar = a),
                                       innov = rnorm(t_len, 0, sigma),
                                       n = t_len)
        nu       <- rmvnorm(t_len + 10, mean = c(0, 0, 0), sigma = phi_matrix)
        x_matrix <- matrix(0, 3, t_len + 10)
        
        for (t in 2:(t_len + 10)){
          x_matrix[, t] <- a_matrix %*% x_matrix[, t - 1] + nu[t, ]
        }
        x_matrix <- t(x_matrix[, -(1:10)])
        
        y_matrix[, i] <- alpha_vec[i] + m_matrix[, i] + beta %*% t(x_matrix) + error_matrix[, i]

        #First differences
        y_diff_tmp <- y_matrix[, i] - dplyr::lag(y_matrix[, i], n = 1, default = NA)
        y_diff     <- as.matrix(y_diff_tmp)[-1, ]
        x_diff_1   <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
        x_diff_2   <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
        x_diff_3   <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)

        #Estimating beta
        x_diff    <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]
        beta_hat  <- solve(t(x_diff) %*% x_diff) %*% t(x_diff) %*% y_diff
        alpha_hat <- mean(y_matrix[, i] - x_matrix %*% as.vector(beta_hat))
            
        y_augm_matrix[, i] <- y_matrix[, i] - x_matrix %*% as.vector(beta_hat) - alpha_hat
            
        #Estimating the variance
        AR.struc        <- estimate_lrv(data = y_augm_matrix[, i], q = q,
                                           r_bar = r, p = 1)
        sigma_hat_i     <- sqrt(AR.struc$lrv)
        sigmahat_vec[i] <- sigma_hat_i
        sigma_hat_UCB_i <- sqrt(sigma_estimation_UCB(y_ = y_matrix[, i],
                                                     x_matrix_ = x_matrix,
                                                     beta_est_ = beta_hat,
                                                     m_ = m, k_n_ = k_n)) 
        sigmahat_vec_UCB[i] <- sigma_hat_UCB_i
                
        estimated_trend_UCB[, i] <- mapply(UCB_estimation, grid_points,
                                           MoreArgs = list(data_p = y_augm_matrix[, i],
                                                           grid_p = grid_points,
                                                           bw = 5/t_len))
        upper_UCB[, i] <- estimated_trend_UCB[, i] + sigma_hat_UCB_i * quant_UCB
        lower_UCB[, i] <- estimated_trend_UCB[, i] - sigma_hat_UCB_i * quant_UCB
        # plot(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
        #      y = estimated_trend_UCB[, i], ylim = c(-1.6, 1.6),
        #      xlab = "", ylab = "", main = NULL,
        #      type = 'l', cex = 0.8)
        # lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
        #       y = upper_UCB[, i], type = "l",
        #       col = "red")
        # lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
        #       y = lower_UCB[, i], type = "l",
        #       col = "red")        
      }
      #MULTISCALE TEST
      psi <- compute_statistics(data = y_augm_matrix,
                                sigma_vec = sigmahat_vec,
                                n_ts = n_ts, grid = grid)    
      simulated_pairwise_statistics[, val] <- as.vector(psi$stat_pairwise)
      result_UCB <- c(result_UCB, (sum((lower_UCB[, 1] < upper_UCB[, 2]) & (lower_UCB[, 2] < upper_UCB[, 1])) == 0))
    }
    
    simulated_statistic <- apply(simulated_pairwise_statistics[1:(n_ts * n_ts), ], 2, max)
    
    size_and_power_vec     <- c()
    size_and_power_UCB_vec <- c()
    
    num_of_rej         <- sum(simulated_statistic > quant)/n_rep
    size_and_power_vec <- c(size_and_power_vec, num_of_rej) 
      
    cat("Ratio of rejection is ", num_of_rej, "with b = ", b,
        ", alpha = ", alpha, "and T = ", t_len, "\n")
      
    num_of_rej_UCB         <- sum(result_UCB)/n_rep
    size_and_power_UCB_vec <- c(size_and_power_UCB_vec, num_of_rej_UCB) 
      
    cat("Ratio of rejection for UCB is ", num_of_rej_UCB, "with b = ", b,
        ", alpha = ", alpha, "and T = ", t_len, "\n")
    
    #Storing the results in a 3D array
    l <- match(b, different_b)
    size_and_power_array[k, l, ]     <- size_and_power_vec
    size_and_power_UCB_array[k, l, ] <- size_and_power_UCB_vec
  }
} 


#######################
#Output of the results#
#######################

for (b in different_b){
  l   <- match(b, different_b)
  tmp <- matrix(NA, nrow = length(different_T), ncol = 2 * length(different_alpha))
  for (i in 1:length(different_alpha)){
    tmp[, 2 * i - 1] <- as.vector(size_and_power_array[, l, i])
    tmp[, 2 * i]     <- as.vector(size_and_power_sizer_array[, l, i])
  }
  
  row.names(tmp) <- paste0("$T = ", row.names(as.matrix(size_and_power_array[, l, ])), "$")
  
  tmp2 <- as.matrix(size_and_power_array[, l, ])
  tmp3 <- as.matrix(size_and_power_sizer_array[, l, ])
  
  if (b == 0){
    filename = paste0("output/revision/", n_ts, "_ts_size_UCB_comparison.tex")
  } else {
    filename = paste0("output/revision/", n_ts, "_ts_power_b_",
                      b * 100, "_UCB_comparison.tex")
  }
  output_matrix(tmp, filename, numcols_ = 7)
  line <- paste0("%This simulation was done for the following values of the parameters: n_ts = ", n_ts,
                 ", with ", n_rep, " simulations for calculating size and power and ", sim_runs,
                 " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
                 a, " and sigma = ", sigma,
                 ". There are no fixed effects. The grid is normal.")
  write(line, file = filename, append = TRUE)
}