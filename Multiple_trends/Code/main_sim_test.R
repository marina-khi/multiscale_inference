rm(list=ls())

library(multiscale)
library(car)
library(dplyr)
library(Matrix)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/functions.R")

##############################
#Defining necessary constants#
##############################

n_ts     <- 15 #number of different time series for simulation
n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles

different_T     <- c(100, 250, 500) #Different lengths of time series
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
different_b     <- c(0, 0.75, 1.00, 1.25) #Zero is for calculating the size

#For the error process
a_hat <- 0.25 
sigma <- 0.5

#For the covariate process
beta    <- 1
a_hat_x <- 0.5
sigma_x <- 0.5

#Parameters for the estimation of long-run-variance
q     <- 25 
r     <- 10


################################
#Calculating the size and power#
################################

size_and_power_array <- array(NA, dim = c(length(different_T),
                                          length(different_b),
                                          length(different_alpha)),
                              dimnames = list(t = different_T,
                                              b = different_b,
                                              alpha = different_alpha))
#Constructing the set of pairwise comparisons
ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
ijset <- ijset[ijset$i < ijset$j, ]

for (t_len in different_T){
  k <- match(t_len, different_T)

  #Constructing the grid
  u_grid <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len)
  
  #Calculating Gaussian quantiles
  quantiles <- compute_quantiles(t_len = t_len, grid = grid, n_ts = n_ts,
                                 ijset = ijset, sim_runs = sim_runs,
                                 correction = TRUE)
  probs  <- as.vector(quantiles$quant[1, ])
  quants <- as.vector(quantiles$quant[2, ])
  
  for (b in different_b){
    y_matrix           <- matrix(NA, nrow = t_len, ncol = n_ts)
    y_augm_matrix      <- matrix(NA, nrow = t_len, ncol = n_ts)
    x_matrix           <- matrix(NA, nrow = t_len, ncol = n_ts)
    error_matrix       <- matrix(NA, nrow = t_len, ncol = n_ts)
    m_matrix           <- matrix(0, nrow = t_len, ncol = n_ts)
    #Only the first trend function is non-zero:
    m_matrix[, 1]      <- (1:t_len - 0.5 * t_len) * (b / t_len)
    colnames(y_matrix) <- 1:n_ts
    
    
    simulated_statistic = replicate(n_rep, {
      #First we simulate the data
      for (i in 1:n_ts){
        x_matrix[, i]     <- arima.sim(model = list(ar = a_hat_x),
                                       innov = rnorm(t_len, 0, sigma_x),
                                       n = t_len)
        error_matrix[, i] <- arima.sim(model = list(ar = a_hat),
                                       innov = rnorm(t_len, 0, sigma),
                                       n = t_len)
        y_matrix[, i]     <- m_matrix[, i] + beta * x_matrix[, i] + error_matrix[, i]

      }
      
      sigmahat_vector <- c()
      beta_hat        <- c()
      alpha_hat       <- c()
      
      #Now we estimate the parameters
      for (i in 1:n_ts){
        #First differences
        x_diff    <- x_matrix[, i]- dplyr::lag(x_matrix[, i], n = 1, default = NA)
        y_diff    <- y_matrix[, i]- dplyr::lag(y_matrix[, i], n = 1, default = NA)
        
        #Estimating beta
        x_diff_tmp <- as.matrix(x_diff)[-1, ]
        y_diff_tmp <- as.matrix(y_diff)[-1, ]
        
        beta_hat_tmp  <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
        beta_hat      <- c(beta_hat, as.vector(beta_hat_tmp))
        
        #Estimating alpha_i
        alpha_hat_tmp <- mean(y_matrix[, i] - x_matrix[, i] * as.vector(beta_hat_tmp))
        alpha_hat     <- c(alpha_hat, alpha_hat_tmp)
        
        y_augm_matrix[, i]        <- y_matrix[, i] - x_matrix[, i] * as.vector(beta_hat_tmp) - alpha_hat_tmp
        AR.struc            <- estimate_lrv(data = y_augm_matrix[, i], q = q,
                                            r_bar = r, p = 1)
        sigma_hat_i         <- sqrt(AR.struc$lrv)
        sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)      
      }
      
      psi     <- compute_statistics(data = y_augm_matrix,
                                    sigma_vec = sigmahat_vector,
                                    n_ts = n_ts, grid = grid)
      results <- max(psi$stat_pairwise)
      results
    })
    
    size_and_power_vec <- c()
    for (alpha in different_alpha){
      if (sum(probs == (1 - alpha)) == 0)
        pos <- which.min(abs(probs - (1 - alpha)))
      if (sum(probs == (1 - alpha)) != 0)
        pos <- which.max(probs == (1 - alpha))    
      quant <- quants[pos]
      
      num_of_rej         <- sum(simulated_statistic > quant)/n_rep
      size_and_power_vec <- c(size_and_power_vec, num_of_rej) 
      
      cat("Ratio of rejection is ", num_of_rej, "with b = ", b,
          ", alpha = ", alpha, "and T = ", t_len, "\n")
    }
    
    #Storing the results in a 3D array
    l <- match(b, different_b)
    size_and_power_array[k, l, ] <- size_and_power_vec
  }
}

#######################
#Output of the results#
#######################

for (b in different_b){
  l   <- match(b, different_b)
  tmp <- as.matrix(size_and_power_array[, l, ])
  if (b == 0){
    filename = paste0("output/tables/", n_ts, "_ts_size.tex")
  } else {
    filename = paste0("output/tables/", n_ts, "_ts_power_b_",
                      b * 100, ".tex")
  }
  output_matrix(tmp, filename)
}