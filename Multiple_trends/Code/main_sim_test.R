rm(list=ls())

library(multiscale)
library(car)
library(Matrix)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/functions.R")

##############################
#Defining necessary constants#
##############################

n_ts     <- 15 #number of different time series for simulation
n_rep    <- 1000 #number of simulations for calculating size and power
sim_runs <- 1000 #number of simulations to calculate the Gaussian quantiles

different_T     <- c(100, 250, 500) #Different lengths of time series
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
different_beta  <- c(0, 0.75, 1.00, 1.25) #Zero is for calculating the size

a_hat <- 0.25 
sigma <- 0.5
q     <- 25 #Parameters for the estimation of long-run-variance
r     <- 10


################################
#Calculating the size and power#
################################

size_and_power_array <- array(NA, dim = c(length(different_T),
                                          length(different_beta),
                                          length(different_alpha)),
                              dimnames = list(t = different_T,
                                              beta = different_beta,
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
  
  for (beta in different_beta){
    simulated_data           <- matrix(NA, nrow = t_len, ncol = n_ts)
    colnames(simulated_data) <- 1:n_ts
    
    m1 <- (1:t_len - 0.5 * t_len) * (beta / t_len)
    
    simulated_statistic = replicate(n_rep, {
      sigmahat_vector <- c()
      simulated_data[, 1] <- m1 + arima.sim(model = list(ar = a_hat),
                                            innov = rnorm(t_len, 0, sigma),
                                            n = t_len)
      AR.struc            <- estimate_lrv(data = simulated_data[, 1], q = q,
                                          r_bar = r, p = 1)
      sigma_hat_i         <- sqrt(AR.struc$lrv)
      sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
      for (i in 2:n_ts){
        simulated_data[, i] <- arima.sim(model = list(ar = a_hat),
                                         innov = rnorm(t_len, 0, sigma),
                                         n = t_len)
        AR.struc            <- estimate_lrv(data = simulated_data[, i], q = q,
                                            r_bar = r, p = 1)
        sigma_hat_i         <- sqrt(AR.struc$lrv)
        sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
      }
      
      #Subtracting the column means
      simulated_data <- simulated_data - colMeans(simulated_data)[col(simulated_data)] 
      
      psi     <- compute_statistics(data = simulated_data,
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
      
      cat("Ratio of rejection is ", num_of_rej, "with beta = ", beta,
          ", alpha = ", alpha, "and T = ", t_len, "\n")
    }
    
    #Storing the results in a 3D array
    l <- match(beta, different_beta)
    size_and_power_array[k, l, ] <- size_and_power_vec
  }
}

#######################
#Output of the results#
#######################

for (beta in different_beta){
  l   <- match(beta, different_beta)
  tmp <- size_and_power_array[, l, ]
  if (beta == 0){
    filename = paste0("output/tables/", n_ts, "_ts_size.tex")
  } else {
    filename = paste0("output/tables/", n_ts, "_ts_power_beta_",
                      beta * 100, ".tex")
  }
  output_matrix(tmp, filename)
}