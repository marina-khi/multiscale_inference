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
n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles

different_T     <- c(250, 500, 1000) #Different lengths of time series
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
len             <- length(different_alpha) #Necessary for output of the results
different_beta  <- c(0.75, 1, 1.25)

a_hat <- 0.5 
sigma <- 0.5
q     <- 25 #Parameters for the estimation of long-run-variance
r     <- 10

######################
#Calculating the size#
######################

size_vec  <- c()

#Constructing the set of pairwise comparisons
ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
ijset <- ijset[ijset$i < ijset$j, ]

for (t_len in different_T){
  simulated_data           <- matrix(NA, nrow = t_len, ncol = n_ts)
  colnames(simulated_data) <- 1:n_ts
  
  #Constructing the grid
  u_grid <- seq(from = 5 / t_len, to = 1, by = 1 / 50)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len)
  
  simulated_statistic = replicate(n_rep, {
    sigmahat_vector <- c()
    for (i in 1:n_ts){
      simulated_data[, i] <- arima.sim(model = list(ar = a_hat),
                                       innov = rnorm(t_len, 0, sigma),
                                       n = t_len)
      simulated_data[, i] <- simulated_data[, i] - mean(simulated_data[, i])
      AR.struc            <- estimate_lrv(data = simulated_data[, i], q = q,
                                          r_bar = r, p = 1)
      sigma_hat_i         <- sqrt(AR.struc$lrv)
      sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
    }
    psi     <- compute_statistics(data = simulated_data,
                                  sigma_vec = sigmahat_vector,
                                  n_ts = n_ts, grid = grid)
    results <- max(psi$stat_pairwise)
    results
  })
  
  #Calculating Gaussian quantiles
  quantiles <- compute_quantiles(t_len = t_len, grid = grid, n_ts = n_ts,
                                 ijset = ijset, sim_runs = sim_runs,
                                 correction = TRUE)
  probs  <- as.vector(quantiles$quant[1, ])
  quants <- as.vector(quantiles$quant[2, ])
  
 
}


#######################
#Output of the results#
#######################

filename = paste0("output/tables/", n_ts, "_ts_size.tex")
creating_matrix_and_texing(size_vec, different_T, different_alpha, filename)


#######################
#Calculating the power#
#######################

power_matrix  <- matrix(NA, nrow = length(different_T),
                        ncol = len * length(different_beta))

#Constructing the set of pairwise comparisons
ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
ijset <- ijset[ijset$i < ijset$j, ]

for (t_len in different_T){
  #Constructing the grid
  u_grid      <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
  h_grid      <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid      <- h_grid[h_grid > log(t_len) / t_len]
  grid  <- construct_grid(t = t_len)
  
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
    
    power_vec <- c()
    for (alpha in different_alpha){
      if (sum(probs == (1 - alpha)) == 0)
        pos <- which.min(abs(probs - (1 - alpha)))
      if (sum(probs == (1 - alpha)) != 0)
        pos <- which.max(probs == (1 - alpha))    
      quant <- quants[pos]
      
      power     <- sum(simulated_statistic > quant)/n_rep
      power_vec <- c(power_vec, power) 
      
      cat("Ratio of rejection under H1 is ", power, "with beta = ", beta,
          ", alpha = ", alpha, "and T = ", t_len, "\n")
    }
    
    #Storing the results in a matrix
    k <- match(t_len, different_T)
    l <- match(beta, different_beta)
    power_matrix[k, (len * (l - 1) + 1):(len * l)] <- power_vec
  }
}

#######################
#Output of the results#
#######################

for (beta in different_beta){
  l <- match(beta, different_beta)
  tmp <- power_matrix[, (len * (l - 1) + 1):(len * l)]
  filename = paste0("output/tables/", n_ts, "_ts_power_beta_",
                    beta * 10, ".tex")
  creating_matrix_and_texing(tmp, different_T, different_alpha, filename)
}