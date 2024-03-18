rm(list=ls())

library(MSinference)
library(haven)
library(car)
library(dplyr)
library(Matrix)
library(foreach)
library(parallel)
library(doParallel)
library(xtable)
library(tictoc)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


##############################
#Defining necessary constants#
##############################


n_ts <- 15 #Number of time series

#For the covariate process
#beta    <- c(1, 1, 1)
beta    <- NULL
a_x_vec <- c(0.25, 0.25, 0.25) #VAR(1) coefficients
phi     <- 0.1                 #dependence between the innovations

#For the error process
a     <- 0.25
sigma <- 0.25

#For the fixed effects
rho      <- 0.1 #covariance between the fixed effects
n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles

if (is.null(beta)){
  different_T     <- c(100, 250, 500, 750, 1000) #Different lengths of time series
  seed <- 113355
} else {
  different_T     <- c(100, 250, 500) #Different lengths of time series  
  seed <- 111333555
}

different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
different_b     <- c(0) #Zero is for calculating the size


#Parameters for the estimation of long-run-variance
q <- 25 
r <- 10

#For parallel computation
numCores  <- round(parallel::detectCores() * .70)



################################
#Calculating the size and power#
################################

source("functions/functions.R")

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
  set.seed(seed)
  #set.seed(13579135)
  k <- match(t_len, different_T)

  #Constructing the grid according to the application
  
  if (is.null(beta)){
    u_grid <- seq(from = 10 / t_len, to = 1, by = 10 / t_len)
    h_grid <- seq(from = 7 / t_len, to = 1 / 4, by = 10 / t_len)
  } else {
    u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
    h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  }
  
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

  #Calculating the Gaussian quantiles in parallel
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:sim_runs, .combine = "cbind") %dopar% {
    repl_revision(rep_ = val, n_ts_ = n_ts, t_len_ = t_len, grid_ = grid,
                   gaussian_sim = TRUE)
    # Loop one-by-one using foreach
  } -> simulated_pairwise_gaussian
  stopCluster(cl)
  toc()

  simulated_gaussian <- apply(simulated_pairwise_gaussian, 2, max)
  
  probs      <- seq(0.5, 0.995, by = 0.005)
  quantiles  <- as.vector(quantile(simulated_gaussian, probs = probs))
  quantiles  <- rbind(probs, quantiles)
  
  colnames(quantiles) <- NULL
  rownames(quantiles) <- NULL
  
  quants <- as.vector(quantiles[2, ])
  
  for (b in different_b){
    
    m_matrix      <- matrix(0, nrow = t_len, ncol = n_ts)
    if (b == 0) {
      cat("SIZE SIMULATIONS\n")
    } else {
      cat("POWER SIMULATIONS WITH b = ", b, "\n")
      
      #Only the first trend function is non-zero:
      #m_matrix[, 1] <- (1:t_len - 0.5 * t_len) * (b / t_len)
      #trend function which is a bump
      m_matrix[, 1] <- bump((1:t_len)/t_len) * b
    }
    
    tic()
    cl <- makePSOCKcluster(numCores)
    registerDoParallel(cl)
    foreach (val = 1:n_rep, .combine = "cbind") %dopar% {
      repl_revision(rep_ = val, n_ts_ = n_ts, t_len_ = t_len, grid_ = grid,
                     a_ = a, sigma_ = sigma,
                     beta_ = beta,
                     a_x_vec_ = a_x_vec, phi_ = phi,
                     rho_ = rho, m_matrix_ = m_matrix, q_ = q, r_ = r)
      # Loop one-by-one using foreach
    } -> simulated_pairwise_statistics
    stopCluster(cl)
    toc()
        
    simulated_statistic <- apply(simulated_pairwise_statistics[1:(n_ts * n_ts), ], 2, max)
    
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
  if (is.null(beta)){
    if (b == 0){
      filename = paste0("output/revision/", n_ts, "_ts_no_covariates_", rho * 100, "_size.tex")
    } else {
      filename = paste0("output/revision/", n_ts, "_ts_no_covariates_", rho * 100, "_power_b_",
                        b * 100, ".tex")
    }
  } else {
    if (b == 0){
      filename = paste0("output/revision/", n_ts, "_ts_", phi*100, "_", rho * 100, "_size.tex")
    } else {
      filename = paste0("output/revision/", n_ts, "_ts_", phi*100, "_", rho * 100, "_power_b_",
                        b * 100, ".tex")
    }
  }
  output_matrix(tmp, filename)
  if (is.null(beta)){
    line <- paste0("%This simulation was done for the seed ", seed,
                   ", for the following values of the parameters: n_ts = ", n_ts,
                   ", with ", n_rep, " simulations for calculating size and power and ", sim_runs,
                   " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
                   a, " and sigma = ", sigma,
                   ". We do not have any covarites. For the fixed effect, we have rho = ", rho,
                   ". The grid is fine (growing with the sample size)")    
  } else {
    line <- paste0("%This simulation was done for the seed ", seed,
                   ", for the following values of the parameters: n_ts = ", n_ts,
                   ", with ", n_rep, " simulations for calculating size and power and ", sim_runs,
                   " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
                   a, " and sigma = ", sigma, 
                   ". For the covariate process a_1 = a_2 = a_3 = ", a_x_vec[1], " and phi = ", phi,
                   ". For the fixed effect, we have rho = ", rho,
                   ". The grid is fine (growing with the sample size)")     
  }
  write(line, file = filename, append = TRUE)
}