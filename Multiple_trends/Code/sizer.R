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

#Load necessary functions  
source("functions/functions.r")
# source("functions/sim.r")
sourceCpp("functions/SiZer_functions.cpp")


##############################
#Defining necessary constants#
##############################

n_ts  <- 2 #Number of time series

n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles for MS test

different_T     <- c(100, 250, 500, 750, 1000) #Different lengths of time series
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
different_b     <- c(0, 0.25, 0.5) #Zero is for calculating the size

#For the error process
a            <- 0.25
sigma        <- 0.25
lrv          <- sigma^2/((1 - a)^2)
sigma_vector <- rep(sqrt(lrv), n_ts)

seed <- 246802468

################################
#Calculating the size and power#
################################

size_and_power_array <- array(NA, dim = c(length(different_T),
                                          length(different_b),
                                          length(different_alpha)),
                              dimnames = list(t = different_T,
                                              b = different_b,
                                              alpha = different_alpha))
size_and_power_sizer_array <- array(NA, dim = c(length(different_T),
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
  k <- match(t_len, different_T)
  
  #Constructing the grid
  u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

  
  #################################
  #Calculating the SiZer quantiles#
  #################################
  
  cat("Calculating the SiZer quantiles\n")
  autocov <- (sigma^2/(1 - a^2)) * (a^seq(0, t_len - 1, by = 1))
  
  sizer.quants <- vector("list", length(different_alpha))
  for (j in 1:length(different_alpha)){
    sizer.quants[[j]] <- SiZer_quantiles(alpha = different_alpha[j],
                                         t_len = t_len,
                                         grid = grid, autocov1 = autocov,
                                         autocov2 = autocov)
  }
  sizer.wghts <- SiZer_weights(t_len = t_len, grid = grid)
  sizer.std    <- SiZer_std(weights = sizer.wghts, autocov1 = autocov,
                            autocov2 = autocov, t_len = t_len)
  #sizer.std1  <- SiZer_std(sigma = sigma, weights_ = sizer.wghts, autocov = autocov, t_len = t_len)
  #sizer.std   <- sqrt(2) * sizer.std1
  
    
  ####################################
  #Calculating the Gaussian quantiles#
  ####################################
  
  cat("Calculating the Gaussian quantiles\n")
  simulated_pairwise_gaussian <- matrix(NA, nrow = n_ts * n_ts, ncol = sim_runs)

  for (val in 1:sim_runs){
    z_matrix      <- matrix(NA, nrow = t_len, ncol = n_ts)
    z_augm_matrix <- matrix(NA, nrow = t_len, ncol = n_ts)
    
    for (i in 1:n_ts){
      z_matrix[, i]      <- rnorm(t_len, 0, sqrt(lrv))
      z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    }
    
    psi <- compute_statistics(data = z_augm_matrix,
                              sigma_vec = sigma_vector,
                              n_ts = n_ts, grid = grid)
    simulated_pairwise_gaussian[, val] <- as.vector(psi$stat_pairwise)
  }

  simulated_gaussian <- apply(simulated_pairwise_gaussian, 2, max)
  
  probs      <- seq(0.5, 0.995, by = 0.005)
  quantiles  <- as.vector(quantile(simulated_gaussian, probs = probs))
  quantiles  <- rbind(probs, quantiles)
  
  colnames(quantiles) <- NULL
  rownames(quantiles) <- NULL
  
  quants <- as.vector(quantiles[2, ])

  #################################
  #Testing for different scenarios#
  #################################
  
  for (b in different_b){
    simulated_pairwise_statistics <- matrix(NA, nrow = n_ts * n_ts, ncol = n_rep)
    sizer_results_matrix          <- matrix(NA, nrow = length(different_alpha), ncol = n_rep)
        
    m_matrix      <- matrix(0, nrow = t_len, ncol = n_ts)
    if (b == 0) {
      cat("SIZE SIMULATIONS\n")
    } else {
      cat("POWER SIMULATIONS WITH b = ", b, "\n")
      #Only the first trend function is non-zero:
      m_matrix[, 1] <- bump((1:t_len)/t_len) * b
    }
    
    for (val in 1:n_rep){
      #Simulated data
      y_matrix      <- matrix(NA, nrow = t_len, ncol = n_ts)
      y_augm_matrix <- matrix(NA, nrow = t_len, ncol = n_ts)
      error_matrix  <- matrix(NA, nrow = t_len, ncol = n_ts)
      
      y_sizer_matrix <- matrix(NA, nrow = t_len * n_ts, ncol = 2)
    
      for (i in 1:n_ts){
        error_matrix[, i] <- arima.sim(model = list(ar = a),
                                       innov = rnorm(t_len, 0, sigma),
                                       n = t_len)
        
        y_matrix[, i]     <- m_matrix[, i] + error_matrix[, i]
        
        #Estimating the fixed effects
        alpha_hat_tmp      <- mean(y_matrix[, i])
        y_augm_matrix[, i] <- y_matrix[, i] - alpha_hat_tmp
      }
      
      #MULTISCALE TEST
      psi <- compute_statistics(data = y_augm_matrix,
                                sigma_vec = sigma_vector,
                                n_ts = n_ts, grid = grid)    
      simulated_pairwise_statistics[, val] <- as.vector(psi$stat_pairwise)
      
      #VALUES FOR SIZER TEST
      values1     <- sizer.wghts %*% y_matrix[, 1]
      sizer.vals1 <- as.vector(values1)
      values2     <- sizer.wghts %*% y_matrix[, 2]
      sizer.vals2 <- as.vector(values2)
      
      sizer_results_vec <- c()
      for (alpha in different_alpha){
        j <- match(alpha, different_alpha)
        SiZer_results <- SiZer_test(values1 = sizer.vals1, values2 = sizer.vals2,
                                    std.devs = sizer.std, quants = sizer.quants[[j]],
                                    grid = grid)
        sizer_results_vec <- c(sizer_results_vec, as.integer(sum(abs(SiZer_results$test))>0))
        
        #SiZermap(SiZer_results$ugrid, SiZer_results$hgrid, SiZer_results$test, plot.title = NA)
      }
      sizer_results_matrix[, val] <- as.vector(sizer_results_vec)
    }
    
    simulated_statistic <- apply(simulated_pairwise_statistics[1:(n_ts * n_ts), ], 2, max)
    
    size_and_power_vec       <- c()
    size_and_power_sizer_vec <- c()
    
    for (alpha in different_alpha){
      j <- match(alpha, different_alpha)
      if (sum(probs == (1 - alpha)) == 0)
        pos <- which.min(abs(probs - (1 - alpha)))
      if (sum(probs == (1 - alpha)) != 0)
        pos <- which.max(probs == (1 - alpha))    
      quant <- quants[pos]
      
      num_of_rej         <- sum(simulated_statistic > quant)/n_rep
      size_and_power_vec <- c(size_and_power_vec, num_of_rej) 
      
      cat("Ratio of rejection is ", num_of_rej, "with b = ", b,
          ", alpha = ", alpha, "and T = ", t_len, "\n")
      
      num_of_rej_sizer         <- sum(sizer_results_matrix[j, ])/n_rep
      size_and_power_sizer_vec <- c(size_and_power_sizer_vec, num_of_rej_sizer) 
      
      cat("Ratio of rejection for SiZer is ", num_of_rej_sizer, "with b = ", b,
          ", alpha = ", alpha, "and T = ", t_len, "\n")
    }
    
    #Storing the results in a 3D array
    l <- match(b, different_b)
    size_and_power_array[k, l, ]       <- size_and_power_vec
    size_and_power_sizer_array[k, l, ] <- size_and_power_sizer_vec
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
    filename = paste0("output/revision/", n_ts, "_ts_size_sizer_comparison.tex")
  } else {
    filename = paste0("output/revision/", n_ts, "_ts_power_b_",
                      b * 100, "_sizer_comparison.tex")
  }
  output_matrix2(tmp, filename)
  line <- paste0("%This simulation was done for the following values of the parameters: n_ts = ", n_ts,
                 ", with ", n_rep, " simulations for calculating size and power and ", sim_runs,
                 " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
                 a, " and sigma = ", sigma,
                 ". There are no fixed effects. The grid is standard.")
  write(line, file = filename, append = TRUE)
}