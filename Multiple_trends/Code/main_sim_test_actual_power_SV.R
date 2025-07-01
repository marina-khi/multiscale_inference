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

source("functions/functions.R")

##############################
#Defining necessary constants#
##############################
seed <- 543212345

n_ts <- 15 #Number of time series

#For the covariate process
beta    <- c(1, 1, 1)
a_x_vec <- c(0.25, 0.25, 0.25) #VAR(1) coefficients
phi     <- 0.25                 #dependence between the innovations

#For the error process
a     <- 0.25
#sigma <- 0.25

#For the fixed effects
rho      <- 0.25 #covariance between the fixed effects
n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles

#Different parameters
different_T     <- c(100, 250, 500) #Different lengths of time series  
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
different_b     <- c(0.25, 0.5, 0.75) #Zero is for calculating the size

#Parameters for the estimation of long-run-variance
q <- 25 
r <- 10

#For parallel computation
numCores  <- round(parallel::detectCores() * .80)

#Calculating actual power
u.lower1 <- 0.2
u.upper1 <- 0.4
u.lower2 <- 0.6
u.upper2 <- 0.8

##################################################
#Calculating the size and power for a normal grid#
##################################################

actual_power_array_SV <- array(NA, dim = c(length(different_T),
                                        length(different_b),
                                        length(different_alpha)),
                            dimnames = list(t = different_T,
                                            b = different_b,
                                            alpha = different_alpha))

majority_power_array_SV <- array(NA, dim = c(length(different_T),
                                          length(different_b),
                                          length(different_alpha)),
                              dimnames = list(t = different_T,
                                              b = different_b,
                                              alpha = different_alpha))

full_power_array_SV <- array(NA, dim = c(length(different_T),
                                      length(different_b),
                                      length(different_alpha)),
                          dimnames = list(t = different_T,
                                          b = different_b,
                                          alpha = different_alpha))

for (t_len in different_T){
  set.seed(seed)
  k <- match(t_len, different_T)
  #Constructing the full grid for calculating the Gaussian quantiles
  u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
  
  #Calculating the Gaussian quantiles in parallel
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:sim_runs, .combine = "cbind") %dopar% {
    source("functions/functions.R")
    repl_SV(rep_ = val, n_ts_ = n_ts, t_len_ = t_len, grid_ = grid,
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
  
  #Restricting the grid to only look at the actual power
  gset_pos      <- grid$gset
  deletions     <- (((u.lower1 <= gset_pos$u + gset_pos$h) & (gset_pos$u - gset_pos$h <= u.upper1)) | ((u.lower2 <= gset_pos$u + gset_pos$h) & (gset_pos$u - gset_pos$h <= u.upper2)))
  grid_actual   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid, deletions = deletions)
  
  #Calculating the true test statistics
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% {
    source("functions/functions.R")
    repl_SV(rep_ = val, n_ts_ = n_ts, t_len_ = t_len,
         grid_ = grid_actual, #ijset_ = ijset,
         a_ = a, beta_ = beta, a_x_vec_ = a_x_vec, phi_ = phi, rho_ = rho,
         different_b_ = different_b,
         q_ = q, r_ = r)
    # Loop one-by-one using foreach
  } -> simulated_pairwise_statistics
  stopCluster(cl)
  toc()
  
  for (j in 1:length(different_b)){
    statistic_values <- simulated_pairwise_statistics[((j - 1) * n_ts * n_ts + 1):(j * n_ts * n_ts), ]
    
    actual_power_vec   <- c()
    majority_power_vec <- c()
    full_power_vec     <- c()
    
    for (alpha in different_alpha){
      if (sum(probs == (1 - alpha)) == 0)
        pos <- which.min(abs(probs - (1 - alpha)))
      if (sum(probs == (1 - alpha)) != 0)
        pos <- which.max(probs == (1 - alpha))    
      quant <- quants[pos]
      
      num_of_actual_rej <- 0
      num_of_majority_rej <- 0
      num_of_full_rej <- 0
      
      for (val in 1:n_rep){
        tmp         <- matrix(statistic_values[, val], nrow = n_ts, ncol = n_ts)
        num_of_rej  <- sum(tmp[1, ] > quant)
        if (num_of_rej > 0)   {num_of_actual_rej   <- num_of_actual_rej + 1}
        if (num_of_rej > 6)   {num_of_majority_rej <- num_of_majority_rej + 1}
        if (num_of_rej == 14) {num_of_full_rej     <- num_of_full_rej + 1}
      }
      actual_power_vec <- c(actual_power_vec, num_of_actual_rej/n_rep)
      majority_power_vec <- c(majority_power_vec, num_of_majority_rej/n_rep)
      full_power_vec <- c(full_power_vec, num_of_full_rej/n_rep)
      
      cat("Ratio of correct rejections in at least one case is ",
          num_of_actual_rej/n_rep, "with b = ", different_b[j],
          ", alpha = ", alpha, "and T = ", t_len, "\n")
      cat("Ratio of correct rejections in majority of the cases is ",
          num_of_majority_rej/n_rep, "with b = ", different_b[j],
          ", alpha = ", alpha, "and T = ", t_len, "\n")
      cat("Ratio of correct rejections in all cases is ",
          num_of_full_rej/n_rep, "with b = ", different_b[j],
          ", alpha = ", alpha, "and T = ", t_len, "\n")
    }
    
    #Storing the results in a 3D array
    actual_power_array_SV[k, j, ] <- actual_power_vec
    majority_power_array_SV[k, j, ] <- majority_power_vec
    full_power_array_SV[k, j, ] <- full_power_vec
  }
}


#Output of the results
for (b in different_b){
  l   <- match(b, different_b)
  tmp <- as.matrix(actual_power_array_SV[, l, ])
  filename = paste0("output/revision/", n_ts, "_ts_", phi*100, "_", rho * 100, "_actual_power_b_",
                    b * 100, "_SV.tex")
  output_matrix(tmp, filename, numcols_ = 4)
  line <- paste0("%This simulation was done for the seed ", seed,
                 ", for the following values of the parameters: n_ts = ", n_ts,
                 ", with ", n_rep, " simulations for calculating actual power and ", sim_runs,
                 " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
                 a, " and sigma = ", sigma, 
                 ". For the covariate process a_1 = a_2 = a_3 = ", a_x_vec[1], " and phi = ", phi,
                 ". For the fixed effect, we have rho = ", rho,
                 ". The grid is normal")     
  write(line, file = filename, append = TRUE)
}

for (b in different_b){
  l   <- match(b, different_b)
  tmp <- as.matrix(majority_power_array_SV[, l, ])
  filename = paste0("output/revision/", n_ts, "_ts_", phi*100, "_", rho * 100, "_majority_power_b_",
                    b * 100, "_SV.tex")
  output_matrix(tmp, filename, numcols_ = 4)
  line <- paste0("%This simulation was done for the seed ", seed,
                 ", for the following values of the parameters: n_ts = ", n_ts,
                 ", with ", n_rep, " simulations for calculating majority power and ", sim_runs,
                 " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
                 a, " and sigma = ", sigma, 
                 ". For the covariate process a_1 = a_2 = a_3 = ", a_x_vec[1], " and phi = ", phi,
                 ". For the fixed effect, we have rho = ", rho,
                 ". The grid is normal")     
  write(line, file = filename, append = TRUE)
}

for (b in different_b){
  l   <- match(b, different_b)
  tmp <- as.matrix(full_power_array_SV[, l, ])
  filename = paste0("output/revision/", n_ts, "_ts_", phi*100, "_", rho * 100, "_full_power_b_",
                    b * 100, "_SV.tex")
  output_matrix(tmp, filename, numcols_ = 4)
  line <- paste0("%This simulation was done for the seed ", seed,
                 ", for the following values of the parameters: n_ts = ", n_ts,
                 ", with ", n_rep, " simulations for calculating full power and ", sim_runs,
                 " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
                 a, " and sigma = ", sigma, 
                 ". For the covariate process a_1 = a_2 = a_3 = ", a_x_vec[1], " and phi = ", phi,
                 ". For the fixed effect, we have rho = ", rho,
                 ". The grid is normal")     
  write(line, file = filename, append = TRUE)
}