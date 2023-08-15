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

##############################
#Defining necessary constants#
##############################

n_rep    <- 1000 #number of simulations for calculating size and power
sim_runs <- 1000 #number of simulations to calculate the Gaussian quantiles

n_ts <- 15

different_T     <- c(100, 250, 500) #Different lengths of time series
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
different_b     <- c(0) #Zero is for calculating the size

# different_b     <- c(0, 0.75, 1.00, 1.25) #Zero is for calculating the size

#For the covariate process
beta        <- c(1, 1, 1)
a_x_vec     <- c(0.5, 0.5, 0.5)
phi         <- 0.1

#For the error process
a <- 0.25

#For the fixed effects
rho <- 0.1 

#Parameters for the estimation of long-run-variance
q <- 15 
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
  k <- match(t_len, different_T)

  #Constructing the grid according to the application
  u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
  h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
  h_grid <- h_grid[h_grid > log(t_len) / t_len]
  grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
  

  #Calculating the Gaussian quantiles in parallel
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:sim_runs, .combine = "cbind") %dopar% {
    repl_revision2(rep_ = val, n_ts_ = n_ts, t_len_ = t_len, grid_ = grid,
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
    #Only the first trend function is non-zero:
    m_matrix[, 1] <- (1:t_len - 0.5 * t_len) * (b / t_len)

    tic()
    
    cl <- makePSOCKcluster(numCores)
    registerDoParallel(cl)
    foreach (val = 1:n_rep, .combine = "cbind") %dopar% {
      repl_revision2(rep_ = val, n_ts_ = n_ts, t_len_ = t_len, grid_ = grid,
                     a_ = a, beta_ = beta, a_x_vec = a_x_vec, phi_ = phi,
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
  if (b == 0){
    filename = paste0("output/revision/", n_ts, "_ts_size.tex")
  } else {
    filename = paste0("output/revision/", n_ts, "_ts_power_b_",
                      b * 100, ".tex")
  }
  output_matrix(tmp, filename)
}