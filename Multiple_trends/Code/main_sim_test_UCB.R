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
library(dplyr)

#Load necessary functions  
source("functions/functions.r")
source("functions/functions_other.r")


##############################
#Defining necessary constants#
##############################

n_ts <- 15 #Number of time series

n_rep    <- 1000 #number of simulations for calculating size and power
sim_runs <- 1000 #number of simulations to calculate the Gaussian quantiles for MS test

different_T <- c(100, 250, 500) #Different lengths of time series
alpha       <- 0.05 #Confidence levels
different_b <- c(0, 0.25, 0.5, 0.75) #Zero is for calculating the size

#For the error process
a     <- 0.25
sigma <- 0.25

#For the covariate process
beta    <- c(1, 1, 1)
a_x_vec <- c(0.25, 0.25, 0.25) #VAR(1) coefficients
phi     <- 0.25                 #dependence between the innovations

#Parameters for the estimation of long-run-variance
q <- 25 
r <- 10

#seed <- 246802468
bw <- 0.2

#For parallel computation
numCores  <- round(parallel::detectCores() * .80)


##################################################
#Cross-validation to obtain the optimal bandwidth#
##################################################

# opt_bw_vec<- c()
# bw <- 0.1
# 
# for (t_len in different_T){
#   m_vec <- bump((1:t_len)/t_len) * 0.5
#   grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
#   
#   for (val in 1:n_rep){
#     #UNIFORM CONFIDENCE BOUNDS
#     estimated_trend_UCB <- c()
#     
#     error    <- arima.sim(model = list(ar = a),
#                           innov = rnorm(t_len, 0, sigma),
#                           n = t_len)
#     nu       <- rmvnorm(t_len + 10, mean = c(0, 0, 0), sigma = phi_matrix)
#     x_matrix <- matrix(0, 3, t_len + 10)
#         
#     for (t in 2:(t_len + 10)){
#       x_matrix[, t] <- a_matrix %*% x_matrix[, t - 1] + nu[t, ]
#     }
#     x_matrix <- t(x_matrix[, -(1:10)])
#         
#     y <- m_vec + beta %*% t(x_matrix) + error[, i]
#         
#     #First differences
#     y_diff_tmp <- y - dplyr::lag(y, n = 1, default = NA)
#     y_diff     <- as.matrix(y_diff_tmp)[-1, ]
#     x_diff_1   <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
#     x_diff_2   <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
#     x_diff_3   <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)
#         
#     #Estimating beta
#     x_diff    <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]
#     beta_hat  <- solve(t(x_diff) %*% x_diff) %*% t(x_diff) %*% y_diff
# 
#     y_augm <- y - x_matrix %*% as.vector(beta_hat)
#         
#     estimated_trend_UCB <- mapply(UCB_estimation, grid_points,
#                                   MoreArgs = list(data_p = y_augm,
#                                                   grid_p = grid_points,
#                                                   bw = bw))
#     y_fitted <- x_matrix %*% as.vector(beta_hat) + estimated_trend_UCB
#     
#     h_matrix <- matrix(NA, ncol = t_len, nrow = t_len)
#     w_matrix <- matrix(NA, ncol = t_len, nrow = t_len)
#     for (t in 1:t_len){
#       for (s in 1:t_len){
#         s_t_2_value1 = s_t_2_UCB(x = s/t_len, h = bw, T_size = t_len, x_vec = grid_points)
#         s_t_1_value1 = s_t_1_UCB(x = s/t_len, h = bw, T_size = t_len, x_vec = grid_points)
#         s_t_0_value1 = s_t_0_UCB(x = s/t_len, h = bw, T_size = t_len, x_vec = grid_points)
#         num1 = s_t_2_value1 * s_t_0_value1 - s_t_1_value1^2
#         s_t_2_value2 = s_t_2_UCB(x = s/t_len, h = bw * sqrt(2), T_size = t_len, x_vec = grid_points)
#         s_t_1_value2 = s_t_1_UCB(x = s/t_len, h = bw * sqrt(2), T_size = t_len, x_vec = grid_points)
#         s_t_0_value2 = s_t_0_UCB(x = s/t_len, h = bw * sqrt(2), T_size = t_len, x_vec = grid_points)
#         num2 = s_t_2_value2 * s_t_0_value2 - s_t_1_value2^2
#         u = t/t_len - s/t_len
#         denom1 = (s_t_2_value1 - s_t_1_value1 * u) * epanechnikov_kernel(u / bw)
#         denom2 = (s_t_2_value2 - s_t_1_value2 * u) * epanechnikov_kernel(u / (bw * sqrt(2)))
#         w_matrix[t, s] = 2 * denom1 / num1 - denom2 / num2
#         h_matrix[t, s] <- w_matrix[t, s] + t(x_matrix[t, ]) %*% solve(t(x_diff) %*% x_diff)        
#       }
#     }
#   }
#   opt_bw_vec <- c(opt_bw_vec, opt_bw)
# }

#######################
#Plotting one instance#
#######################

different_T_plotting <- c(500)

for (t_len in different_T_plotting){
  #Calculating the Gaussian quantiles

  cat("Calculating the Gaussian quantiles\n")

  #Calculating the Gaussian quantiles for UCB in parallel
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:sim_runs, .combine = "cbind") %dopar% {
    source("functions/functions_other.R")
    repl_UCB(rep_ = val, n_ts_ = n_ts, t_len_ = t_len, bw_ = bw,
             gaussian_sim = TRUE)
    # Loop one-by-one using foreach
  } -> simulated_gaussian_UCB
  stopCluster(cl)
  toc()

  probs <- seq(0.5, 0.995, by = 0.005)

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

  quant_UCB <- quants_UCB[pos]

  #Calculating UCB

  grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)

  k_n            <- floor(t_len^(1/3))
  m              <- floor(t_len / k_n)
  lower_boundary <- ceiling(0.05 * t_len)
  upper_boundary <- floor(0.95 * t_len)
  grid           <- grid_points[lower_boundary:upper_boundary] #For plotting

  m_matrices          <- list()
  y_matrices          <- list()
  y_augm_matrices     <- list()
  estimated_trend_UCB <- list()
  upper_UCB           <- list()
  lower_UCB           <- list()

  for (k in 1:length(different_b)){
    y_matrices[[k]]          <- matrix(NA, nrow = t_len, ncol = n_ts)
    m_matrices[[k]]          <- matrix(0, nrow = t_len, ncol = n_ts)
    m_matrices[[k]][, 1]     <- bump((1:t_len)/t_len) * different_b[k]
    y_augm_matrices[[k]]     <- matrix(NA, nrow = t_len, ncol = n_ts)
    estimated_trend_UCB[[k]] <- matrix(NA, nrow = t_len, ncol = n_ts)
    upper_UCB[[k]]           <- matrix(NA, nrow = t_len, ncol = n_ts)
    lower_UCB[[k]]           <- matrix(NA, nrow = t_len, ncol = n_ts)
  }

  error_matrix  <- matrix(NA, nrow = t_len, ncol = n_ts)

  phi_matrix       <- matrix(phi, nrow = 3, ncol = 3)
  diag(phi_matrix) <- 1
  a_matrix         <- diag(a_x_vec)

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

    x_diff_1 <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
    x_diff_2 <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
    x_diff_3 <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)

    #Estimating beta
    x_diff_tmp <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]

    k <- 1
    for (b in different_b){
      y_matrices[[k]][, i] <- m_matrices[[k]][, i] + beta %*% t(x_matrix) + error_matrix[, i]

      #First differences
      y_diff     <- y_matrices[[k]][, i] - dplyr::lag(y_matrices[[k]][, i], n = 1, default = NA)
      y_diff_tmp <- as.matrix(y_diff)[-1, ]

      beta_hat_tmp <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp

      y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp)

      sigma_hat_UCB_i <- sqrt(sigma_estimation_UCB(y_ = y_matrices[[k]][, i],
                                                   x_matrix_ = x_matrix,
                                                   beta_est_ = beta_hat_tmp,
                                                   m_ = m, k_n_ = k_n))

      estimated_trend_UCB[[k]][, i] <- mapply(UCB_estimation, grid_points,
                                              MoreArgs = list(data_p = y_augm_matrices[[k]][, i],
                                                              grid_p = grid_points,
                                                              bw = bw))
      upper_UCB[[k]][, i] <- estimated_trend_UCB[[k]][, i] + sigma_hat_UCB_i * quant_UCB
      lower_UCB[[k]][, i] <- estimated_trend_UCB[[k]][, i] - sigma_hat_UCB_i * quant_UCB

      k <- k + 1
    }
  }

  for (k in 1:length(different_b)){

    #Ignoring the boundary issues
    tmp_u <- upper_UCB[[k]][lower_boundary:upper_boundary, ]
    tmp_l <- lower_UCB[[k]][lower_boundary:upper_boundary, ]

    filename = paste0("output/revision/UCB_plot_with_T_", t_len, "_and_b_", different_b[k]*100, ".pdf")
    pdf(filename, width = 5, height = 3.5, paper="special")

    #Setting the layout of the graphs
    par(cex = 1, tck = -0.025)
    par(mar = c(0, 0, 0, 0)) #Margins for each plot
    par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

    plot(x = grid, xlim = c(0, 1),
         y = estimated_trend_UCB[[k]][lower_boundary:upper_boundary, 1], type = 'l', col = 'red', ylim = c(-1.2, 1.2),
         xlab = "", ylab = "", main = NULL, cex = 0.8)
    lines(x = grid,
          y = upper_UCB[[k]][lower_boundary:upper_boundary, 1], type = "l",
          col = "red", lty = 2)
    lines(x = grid,
          y = lower_UCB[[k]][lower_boundary:upper_boundary, 1], type = "l",
          col = "red", lty = 2)
    lines(x = grid,
          y = estimated_trend_UCB[[k]][lower_boundary:upper_boundary, 2], type = 'l', col = 'blue')
    lines(x = grid,
          y = upper_UCB[[k]][lower_boundary:upper_boundary, 2], type = "l",
          col = "blue", lty = 2)
    lines(x = grid,
          y = lower_UCB[[k]][lower_boundary:upper_boundary, 2], type = "l",
          col = "blue", lty = 2)
    dev.off()
  }
}


################################
#Calculating the size and power#
################################

size_and_power_UCB_array <- array(NA, dim = c(length(different_T),
                                              length(different_b),
                                              1),
                                  dimnames = list(t = different_T,
                                                  b = different_b,
                                                  alpha = alpha))

for (t_len in different_T){
#  set.seed(seed)
  k   <- match(t_len, different_T)

  ####################################
  #Calculating the Gaussian quantiles#
  ####################################
  
  cat("Calculating the Gaussian quantiles\n")

  #Calculating the Gaussian quantiles for UCB in parallel
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:sim_runs, .combine = "cbind") %dopar% {
    source("functions/functions_other.R")
    repl_UCB(rep_ = val, n_ts_ = n_ts, t_len_ = t_len, bw_ = 0.1,
             gaussian_sim = TRUE)
    # Loop one-by-one using foreach
  } -> simulated_gaussian_UCB
  stopCluster(cl)
  toc()

  probs <- seq(0.5, 0.995, by = 0.005)

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

  quant_UCB <- quants_UCB[pos]
  
    
  #################################
  #Testing for different scenarios#
  #################################
  
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% {
    source("functions/functions.R")
    repl_UCB(rep_ = val, n_ts_ = n_ts, t_len_ = t_len, bw_ = 0.1,
         a_ = a, sigma_ = sigma,
         beta_ = beta, a_x_vec_ = a_x_vec, phi_ = phi,
         different_b_ = different_b, quant_UCB_ = quant_UCB,
         gaussian_sim = FALSE)
    # Loop one-by-one using foreach
  } -> pairwise_comparison_UCB
  stopCluster(cl)
  toc()
  
  for (j in 1:length(different_b)){
    pairwise_results_UCB <- apply(pairwise_comparison_UCB[((j - 1) * n_ts * n_ts + 1):(j * n_ts * n_ts), ], 2, sum)
    num_of_rej_UCB       <- sum(pairwise_results_UCB != 0)/n_rep

    cat("Ratio of rejection for UCB is ", num_of_rej_UCB, "with b = ", different_b[j],
        ", alpha = ", alpha, "and T = ", t_len, "\n")

    #Storing the results in a 3D array
    size_and_power_UCB_array[k, j, ] <- num_of_rej_UCB
  }
} 

save(size_and_power_UCB_array, file = "output/revision/UCB_simulations.R")
#load(file = "output/revision/UCB_simulations.R")


#######################
#Output of the results#
#######################

tmp <- as.matrix(size_and_power_UCB_array[, , 1])
row.names(tmp) <- paste0("$T = ", row.names(as.matrix(size_and_power_UCB_array[ , , 1])), "$")

  
filename = paste0("output/revision/", n_ts, "_ts_UCB.tex")

#Create a matrix (for size and power table for example) and write them in the tex file
addtorow     <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{4}{c}{different bump height $b$} \\\\\n",
                      "$T$ & 0 & 0.25 & 0.5 & 0.75 \\\\\n") 
print.xtable(xtable(tmp, digits = c(3), align = "ccccc"), type = "latex",
             file = filename, add.to.row = addtorow, include.colnames = FALSE,
             sanitize.text.function=function(x){x})
line <- paste0("%This simulation was done for the following values of the parameters: n_ts = ", n_ts,
               ", with ", n_rep, " simulations for calculating size and power and ", sim_runs,
               " simulations to calculate the Gaussian quantiles. Furthermore, for the error process we have a = ",
               a, " and sigma = ", sigma,
               ". There are no fixed effects. The grid is normal.")
write(line, file = filename, append = TRUE)