rm(list=ls())

library(readxl)
library(haven)
library(tidyr)
library(multiscale)
library(dplyr)
library(zoo)

source("functions/functions.R")


######################
#Necessary parameters#
######################

alpha     <- 0.05
sim_runs  <- 5000
q         <- 15 #Parameters for the estimation of sigma
r         <- 10


##############
#Data loading#
##############

y <- read.csv(file = "misc/Data/volatility_start_change2007.csv")
tmp <- read.csv(file = "misc/Data/augmented_volatility.csv")
y_augm <- as.matrix(tmp)

####################################
#Necessary data-dependent variables#
####################################

#ticks     <- c(1, 31, 61, 91, 121)
dates     <- unique(y$date)
n_ts      <- ncol(y_augm)
t_len     <- nrow(y_augm)
stocks    <- colnames(y_augm) 

#####################
#Estimating variance#
#####################

#Order selection
q_vec <- 10:20
r_vec <- 10:15
order_results <- c()

for (j in 1:n_ts){
  criterion_matrix <- expand.grid(q = q_vec, r = r_vec)
  
  criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
  
  for (i in 1:nrow(criterion_matrix)){
    FPE <- c()
    AIC <- c()
    AICC <- c()
    SIC <- c()
    HQ <- c()
    
    different_orders <- (1:9)
    
    for (order in different_orders){
      AR.struc      <- estimate_lrv(data = y_augm[, j], q = criterion_matrix$q[[i]],
                                    r_bar = criterion_matrix$r[[i]], p = order)
      sigma_eta_hat <- sqrt(AR.struc$vareta)
      FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
      AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
      AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order +2)/t_len))
      SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
      HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
    }
    criterion_matrix$FPE[[i]]  <- which.min(FPE)
    criterion_matrix$AIC[[i]]  <- which.min(AIC)
    criterion_matrix$AICC[[i]] <- which.min(AICC)
    criterion_matrix$SIC[[i]]  <- which.min(SIC)
    criterion_matrix$HQ[[i]]   <- which.min(HQ)
  }
  maxim <- min(criterion_matrix[, 3:7])
  order_results <- c(order_results, maxim)
  # cat("For the stock ", colnames(y_augm)[j],
  #     "for the volatility the results are as follows: ", max(criterion_matrix$FPE), " ",
  #     max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ",
  #     max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
}

#Calculating each sigma_i separately
sigmahat_vector <- c()
ahat_vector <- c()
for (i in 1:n_ts){
  AR.struc        <- estimate_lrv(data = y_augm[, i], q = q, r_bar = r,
                                  p = order_results[i])  
                                  #p = 1)
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
  ahat_vector <- c(ahat_vector, AR.struc$ahat)
}


#########
#Testing#
#########

#Constructing the grid
u_grid      <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
h_grid      <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
h_grid      <- h_grid[h_grid > log(t_len) / t_len]
grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

result <- multiscale_test(data = y_augm, sigma_vec = sigmahat_vector,
                          alpha = alpha,  n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

statistic_vector <- result
statistic_value  <- max(statistic_vector)
if (statistic_value > quant) {
  statistic_matrix  <- matrix(statistic_vector, ncol = n_ts, nrow =  n_ts, byrow = FALSE)
  statistic_matrix  <- forceSymmetric(statistic_matrix, uplo = "U")
  statistic_matrix  <- as.dist(statistic_matrix)
  clustering        <- hclust(statistic_matrix, method = "complete")
  groups            <- cutree(clustering, h = quant)
  number_of_groups  <- max(groups)
} else {
  number_of_groups <- 1
  groups           <- rep(1, n_ts)
}



#####################
#Plots for the paper#
#####################

#Producing the smoothed curves using local linear estimator

produce_smoothed_plots(matrix = hp_log,
                       pdfname = "plots/smoothed_hp_data.pdf",
                       y_min = min(hp_log) + 0.02, y_max = max(hp_log) - 0.01,
                       ticks_at =  ticks, ticks_labels = dates[ticks],
                       yaxp_ = c(0, 6, 6))

produce_smoothed_plots(matrix = hp_log_augm,
                       pdfname = "plots/smoothed_hp_data_augmented.pdf",
                       y_min = min(hp_log_augm) + 0.01,
                       y_max = max(hp_log_augm) - 0.01,
                       ticks_at =  ticks, ticks_labels = dates[ticks],
                       yaxp_ = c(-3, 3, 6))


#Producing plots with the final results
for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  if (result$stat_pairwise[i, j] > result$quant){
    produce_plots_hp(results = result, data_i = hp_log_augm[, i],
                     data_j = hp_log_augm[, j], at_ = ticks,
                     labels_ = dates[ticks], name_i = countries[i],
                     name_j = countries[j])
  }
}

# ################################################################################
# #################  GROWTH RATE of the log of house prices  #####################
# ################################################################################
# 
# #####################
# #Estimating variance#
# #####################
# 
# #Order selection
# q_vec <- 10:20
# r_vec <- 10:15
# order_results_2 <- c()
# 
# for (j in 1:n_ts){
#   criterion_matrix <- expand.grid(q = q_vec, r = r_vec)
#   
#   criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
#   
#   for (i in 1:nrow(criterion_matrix)){
#     FPE <- c()
#     AIC <- c()
#     AICC <- c()
#     SIC <- c()
#     HQ <- c()
#     
#     different_orders <- (1:9)
#     
#     for (order in different_orders){
#       AR.struc      <- estimate_lrv(data = hp_growth_rate_augm[, j],
#                                     q = criterion_matrix$q[[i]],
#                                     r_bar = criterion_matrix$r[[i]], p = order)
#       sigma_eta_hat <- sqrt(AR.struc$vareta)
#       FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
#       AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
#       AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order + 2)/t_len))
#       SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
#       HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
#     }
#     criterion_matrix$FPE[[i]]  <- which.min(FPE)
#     criterion_matrix$AIC[[i]]  <- which.min(AIC)
#     criterion_matrix$AICC[[i]] <- which.min(AICC)
#     criterion_matrix$SIC[[i]]  <- which.min(SIC)
#     criterion_matrix$HQ[[i]]   <- which.min(HQ)
#   }
#   maxim <- max(criterion_matrix[, 3:7])
#   order_results_2 <- c(order_results_2, maxim)
#   cat("For the country ", colnames(hp_growth_rate_augm)[j],
#       "for the growth rate of house prices the results are as follows: ",
#       max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ",
#       max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ",
#       max(criterion_matrix$HQ), " \n")
# }
# 
# #Calculating each sigma_i separately
# sigmahat_vector_2 <- c()
# for (i in 1:n_ts){
#   AR.struc        <- estimate_lrv(data = hp_log_augm[, i], q = q, r_bar = r,
#                                   #p = order_results_2[i])  
#                                   p = 1)
#   sigma_hat_i     <- sqrt(AR.struc$lrv)
#   sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_hat_i)
# }
# 
# 
# #########
# #Testing#
# #########
# 
# #Constructing the grid
# u_grid      <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
# h_grid      <- seq(from = 5 / t_len, to = 1 / 4, by = 1 / t_len)
# h_grid      <- h_grid[h_grid > log(t_len) / t_len]
# grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
# 
# result <- multiscale_test(data = hp_growth_rate_augm,
#                           sigma_vec = sigmahat_vector_2, alpha = alpha,
#                           n_ts = n_ts, grid = grid, sim_runs = sim_runs,
#                           epidem = FALSE)
# 
# 
# #####################
# #Plots for the paper#
# #####################
# 
# #Producing the smoothed curves using local linear estimator
# 
# 
# produce_smoothed_plots(matrix = hp_growth_rate,
#                        pdfname = "plots/growth_rate/smoothed_hp_data.pdf",
#                        y_min = min(hp_growth_rate), y_max = max(hp_growth_rate),
#                        ticks_at =  ticks, ticks_labels = dates[ticks])
# 
# produce_smoothed_plots(matrix = hp_growth_rate_augm,
#                        pdfname = "plots/growth_rate/smoothed_hp_data_augmented.pdf",
#                        y_min = min(hp_growth_rate_augm),
#                        y_max = max(hp_growth_rate_augm),
#                        ticks_at =  ticks, ticks_labels = dates[ticks])
# 
# 
# #Producing plots with the final results
# for (l in seq_len(nrow(result$ijset))){
#   i <- result$ijset[l, 1]
#   j <- result$ijset[l, 2]
#   if (result$stat_pairwise[i, j] > result$quant){
#     produce_plots(results = result, data_i = hp_growth_rate_augm[, i],
#                   data_j = hp_growth_rate_augm[, j], at_ = ticks,
#                   labels_ = dates[ticks], name_i = countries[i],
#                   name_j = countries[j])
#   }
# }
# 
