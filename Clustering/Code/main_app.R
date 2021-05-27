#####################################################
#Analysis of covid data (provided by John Hopkins U)#
#####################################################
rm(list=ls())

library(tidyr)
library(aweek)
library(dendextend)
library(Rcpp)

require(rworldmap)

source("functions.R")
Rcpp::sourceCpp("integration.cpp")

#Defining necessary constants
b_bar  <- 2
b_grid <- seq(1, b_bar, by = 0.05)

#Half of the bandwidth window (in absolute terms)
bw_abs <- 7

#Loading the data
data  <- data_load()
t_len <- data$t_len
n_ts  <- data$n_ts
covid_mat <- data$covid_mat
countries <- names(data$covid_list)

#Grid for b and for smoothing
grid_points <- seq(1/t_len, 1, by = 1/t_len)

# #And finally calculating the distance matrix
# Delta_hat_tmp <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
# b_res         <- matrix(data = rep(NA, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
# 
# for (b in b_grid){
#   norm_b <- c()
#   norm   <- c()
#   for (k in 1:n_ts){
#     norm_b <- c(norm_b, integrate1_cpp(b = b, data_points = covid_mat[, k],
#                                        grid_points = grid_points,
#                                        bw = bw_abs/t_len, subdiv = 2000)$res)
#     norm <- c(norm, integrate1_cpp(b = 1.0, data_points = covid_mat[, k],
#                                    grid_points = grid_points,
#                                    bw = bw_abs/t_len, subdiv = 2000)$res)
#   }
#   for (i in 1:(n_ts - 1)){
#     for (j in (i + 1):n_ts){
#       delta_ij <- 1/b * integrate2_cpp(b = b, data_points_1 = covid_mat[, i],
#                                        data_points_2 = covid_mat[, j],
#                                        norm_1 = norm_b[i], norm_2 = norm[j],
#                                        grid_points = grid_points, bw = bw_abs/t_len,
#                                        subdiv=2000)$res
#       delta_ji <- 1/b * integrate2_cpp(b = b, data_points_1 = covid_mat[, j],
#                                        data_points_2 = covid_mat[, i],
#                                        norm_1 = norm_b[j], norm_2 = norm[i],
#                                        grid_points = grid_points, bw = bw_abs/t_len,
#                                        subdiv=2000)$res
#       if (b == 1) {
#         Delta_hat_tmp[i, j] <- delta_ij
#         Delta_hat_tmp[j, i] <- delta_ji
#         b_res[i, j] <- 1
#         b_res[j, i] <- 1
#       } else {
#         if (delta_ij < Delta_hat_tmp[i, j]) {
#           Delta_hat_tmp[i, j] <- delta_ij
#           b_res[i, j] <- b
#           b_res[j, i] <- 1
#         } 
#         if (delta_ji < Delta_hat_tmp[j, i]) {
#           Delta_hat_tmp[j, i] <- delta_ji
#           b_res[j, i] <- b
#           b_res[i, j] <- 1          
#         }
#       }
#     }
#   }
#   cat("b = ", b, " - success\n")
# }
# 
# #Delta_hat_tmp was a temporary non-symmetrical matrix,
# #for the distance matrix we need a symmetrical one
# Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
# for (i in 1:(n_ts - 1)){
#   for (j in (i + 1):n_ts){
#     Delta_hat[i, j] <- min(Delta_hat_tmp[i, j], Delta_hat_tmp[j, i])
#     Delta_hat[j, i] <- Delta_hat[i, j]
#   }
# }
# 
# colnames(Delta_hat) <- countries
# rownames(Delta_hat) <- countries
# colnames(b_res) <- countries
# rownames(b_res) <- countries
# 
# save(Delta_hat, b_res, file = "results_14days.RData")
load("results_14days.RData")

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)

n_cl       <- 12
#results_output(res, covid_mat, Delta_hat, b_res, n_cl, countries,
#               path = "plots/", bw_abs, grid_points, t_len)

