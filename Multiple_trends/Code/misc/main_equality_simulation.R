library(car)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Equality/C_code/psihat_statistic.R")
dyn.load("Equality/C_code/psihat_statistic_ll.dll")
dyn.load("Equality/C_code/psihat_statistic_nw.dll")
source("Equality/C_code/estimating_sigma.R")
dyn.load("Equality/C_code/estimating_sigma.dll")
source("Equality/functions.R")
source("Equality/simulations_based_on_data.R")


##############################
#Defining necessary constants#
##############################

N_ts_sim <- 15 #number of different time series for simulation
N_rep    <- 1000 #number of simulations for calculating size and power

different_T     <- c(250, 350, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported

a_hat <- 0.267 #Parameters that are used in simulations but were estimated beforehand in apllication (shape)
sigma <- 0.59 #Variance = 0.35


############################
#Calculating size and power#
############################

results_size     <- simulations_size(a_hat, sigma, N_ts_sim, N_rep, different_T, different_alpha, kernel_method)
results_power    <- simulations_power(a_hat, sigma, N_ts_sim, N_rep, different_T, different_alpha, kernel_method)


###################################
#Implementing clustering algorithm#
###################################

simulations_clustering(a_hat, sigma, N_ts_sim, N_rep, different_T, different_alpha, kernel_method)
results_clusters <- clustering_analysis(N_ts_sim, different_T, different_alpha)