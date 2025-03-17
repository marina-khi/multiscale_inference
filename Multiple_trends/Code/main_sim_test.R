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

source("functions/size_and_power.R")

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
sigma <- 0.25

#For the fixed effects
rho      <- 0.25 #covariance between the fixed effects
n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles

#Different parameters
different_T     <- c(100, 250, 500) #Different lengths of time series  
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
#different_b     <- c(0, 0.25, 0.5, 0.75) #Zero is for calculating the size
different_b     <- c(0) #Zero is for calculating the size


#Parameters for the estimation of long-run-variance
q <- 25 
r <- 10

#For parallel computation
numCores  <- round(parallel::detectCores() * .80)

# 
# ############################
# #Plotting the bump function#
# ############################
# t_len <- 250
# pdf(paste0("output/revision/bump_function.pdf"),
#     width = 12, height = 8, paper="special")
# par(mfrow = c(2, 2))
# par(mar = c(4, 3, 0.5, 0)) #Margins for each plot
# par(oma = c(0.5, 0.5, 0.5, 0.2)) #Outer margins
# 
# for (b in different_b){
#   errors <- arima.sim(model = list(ar = a),
#                       innov = rnorm(t_len, 0, sigma),
#                       n = t_len)
#   plot(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
#        y = (bump((1:t_len)/t_len) * b), ylim = c(-1.6, 1.6),
#        xlab = "", ylab = "", main = NULL,
#        type = 'l', cex = 0.8)
#   lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
#         y = (bump((1:t_len)/t_len) * b) + errors, type = "l",
#         col = "red")
#   mtext(side = 1, text = paste0("b = ", b), line = 2.3, cex = 1)
# }
# dev.off()
# 

##################################################
#Calculating the size and power for a normal grid#
##################################################

size_and_power_calculations(grid_type_ = "normal", type_of_m_ = "",
                            seed_ = seed, n_ts_ = n_ts,
                            beta_ = beta, a_x_vec_ = a_x_vec, phi_ = phi,
                            a_ = a, sigma_ = sigma,
                            rho_ = rho,
                            n_rep_ = n_rep, sim_runs_ = sim_runs,
                            different_T_ = different_T,
                            different_alpha_ = different_alpha,
                            different_b_ = different_b,
                            q_ = q, r_ = r, numCores_ = numCores,
                            filename_ext_ = "")

####################################################################################
#Calculating the size and power for a normal grid and bump functions under the null#
####################################################################################

seed <- 111222333
size_and_power_calculations(grid_type_ = "normal", type_of_m_ = "bump",
                            seed_ = seed, n_ts_ = n_ts,
                            beta_ = beta, a_x_vec_ = a_x_vec, phi_ = phi,
                            a_ = a, sigma_ = sigma,
                            rho_ = rho,
                            n_rep_ = n_rep, sim_runs_ = sim_runs,
                            different_T_ = different_T,
                            different_alpha_ = different_alpha,
                            different_b_ = different_b,
                            q_ = q, r_ = r, numCores_ = numCores,
                            filename_ext_ = "_bump_null")

##################################################
#Once more for dense grid                        #
##################################################

n_rep    <- 1000 #number of simulations for calculating size and power
sim_runs <- 1000 #number of simulations to calculate the Gaussian quantiles

size_and_power_calculations(grid_type_ = "dense", type_of_m_ = "",
                            seed_ = seed, n_ts_ = n_ts,
                            beta_ = beta, a_x_vec_ = a_x_vec, phi_ = phi,
                            a_ = a, sigma_ = sigma,
                            rho_ = rho,
                            n_rep_ = n_rep, sim_runs_ = sim_runs,
                            different_T_ = different_T,
                            different_alpha_ = different_alpha,
                            different_b_ = different_b,
                            q_ = q, r_ = r, numCores_ = numCores,
                            filename_ext_ = "_dense_grid")


##################################################
#Once more for dyadic grid                       #
##################################################

n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles

size_and_power_calculations(grid_type_ = "dyadic", type_of_m_ = "",
                            seed_ = seed, n_ts_ = n_ts,
                            beta_ = beta, a_x_vec_ = a_x_vec, phi_ = phi,
                            a_ = a, sigma_ = sigma,
                            rho_ = rho,
                            n_rep_ = n_rep, sim_runs_ = sim_runs,
                            different_T_ = different_T,
                            different_alpha_ = different_alpha,
                            different_b_ = different_b,
                            q_ = q, r_ = r, numCores_ = numCores,
                            filename_ext_ = "_dyadic_grid")
