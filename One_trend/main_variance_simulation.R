library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/simulations_variance.R")

source("Shape/C_code/estimating_sigma.R")
dyn.load("Shape/C_code/estimating_sigma.dll")


###############################
#Defining necessary parameters#
###############################

p                <- 1 #Order of AR(p)
N_rep            <- 1000 #Number of replications for comparison of the estimates
sigma_eta        <- 1 # Sqrt root of the variance of the innovation \eta_t
#different_T      <- c(100) #Different lengths of time series for which we compare the estimates
#different_a      <- c(-0.95)
#different_slopes <- c(4)

#true_sigma <- sqrt(sigma_eta^2/((1 - a_1)^2))
a_1    <- -0.95
T_size <- 500
slope <- 4

#######################################################
#Calculating histograms for extreme negative value a_1#
#######################################################


histograms_for_variance_estimators(a_1, sigma_eta, T_size, p, slope, N_rep, K1, K2, L1, L2)