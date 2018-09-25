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
different_a      <- c(-0.95, -0.5, -0.25, 0.25, 0.5, 0.95)
different_slopes <- c(0, 1, 4, 10)

#true_sigma <- sqrt(sigma_eta^2/((1 - a_1)^2))
#a_1    <- -0.95
T_size <- 500
#slope  <- 4

L1 <- 20
L2 <- 30
K1 <- p + 1
K2 <- 10

#######################################################
#Calculating histograms for extreme negative value a_1#
#######################################################

for (a_1 in different_a){
  for (slope in different_slopes){
    histograms_for_variance_estimators(a_1, sigma_eta, T_size, p, slope, N_rep, K1, K2, L1, L2)
  }
}
