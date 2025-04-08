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
sigma <- 0.25

#For the fixed effects
rho      <- 0.25 #covariance between the fixed effects
n_rep    <- 5000 #number of simulations for calculating size and power
sim_runs <- 5000 #number of simulations to calculate the Gaussian quantiles

#Different parameters
different_T     <- c(100, 250, 500) #Different lengths of time series  
different_b     <- c(0.25, 0.5, 0.75) #Zero is for calculating the size

#Parameters for the estimation of long-run-variance
q <- 25 
r <- 10

#For parallel computation
numCores  <- round(parallel::detectCores() * .80)

####################################
#Calculating the coefficient values#
####################################


for (t_len in different_T){
  set.seed(seed)
  k <- match(t_len, different_T)

  #Calculating the coefficients
  tic()
  cl <- makePSOCKcluster(numCores)
  registerDoParallel(cl)
  foreach (val = 1:n_rep, .combine = "cbind") %dopar% {
    source("functions/functions.R")
    repl_precision(rep_ = val, n_ts_ = n_ts, t_len_ = t_len,
                  a_ = a, sigma_ = sigma,
         beta_ = beta, a_x_vec_ = a_x_vec, phi_ = phi, rho_ = rho,
         different_b_ = different_b,
         q_ = q, r_ = r)
    # Loop one-by-one using foreach
  } -> simulated_beta_and_sigma
  stopCluster(cl)
  toc()
  
  for (j in 1:length(different_b)){
    simulated_beta <- simulated_beta_and_sigma[((j - 1) * (length(different_b) + 1) + 1):(j * (length(different_b) + 1) - 1), ]
    simulated_sigma <- simulated_beta_and_sigma[j * (length(different_b) + 1), ]
    

    filename = paste0("output/revision/sigma_histogram_T_", t_len, "_b_", different_b[j] * 100, ".pdf")
    pdf(filename, width = 7, height = 5, paper="special")
#    layout(matrix(c(1, 2), ncol=1), widths=c(2.4, 2.4),
#           heights=c(1.5, 1.8), TRUE)
    
    #Setting the layout of the graphs
    par(cex = 1, tck = -0.025)
    par(mar = c(2.5, 4.2, 2, 0)) #Margins for each plot
#    par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
    
    hist(simulated_sigma, main = "")
    abline(v = 1/3, col = "red")
    title(main = "Histiogram of the estimated square root of the long-run variance", font.main = 1, line = 0.5)
    dev.off()
    
    filename = paste0("output/revision/beta_histogram_T_", t_len, "_b_", different_b[j] * 100, ".pdf")
    pdf(filename, width = 21, height = 5, paper="special")
    layout(matrix(c(1, 2, 3), ncol = 3), widths=c(6.6, 6.6, 6.6),
           heights=c(5, 5, 5), TRUE)
    
    #Setting the layout of the graphs
    par(cex = 1, tck = -0.025)
    par(mar = c(2.5, 0.5, 2, 0)) #Margins for each plot
    par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
    
    hist(simulated_beta[1, ], main = "", breaks = seq(from = 0.95, to = 1.05, by= 0.01))
    abline(v = 1, col = "red")
    title(main = "Histiogram of the estimated first coefficient", font.main = 1, line = 0.5)

    hist(simulated_beta[2, ], main = "", breaks = seq(from = 0.95, to = 1.05, by= 0.01))
    abline(v = 1, col = "red")
    title(main = "Histiogram of the estimated second coefficient", font.main = 1, line = 0.5)
    
    hist(simulated_beta[3, ], main = "", breaks = seq(from = 0.95, to = 1.05, by= 0.01))
    abline(v = 1, col = "red")
    title(main = "Histiogram of the estimated third coefficient", font.main = 1, line = 0.5)
    
    dev.off()
    
  }
}
