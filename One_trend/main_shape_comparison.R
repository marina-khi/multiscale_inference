#library(xtable)
#options(xtable.floating = FALSE)
#options(xtable.timestamp = "")
source("Shape/functions.R")

#source("Shape/C_code/estimating_sigma.R")
#dyn.load("Shape/C_code/estimating_sigma.dll")
dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/SiZer_simulations.R")


###############################
#Defining necessary parameters#
###############################

N_rep           <- 1000 #Number of replications for calculating the size and the power of the test
different_T     <- c(500) #Different lengths of time series for which we calculate size and power
#different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

a_1         <- 0.5
sigma_eta   <- 1
alpha       <- 0.05

for (T_size in different_T){
  PDFPath = paste0("Paper/Plots/SiZer_comparison_", T_size, ".pdf")
  SiZer_simulations(T_size, a_1, sigma_eta, alpha, N_rep, PDFpath)
}




