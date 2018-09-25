library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/estimating_sigma_new.R")

dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/simulations_based_on_data.R")


###############################
#Defining necessary parameters#
###############################

N_rep           <- 1000 #Number of replications for calculating the size and the power of the test
different_T     <- c(250, 350, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power
different_a1    <- c(-0.5, -0.25, 0.25, 0.5)
sigma_eta       <- 1

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported
test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 

L1 <- 20
L2 <- 30  

######################################################
#Calculating the power and size of the test for AR(1)#
######################################################

# for (a_hat in different_a1){
#   PDFpartialPath = paste0("Paper/Plots/finite_sample_properties_a1_", a_hat*100)
#   simulations_general(N_rep, different_T, different_alpha, a_hat, sigma_eta, order = 1, test_problem, kernel_method, filename= PDFpartialPath, L1, L2, K1 = 1 + 1, K2 = 10)
# }

################################################################################
#Calculating the power and size of the test for one specification based on data#
################################################################################


temperature             <- read.table("Shape/data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr            <- temperature[temperature$YEAR > -99, 'YEAR']

simulations_based_on_data(N_rep, yearly_tempr, different_alpha, different_slopes = c(0, 1.5, 1.875, 2.5), order = 2, test_problem, L1, L2, K1 = 2 + 1, K2 = 10)