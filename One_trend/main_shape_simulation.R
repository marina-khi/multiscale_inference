library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/estimating_sigma_new.R")

source("Shape/C_code/estimating_sigma.R")
dyn.load("Shape/C_code/estimating_sigma.dll")
dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/simulations_based_on_data.R")


###############################
#Defining necessary parameters#
###############################

N               <- 1000 #Number of replications for calculating the size and the power of the test
different_T     <- c(250, 350, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported
test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 

###################################################################################
#Estimating AR(1) parameters from the data in order to use them in the simulations#
###################################################################################

temperature             <- read.table("Shape/data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr            <- temperature[temperature$YEAR > -99, 'YEAR']

#Tuning parameters
T_tempr <- length(yearly_tempr)
L1      <- floor(sqrt(T_tempr))
L2      <- floor(2 * sqrt(T_tempr))

result    <- estimating_sigma_for_AR1(yearly_tempr, L1, L2)
a_hat     <- result[[2]] #Estimation of the AR coefficient
sigma_eta <- result[[3]] #Estimation of the sqrt of the variance of the innovation 


############################################
#Calculating the power and size of the test#
############################################
simulations_based_on_data(N, different_T, different_alpha, a_hat, sigma_eta, test_problem, kernel_method)