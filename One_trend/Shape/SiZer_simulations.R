library(xtable)
library(matlib)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")


source("Shape/C_code/estimating_sigma.R")
dyn.load("Shape/C_code/estimating_sigma.dll")
dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/simulations_based_on_data.R")


###############################
#Defining necessary parameters#
###############################

a_1 <- 0.9
sigma_eta <- 1
T_size <- 500


####################################
#Estimating autocovariance function#
####################################

gamma = c()
for (i in 0:(T_size-1)){gamma = c(gamma, autocovariance_function_AR1(i, a_1, sigma_eta))}

i = 0.5
h = 5 / T_size
y_data_ar_1  <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
result       <- calculating_estimator_and_variance(T_size, i, h, gamma, y_data_ar_1)
variance_hat <- result[[2]]
estimator    <- result[[1]]
m_hat_prime  <- estimator[2,1]
sd           <- variance_hat[2, 2]


