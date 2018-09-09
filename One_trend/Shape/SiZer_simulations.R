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
y_data_ar_1    <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
g_t_set        <- creating_g_set(T_size, "nw") #Here we use NW method to create the set of locations and bandwidth because we want the restriction [u-h, u+h] \in [0,1] 
g_t_set$result <-numeric(nrow(g_t_set)) # Setting the results of SiZer method to be zero for each location and bandwidth

g_t_set$result <- with(g_t_set, SiZer_single_estimation(u, h, y_data_ar_1, gamma))

SiZer_single_estimation <- function(i, h, y_data, gamma_estimated){
  T_data <- length(y_data)
  result         <- calculating_estimator_and_variance(T_data, i, h, gamma_estimated, y_data)
  m_hat_prime    <- result[[1]][2, 1]
  sd_m_hat_prime <- result[[2]][2, 2]
  
}

