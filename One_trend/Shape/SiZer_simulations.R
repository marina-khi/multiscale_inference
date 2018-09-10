library(xtable)
library(matlib)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")

############################################
#These are all the SiZer-specific functions#
############################################

#Implementation of \widehat{Var}(\bar{Y}) from "SiZer for time series" paper
estimating_variance_of_data <- function(data){
  T_size = length(data)
  divided_sample = split(data, cut(seq_along(data), sqrt(T_size), labels = FALSE)) 
  divided_sample_mean = sapply(divided_sample, FUN = mean)
  result = sum((divided_sample_mean - mean(divided_sample_mean))^2)/ ((length(divided_sample_mean) - 1) * length(divided_sample_mean))
  return(result)
}

#Implementation of calculating q(h) from "SiZer for time series" paper
gaussian_quantile_for_sizer <- function(i_0, h, y_data, var_epsilon, var_y_bar, ESS, alpha){
  T_data   <- length(y_data)
  T_star   <- var_epsilon/var_y_bar
  ESS_star <- (T_star/T_data) * ESS 
  l        <- T_data / ESS
  return(list(qnorm((1 + (1 - alpha)^(1 / l))/2), ESS_star))
}



#Calculate autocovariance function for AR(1) model \varepsilon_t = a_1 \varepsilon_{t-1} + \eta_t based on the coefficients of the model
autocovariance_function_AR1 <- function(k, a_1, sigma_eta){
  if (k%%1==0)
  {
    result = sigma_eta * sigma_eta * a_1^(abs(k)) / (1 - a_1 * a_1)
  } else {
    print('Check the input: k is not integer')
  }
  return(result)
}

calculating_estimator_and_variance <- function(T_size, i_0, h, gamma, data){
  sigmahat_matrix <- matrix(data = NA, nrow =T_size, ncol =T_size)
  w_vector = c()
  for (i in 1:T_size){
    for (j in 1:T_size){
      sigmahat_matrix[i,j] = gamma[abs(i - j) + 1] * epanechnikov_kernel((i/T_size - i_0)/h) * epanechnikov_kernel((j/T_size - i_0)/h)/(h^2)
    }
    w_vector = c(w_vector, epanechnikov_kernel((i/T_size - i_0)/h) / h)
  }
  w_matrix = diag(w_vector)
  x_matrix =  matrix(c(rep(1, T_size), seq(1/T_size - i_0, 1-i_0, by = 1/T_size)), nrow=T_size, byrow=FALSE)
  inverse_matrix = tryCatch({inv(t(x_matrix) %*% w_matrix %*% x_matrix)},
                            error = function(e) {print("Something is wrong, the matrix can not be inverted")})
  variance_matrix = inverse_matrix %*% (t(x_matrix) %*% sigmahat_matrix %*% x_matrix) %*% inverse_matrix
  estimator = inverse_matrix %*% (t(x_matrix) %*% w_matrix %*% data)
  return(list(estimator, variance_matrix))
}

SiZer_single_estimation <- function(i, h, y_data, gamma_estimated, q_h){
  T_data <- length(y_data)
  
  result         <- calculating_estimator_and_variance(T_data, i, h, gamma_estimated, y_data)
  m_hat_prime    <- result[[1]][2, 1]
  sd_m_hat_prime <- sqrt(result[[2]][2, 2])

  if (m_hat_prime - q_h * sd_m_hat_prime > 0){
    result = 1
  } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
    result = -1
  } else {
    result = 0
  }
  #cat("i = ", i, " and h = ", h, ": interval = [", m_hat_prime - q_h * sd_m_hat_prime, ",",m_hat_prime + q_h * sd_m_hat_prime, "]\n")
  return(result)
}
###############################
#Defining necessary parameters#
###############################

a_1       <- 0.9
sigma_eta <- 1
T_size    <- 500
alpha     <- 0.05

####################################
#Estimating autocovariance function#
####################################

gamma = c()
for (k in 0:(T_size-1)){gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))}

i_0  <- round(T_size/2) /T_size #Auxiliary data point to calculate Effective Sample Size
h    <- 0.1

y_data_ar_1       <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
variance_of_y_bar <- estimating_variance_of_data(y_data_ar_1)
ESS               <- sum(sapply((i_0 - seq(0, 1, by = 1/T_size))/h, epanechnikov_kernel))/epanechnikov_kernel(0)

g_t_set        <- creating_g_set(T_size, "nw") #Here we use NW method to create the set of locations and bandwidth because we want the restriction [u-h, u+h] \in [0,1] 
g_t_set$result <- numeric(nrow(g_t_set)) # Setting the results of SiZer method to be zero for each location and bandwidth

true_var <- autocovariance_function_AR1(0, a_1, sigma_eta) / T_size
for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * autocovariance_function_AR1(k, a_1, sigma_eta)}


system.time({
  result_2 <- gaussian_quantile_for_sizer(i_0, h, y_data_ar_1, gamma[1], variance_of_y_bar, ESS, alpha) #Here we estimate sigma^2 (variance of epsilon from SiZer paper) by \widehat{gamma}(0)  
  q_h      <- result_2[[1]]
  ESS_star <- result_2[[2]]
  if (ESS_star >= 5){
    for (i in 0:T_size) {
      SiZer_single_estimation(i/T_size, h, y_data_ar_1, gamma, q_h)
    }
  }
})

#g_t_set$result <- with(g_t_set, SiZer_single_estimation(u, i_0, h, y_data_ar_1, gamma, variance_of_y_bar, ESS, alpha))



