library(xtable)
library(matlib)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("Shape/functions.R")
source("Shape/C_code/estimating_sigma.R")
dyn.load("Shape/C_code/estimating_sigma.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")

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

SiZer_single_estimation <- function(i, h, y_data, gamma_estimated, q_h, X_matrix, W_vector, XtW, XtWX_inv){
  T_data <- length(y_data)
  
  result         <- calculating_estimator_and_variance(i, h, gamma_estimated, y_data, X_matrix, W_vector, XtW, XtWX_inv)
  m_hat_prime    <- result[[1]][2, 1]
  sd_m_hat_prime <- sqrt(result[[2]][2])

  return(result)
}
###############################
#Defining necessary parameters#
###############################

a_1         <- 0.5
sigma_eta   <- 1
T_size      <- 200
alpha       <- 0.05
different_i <- seq(from = 1/T_size, to = 1, by = 1/T_size)
different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
N_rep       <- 100  

i_0  <- round(T_size/2) /T_size #Auxiliary data point to calculate Effective Sample Size
L1 = floor(sqrt(T_size))
L2 = floor(2 * sqrt(T_size))
kernel_ind = 2

############################################################
#Calculating everything that does not depend on Y for SiZer#
############################################################
gamma = c()
for (k in 0:(T_size-1)){                                            #gamma(k) = sigma_eta^2 * a_1^|k| / (1 - a_1^2)
  gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))
}

#Calculating \Var(\bar{Y}) based on the true values of gamma(k)
true_var <- autocovariance_function_AR1(0, a_1, sigma_eta) / T_size
for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * autocovariance_function_AR1(k, a_1, sigma_eta)}

T_star   <- gamma[1]/true_var

SiZer_matrix          <- expand.grid(u = different_i, h = different_h) #Creating a dataframe with all possible combination of i and h
SiZer_matrix$values   <- numeric(nrow(SiZer_matrix)) # Setting the values of the statistic to be zero
SiZer_matrix$sd       <- numeric(nrow(SiZer_matrix)) # Setting the values of standard deviation to be zero
SiZer_matrix$ESS_star <- numeric(nrow(SiZer_matrix)) # Setting the values of ESS* to be zero
SiZer_matrix$q_h      <- numeric(nrow(SiZer_matrix)) # Setting the values of the gaussian quantile to be zero
SiZer_matrix$small_ESS<- numeric(nrow(SiZer_matrix)) # Later we will delete all the row such that ESS* is too small
SiZer_matrix$lambda   <- lambda(SiZer_matrix[['h']]) #Calculating the lambda(h) in order to speed up the function psistar_statistic


for (row in 1:nrow(SiZer_matrix)){
  
  i = SiZer_matrix[row, 'u']
  h = SiZer_matrix[row, 'h']
  
  ESS      <- sum(sapply((i - seq(1/T_size, 1, by = 1/T_size))/h, epanechnikov_kernel))/epanechnikov_kernel(0)
  ESS_star <- (T_star/T_size) * ESS 
  l        <- T_size / ESS_star
  q_h      <- qnorm((1 + (1 - alpha)^(1 / l))/2)
  
  if (ESS_star <= 5){
    SiZer_matrix[row, 'small_ESS'] <- 1
  } else {
    x_matrix     <- matrix(c(rep(1, T_size), seq(1/T_size - i, 1-i, by = 1/T_size)), nrow=T_size, byrow=FALSE)
    w_vector     <- sapply(seq(1/T_size - i, 1-i, by = 1/T_size)/h, FUN = epanechnikov_kernel)/h
    XtW          <- t(apply(x_matrix,2,function(x){x*w_vector}))   #t(X) %*% diag(w)   computed faster.  :)
    XtWX_inverse <- tryCatch({solve(XtW %*% x_matrix)},  
                            error = function(e) {print("Something is wrong, the matrix can not be inverted")})
  
    XtSigma <- matrix(data = NA, nrow =2, ncol =T_size)
    for (j in 1:T_size){      #t(X) %*% Sigma computed faster in order not to save the whole Sigma matrix (T_size * T_size) in the memory
      XtSigma[1, j] = 0
      XtSigma[2, j] = 0
      for (k in 1:T_size){
        XtSigma[1, j] = XtSigma[1, j] + gamma[abs(k - j) + 1] * w_vector[k] * w_vector[j]
        XtSigma[2, j] = XtSigma[2, j] + gamma[abs(k - j) + 1] * w_vector[k] * w_vector[j] * (k/T_size - i)   
      }
    }
    sd = sqrt((XtWX_inverse %*% XtSigma %*% x_matrix %*% XtWX_inverse)[2,2])
  
    SiZer_matrix[row, 'sd']       = sd
    SiZer_matrix[row, 'q_h']      = q_h
    SiZer_matrix[row, 'ESS_star'] = ESS_star
  }
}

SiZer_matrix <- SiZer_matrix[SiZer_matrix$small_ESS == 0,]

#################################################################
#Calculating everything that does not depend on Y for our method#
#################################################################


#Calculating Gaussian statistic for our own method
gaussian_statistic_distribution <- replicate(1000, {
  z = rnorm(T_size, 0, 1)
  psistar_statistic_ll(z, SiZer_matrix, kernel_ind, 1)
})
gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)





size_of_the_test <- replicate(N_rep, {
  y_data_ar_1 <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
  sigmahat    <- estimating_sigma_for_AR1(y_data_ar_1, L1, L2)[[1]]
  g_t_set     <- psihat_statistic_ll(y_data_ar_1, SiZer_matrix, kernel_ind, sigmahat)[[1]]

  results_our   <- c()
  results_their <- c()
  
  for (row in 1:nrow(g_t_set)){
    i              = g_t_set[row, 'u']
    h              = g_t_set[row, 'h']
    q_h            = g_t_set[row, 'q_h']
    sd_m_hat_prime = g_t_set[row, 'sd']

    x_matrix     <- matrix(c(rep(1, T_size), seq(1/T_size - i, 1-i, by = 1/T_size)), nrow=T_size, byrow=FALSE)
    w_vector     <- sapply(seq(1/T_size - i, 1-i, by = 1/T_size)/h, FUN = epanechnikov_kernel)/h
    XtW          <- t(apply(x_matrix,2,function(x){x*w_vector}))   #t(X) %*% diag(w)   computed faster.  :)
    XtWX_inverse <- tryCatch({solve(XtW %*% x_matrix)},  
                           error = function(e) {print("Something is wrong, the matrix can not be inverted")})
  
    m_hat_prime <- (XtWX_inverse %*% XtW %*% y_data_ar_1)[2]
  
    if (m_hat_prime - q_h * sd_m_hat_prime > 0){
      results_their = c(results_their, 1)
    } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
      results_their = c(results_their, -1)
    } else {
      results_their = c(results_their, 0)
    }
    if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
      results_our = c(results_our, 1)
    } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
      results_our = c(results_our, -1)
    } else {
      results_our = c(results_our, 0)
    }
    
    #cat("i = ", i, " and h = ", h, ": interval = [", m_hat_prime - q_h * sd_m_hat_prime, ",",m_hat_prime + q_h * sd_m_hat_prime, "]\n")
  }
  c(sum(abs(results_their)), sum(abs(results_our)))
})

power_of_the_test <- replicate(N_rep, {
  y_data_ar_1 <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
  m = numeric(T_size)
  for (i in 1:T_size) {m[i] = i - 0.5*T_size}
  y_data_ar_1_with_trend <- y_data_ar_1 + m
  sigmahat    <- estimating_sigma_for_AR1(y_data_ar_1_with_trend, L1, L2)[[1]]
  g_t_set     <- psihat_statistic_ll(y_data_ar_1_with_trend, SiZer_matrix, kernel_ind, sigmahat)[[1]]
  
  results_our   <- c()
  results_their <- c()
  
  for (row in 1:nrow(g_t_set)){
    i              = g_t_set[row, 'u']
    h              = g_t_set[row, 'h']
    q_h            = g_t_set[row, 'q_h']
    sd_m_hat_prime = g_t_set[row, 'sd']
    
    x_matrix     <- matrix(c(rep(1, T_size), seq(1/T_size - i, 1-i, by = 1/T_size)), nrow=T_size, byrow=FALSE)
    w_vector     <- sapply(seq(1/T_size - i, 1-i, by = 1/T_size)/h, FUN = epanechnikov_kernel)/h
    XtW          <- t(apply(x_matrix,2,function(x){x*w_vector}))   #t(X) %*% diag(w)   computed faster.  :)
    XtWX_inverse <- tryCatch({solve(XtW %*% x_matrix)},  
                             error = function(e) {print("Something is wrong, the matrix can not be inverted")})
    
    m_hat_prime <- (XtWX_inverse %*% XtW %*% y_data_ar_1_with_trend)[2]
    
    if (m_hat_prime - q_h * sd_m_hat_prime > 0){
      results_their = c(results_their, 1)
    } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
      results_their = c(results_their, -1)
    } else {
      results_their = c(results_their, 0)
    }
    if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
      results_our = c(results_our, 1)
    } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
      results_our = c(results_our, -1)
    } else {
      results_our = c(results_our, 0)
    }
    
    #cat("i = ", i, " and h = ", h, ": interval = [", m_hat_prime - q_h * sd_m_hat_prime, ",",m_hat_prime + q_h * sd_m_hat_prime, "]\n")
  }
  c(sum(abs(results_their)), sum(abs(results_our)))
})