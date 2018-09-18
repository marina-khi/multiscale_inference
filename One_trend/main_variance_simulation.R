#library(xtable)
#options(xtable.floating = FALSE)
#options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/estimating_sigma_new.R")

source("Shape/C_code/estimating_sigma.R")
dyn.load("Shape/C_code/estimating_sigma.dll")


###############################
#Defining necessary parameters#
###############################

N                <- 1000 #Number of replications for comparison of the estimates
different_T      <- c(250, 350, 500, 1000) #Different lengths of time series for which we compare the estimates
different_a      <- c(-0.5, -0.25, 0.25, 0.5)
different_slopes <- c(0)
sigma_eta        <- 1

for (a_1 in different_a){
  true_sigma <- sqrt(sigma_eta^2/((1 - a_1)^2))
  for (T_size in different_T){
    line_trend  <- numeric(T_size)
#    L1      <- floor(sqrt(T_size))
#    L2      <- floor(2 * sqrt(T_size))
    for (slope in different_slopes){
      for (i in 1:T_size) {line_trend[i] = (i - 0.5*T_size) * slope/T_size}
      y_data_ar_1_with_trend <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta)) + line_trend
#      var_old <- estimating_sigma_for_AR1(y_data_ar_1_with_trend, L1, L2)[[3]]
      a_hat_overidentified <- estimating_sigma_overidentified(y_data_ar_1_with_trend, 1)
      a_hat_oracle         <- estimating_sigma_oracle(y_data_ar_1_with_trend, 1)
      
#      cat("Difference with old method: ", var_old, ", a_1 = ", a_1, ", slope = ", slope, ", T = ", T_size, "\n")
      cat("Estimate of parameter a_1 with overidentified method: ", a_hat_overidentified, ", a_1 = ", a_1, ", slope = ", slope, ", T = ", T_size, "\n")
      cat("Estimate of parameter a_1 with oracle method: ", a_hat_oracle, ", a_1 = ", a_1, ", slope = ", slope, ", T = ", T_size, "\n")
      cat("\n")
    }
  }
}

