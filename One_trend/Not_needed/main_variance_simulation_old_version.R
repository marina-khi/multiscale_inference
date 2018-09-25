library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/estimating_sigma_new.R")

source("Shape/C_code/estimating_sigma.R")
dyn.load("Shape/C_code/estimating_sigma.dll")


###############################
#Defining necessary parameters#
###############################

N_rep            <- 100 #Number of replications for comparison of the estimates
different_T      <- c(500) #Different lengths of time series for which we compare the estimates
different_a      <- c(-0.5, -0.25, 0.25, 0.5)
different_slopes <- c(0)
sigma_eta        <- 1



for (a_1 in different_a){
  true_sigma <- sqrt(sigma_eta^2/((1 - a_1)^2))
  for (T_size in different_T){
    line_trend  <- numeric(T_size)
    L1      <- floor(sqrt(T_size))
    L2      <- floor(2 * sqrt(T_size))
    # for (slope in different_slopes){
    slope = 4
    for (i in 1:T_size) {line_trend[i] = (i - 0.5*T_size) * slope/T_size}
    a_hat_hall_vect           <- c() 
    a_hat_max_likelihood_vect <- c()
    a_hat_old_vect            <- c() 
    a_hat_overidentified_vect <- c()
    a_hat_oracle_vect         <- c()
    for (i in 1:N_rep){
        y_data_ar_1_with_trend <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta)) + line_trend
      a_hat_hall           <- estimating_sigma_for_AR1(y_data_ar_1_with_trend, L1, L2)[[2]]
      a_hat_max_likelihood <- estimating_sigma_max_likelihood(y_data_ar_1_with_trend, 1)
      a_hat_old            <- estimating_sigma_old_simple(y_data_ar_1_with_trend, 1)
      a_hat_overidentified <- estimating_sigma_overidentified(y_data_ar_1_with_trend, 1)
      a_hat_oracle         <- estimating_sigma_oracle(y_data_ar_1_with_trend, 1)
      
#      cat("Difference with old method: ", var_old, ", a_1 = ", a_1, ", slope = ", slope, ", T = ", T_size, "\n")
#      cat("Estimate of parameter a_1 with overidentified method: ", a_hat_overidentified, ", a_1 = ", a_1, ", slope = ", slope, ", T = ", T_size, "\n")
#      cat("Estimate of parameter a_1 with twice overidentified method: ", a_hat_overidentified_twice, ", a_1 = ", a_1, ", slope = ", slope, ", T = ", T_size, "\n")
      
#      cat("Estimate of parameter a_1 with oracle method: ", a_hat_oracle, ", a_1 = ", a_1, ", slope = ", slope, ", T = ", T_size, "\n")
#      cat("\n")
      a_hat_hall_vect           <- c(a_hat_hall_vect, a_hat_hall) 
      a_hat_max_likelihood_vect <- c(a_hat_max_likelihood_vect, a_hat_max_likelihood)
      a_hat_old_vect            <- c(a_hat_old_vect, a_hat_old) 
      a_hat_overidentified_vect <- c(a_hat_overidentified_vect, a_hat_overidentified)
      a_hat_oracle_vect         <- c(a_hat_oracle_vect, a_hat_oracle) 
      
#      vect <- c(vect, T_size, a_1, a_hat_old, a_hat_overidentified, a_hat_overidentified_twice, a_hat_overidentified_3, a_hat_oracle)
    }
    smallest_value <- min(a_hat_hall_vect, a_hat_max_likelihood_vect, a_hat_old_vect, a_hat_overidentified_vect, a_hat_oracle_vect)
    biggest_value <- max(a_hat_hall_vect, a_hat_max_likelihood_vect, a_hat_old_vect, a_hat_overidentified_vect, a_hat_oracle_vect)
    
    hist1 <- hist(a_hat_hall_vect, main=paste0("Hall and Van Keilegom, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), plot = FALSE)
    hist2 <- hist(a_hat_max_likelihood_vect, main = paste0("Maximul likelihood, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), plot = FALSE)
    hist3 <- hist(a_hat_old_vect, main = paste0("Invertible matrix?, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), plot = FALSE)
    hist4 <- hist(a_hat_overidentified_vect, main = paste0("Overidentified, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), plot = FALSE)
    hist5 <- hist(a_hat_oracle_vect, main = paste0("Oracle method, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), plot = FALSE)
    
    highestCount <- max(hist1$counts, hist2$counts, hist3$counts, hist4$counts, hist5$counts)
    
    
    hist(a_hat_hall_vect, main=paste0("Hall and Van Keilegom, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), ylim=c(0,highestCount))
    #hist(a_hat_max_likelihood_vect, main = paste0("Maximul likelihood, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), ylim=c(0,highestCount))
    #hist(a_hat_old_vect, main = paste0("Invertible matrix?, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), ylim=c(0,highestCount))
    hist(a_hat_overidentified_vect, main = paste0("Overidentified, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), ylim=c(0,highestCount))
    #hist(a_hat_oracle_vect, main = paste0("Oracle method, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.05), ylim=c(0,highestCount))
  }
}

# matrix_ <- matrix(vect, nrow = length(different_T) * length(different_a), ncol = 7, byrow = TRUE)
# colnames(matrix_) <- c("T", "a_1", "simple", "overidentified", "overidentified*2", "overidentified*3", "oracle")
# 
# print.xtable(xtable(matrix_, digits = c(3), align = "cccccccc"), type="latex",  file="Paper/Plots/variance_estimators.tex", include.colnames = TRUE, include.rownames =FALSE )
