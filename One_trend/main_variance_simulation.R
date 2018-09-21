library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
#source("Shape/estimating_sigma_new.R")

source("Shape/C_code/estimating_sigma.R")
dyn.load("Shape/C_code/estimating_sigma.dll")


###############################
#Defining necessary parameters#
###############################

p                <- 1 #Order of AR(p)
N_rep            <- 100 #Number of replications for comparison of the estimates
different_T      <- c(100) #Different lengths of time series for which we compare the estimates
different_a      <- c(0.95)
different_slopes <- c(4)
sigma_eta        <- 1



for (a_1 in different_a){
  #true_sigma <- sqrt(sigma_eta^2/((1 - a_1)^2))

  for (T_size in different_T){
    line_trend  <- numeric(T_size)
    
    K1 <- p+1
    K2 <- 10
    L1 <- 20
    L2 <- 30
    
    for (slope in different_slopes){
      line_trend = (1:T_size - 0.5*T_size) * slope/T_size

      a_hat_hall_vect    <- c() 
      a_hat_method1_vect <- c()
      a_hat_method2_vect <- c()
      a_hat_oracle_vect  <- c()
      
      for (i in 1:N_rep){
        set.seed(i)
        
        eps               <-  arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))
        y_data_with_trend <- eps + line_trend
        
        a_hat_method1         <- AR_coefficients(y_data_with_trend, L1, L2, rep(0,L2), p)
        sigma_eta_hat_method1 <- calculating_sigma_eta(y_data_with_trend, a_hat_method1, p)
        sigma_hat_method1     <- sqrt(sigma_eta_hat_method1^2 / (1 - sum(a_hat_method1))^2) 

        corrections_value    <- corrections(a_hat_method1, sigma_eta_hat_method1, K2+1)
        a_hat_method2  <- AR_coefficients(y_data_with_trend, K1, K2, corrections_value, p)
        sigma_eta_hat_method2 <- calculating_sigma_eta(y_data_with_trend, a_hat_method2, p)
        sigma_hat_method2     <- sqrt(sigma_eta_hat_method2^2 / (1 - sum(a_hat_method2))^2) 
        
        a_hat_method1_vect <- c(a_hat_method1_vect, a_hat_method1)
        a_hat_method2_vect <- c(a_hat_method2_vect, a_hat_method2)

        result <- estimating_sigma_for_AR1(y_data_with_trend, L1, L2)
        a_hat_hall_vect   <- c(a_hat_hall_vect, result[[2]])
#          sig_vec_HvK[nb] <- sigma_eta(y_data,a_vec_HvK[nb],p)
#          lrv_vec_HvK[nb] <- sig_vec_HvK[nb]/(1-sum(a_vec_HvK[nb]))^2
          
          res_oracle         <- lm(eps[2:T_size] ~ eps[1:(T_size-1)] - 1)
          a_hat_oracle_vect  <- c(a_hat_oracle_vect, as.vector(res_oracle$coef))
          #sig_vec_oracle[nb] <- mean(res_oracle$residuals^2)
          #lrv_vec_oracle[nb] <- sig_vec_oracle[nb]/(1-sum(a_vec_oracle[nb]))^2 
          
        }
      }
      
      smallest_value <- min(a_hat_hall_vect, a_hat_method1_vect, a_hat_method2_vect, a_hat_oracle_vect)
      biggest_value <- max(a_hat_hall_vect, a_hat_method1_vect, a_hat_method2_vect, a_hat_oracle_vect)
      
      hist1 <- hist(a_hat_oracle_vect, breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), plot = FALSE)
      hist2 <- hist(a_hat_method1_vect, breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), plot = FALSE)
      hist3 <- hist(a_hat_method2_vect, breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), plot = FALSE)
      hist4 <- hist(a_hat_hall_vect, breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), plot = FALSE)
  
      highestCount <- max(hist1$counts, hist2$counts, hist3$counts, hist4$counts)
      dev.new()
      par(mfrow=c(1,4))
      hist(a_hat_oracle_vect, main=paste0("Oracle, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), ylim=c(0,highestCount))
      hist(a_hat_method1_vect, main = paste0("Method 1, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), ylim=c(0,highestCount))
      hist(a_hat_method2_vect, main = paste0("Method 2, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), ylim=c(0,highestCount))
      hist(a_hat_hall_vect, main = paste0("Hall and Van Keilegom, T = ", T_size, ", a_1 = ", a_1), breaks = seq(smallest_value, biggest_value + 0.05, by = 0.01), ylim=c(0,highestCount))
    }
}