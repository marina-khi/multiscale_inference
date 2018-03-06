source("functions.R")
dyn.load("C_code/psihat_statistic.dll")
dyn.load("C_code/estimating_sigma.dll")
source("C_code/psihat_statistic.R")
source("C_code/estimating_sigma.R")

##############################
#Defining necessary constants#
##############################

set.seed(16) #For reproducibility

sigma <- 1 #square root of long run-variance
N <- 1000 #number of simulations
ar1 <- 0.5 #Coefficient for autoregressive model
#T <- 1000 #length 
#alpha <-0.05
#noise_to_signal <- 1
p <- 1 #Order of AR(p) process of the error terms

kernel_f = "epanechnikov" #Only "epanechnikov" and "biweight" kernel functions are currently supported

if (kernel_f == "epanechnikov"){
  kernel_ind = 1
} else if (kernel_f == "biweight"){
  kernel_ind = 2
} else {
  print('Currently only Epanechnikov and Biweight kernel functions are supported')
}

##################################
#Calculating the size of the test#
##################################

size_iid <- c()
size_ar1 <- c()

#We don't have the trend function, the null hypothesis is true and the errors are iid.
for (T in c(250, 500, 1000)){
  g_t_set = creating_g_set(T)
  
  for (alpha in c(0.01, 0.05, 0.1)){
    #Calculating gaussian quantiles for given T and alpha
    gaussian_quantile = calculating_gaussian_quantile(T, g_t_set, kernel_ind, sigma, alpha)
    
    #Replicating test procedure N times
    size_of_the_test_iid = replicate(N, {
      y_data_zero_trend = rnorm(T, 0, sigma)
      result_notrend_iid = psihat_statistic(y_data_zero_trend, g_t_set, kernel_ind, sigma)[[2]]
      if (result_notrend_iid > gaussian_quantile) {d = 1} else {d = 0}
      d
    })
    size_iid <- c(size_iid, sum(size_of_the_test_iid)/N)
    cat("Ratio of rejection in iid case under H0 is ", sum(size_of_the_test_iid)/N, "with T = ", T, "and alpha = ", alpha, "\n")
  }
}
matrix_notrend_iid <- matrix(size_iid, nrow = 3, ncol = 3, byrow = TRUE)

#We don't have the trend function, the null hypothesis is true and the errors are from AR(1)
for (T in c(250, 500, 1000)){
  L1 <- floor(sqrt(T))
  L2 <- floor(2 * sqrt(T))
  g_t_set = creating_g_set(T)
  
  for (alpha in c(0.01, 0.05, 0.1)){
    #Calculating gaussian quantiles for given T and alpha
    gaussian_quantile = calculating_gaussian_quantile(T, g_t_set, kernel_ind, sigma, alpha)

    #Replicating test procedure N times
    size_of_the_test_ar_1 = replicate(N, {
      y_data_ar_1 = arima.sim(model = list(ar = ar1), n = T, innov = rnorm(T, 0, sigma))
      #sigmahat = sigma/(1 - ar1)
      sigmahat = estimating_sigma_for_AR1(y_data_ar_1, L1, L2)
      result_notrend_ar1 = psihat_statistic(y_data_ar_1, g_t_set, kernel_ind, sigmahat)[[2]]
      if (result_notrend_ar1 > gaussian_quantile) {d = 1} else {d = 0}
      d
    })
    size_ar1 <- c(size_ar1, sum(size_of_the_test_ar_1)/N)
    cat("Ratio of rejection in AR(1) under H0 case is ", sum(size_of_the_test_ar_1)/N, "with T = ", T, "and alpha = ", alpha, "\n")
  }
}
matrix_notrend_ar1 <- matrix(size_ar1, nrow = 4, ncol = 3, byrow = TRUE)

###################################
#Calculating the power of the test#
###################################
power_ar1 = c()
#This is for partly linear function + AR(1) case
for (noise_to_signal in c(1, 2, 4)){
  a = sqrt((48 * sigma * sigma)/(5*noise_to_signal))#This is the value of m(1), it depends on Var(e) and Noise to Signal Ratio
  for (T in c(250, 500)){
    L1 = floor(sqrt(T))
    L2 = floor(2 * sqrt(T))
    
    g_t_set = creating_g_set(T)
    
    b = 0.5 * T/a #Constant needed just for calculation
    m = numeric(T)
    for (i in 1:T) {if (i/T < 0.5) {m[i] = 0} else {m[i] = (i - 0.5*T)/b}}
    
    for (alpha in c(0.01, 0.05, 0.1)){
      #Calculating gaussian quantiles for given T and alpha
      gaussian_quantile = calculating_gaussian_quantile(T, g_t_set, kernel_ind, sigma, alpha)
      
      #Replicating test procedure N times
      size_of_the_test_with_trend = replicate(N, {
        #Adding a function that is 0 on the first half and linear on the second half
        y_data_ar_1_with_trend = m + arima.sim(model = list(ar = ar1), innov = rnorm(T, 0, sigma), n = T)
        sigmahat = estimating_sigma_for_AR1(y_data_ar_1_with_trend, L1, L2)
        result_with_trend = psihat_statistic(y_data_ar_1_with_trend, g_t_set, kernel_ind, sigmahat)[[2]]
        if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      power_ar1 = c(power_ar1, sum(size_of_the_test_with_trend)/N)
      cat("Ratio of rejection in AR(1) under H1 case with partially linear trend is", sum(size_of_the_test_with_trend)/N, "with NTS ratio =", noise_to_signal, ", T =", T, "and alpha =", alpha, "\n")
    }
  }
}

power_iid = c()
#This is for partly linear function + AR(1) case
for (noise_to_signal in c(1, 2, 4)){
  a = sqrt((48 * sigma * sigma)/(5 * noise_to_signal))#This is the value of m(1), it depends on Var(e) and Noise to Signal Ratio
  for (T in c(250, 500, 1000)){
    L1 = floor(sqrt(T))
    L2 = floor(2 * sqrt(T))
    
    g_t_set = creating_g_set(T)
    
    b = 0.5 * T/a #Constant needed just for calculation
    m = numeric(T)
    for (i in 1:T) {if (i/T < 0.5) {m[i] = 0} else {m[i] = (i - 0.5*T)/b}}
    
    for (alpha in c(0.01, 0.05, 0.1)){
      #Calculating gaussian quantiles for given T and alpha
      gaussian_quantile = calculating_gaussian_quantile(T, g_t_set, kernel_ind, sigma, alpha)
      
      #Replicating test procedure N times
      size_of_the_test_with_trend_iid = replicate(N, {
        #Adding a function that is 0 on the first half and linear on the second half
        y_data_ar_1_with_trend_iid = m + rnorm(T, 0, sigma)
        result_with_trend_iid = psihat_statistic(y_data_ar_1_with_trend_iid, g_t_set, kernel_ind, sigma)[[2]]
        if (result_with_trend_iid > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      power_iid = c(power_iid, sum(size_of_the_test_with_trend_iid)/N)
      cat("Ratio of rejection in iid case under H1 with partially linear trend is", sum(size_of_the_test_with_trend_iid)/N, "with NTS ratio =", noise_to_signal, ", T =", T, "and alpha =", alpha, "\n")
    }
  }
}

#a_t_set_3 <- plotting_all_rejected_intervals(y_data_3, g_t_set, gaussian_quantile, kernel_ind, sigmahat, "Output/", "nullplot.jpg") #We expect to fail to reject H_0