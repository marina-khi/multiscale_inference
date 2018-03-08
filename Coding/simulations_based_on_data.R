library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("functions.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
dyn.load("C_code/psihat_statistic.dll")
source("C_code/psihat_statistic.R")


N <- 1000 #Number of replications for calculating the size and the power of the test
different_T     <- c(250, 350, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

num_T     <- length(different_T)
num_alpha <- length(different_alpha)

kernel_f = "biweight" #Only "epanechnikov" and "biweight" kernel functions are currently supported

if (kernel_f == "epanechnikov"){
  kernel_ind = 1
} else if (kernel_f == "biweight"){
  kernel_ind = 2
} else {
  print('Currently only Epanechnikov and Biweight kernel functions are supported')
}

#########################################################
#Loading the real data for yearly temperature in England#
#########################################################

temperature             <- read.table("data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr            <- temperature[temperature$YEAR > -99, 'YEAR']
yearly_tempr_normalised <- yearly_tempr - mean(yearly_tempr)

T_tempr <- length(yearly_tempr)

#Tuning parameters
L1 <- floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))

grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating
plot(grid_points, yearly_tempr_normalised, type = 'h')

result <- estimating_sigma_for_AR1(yearly_tempr_normalised, L1, L2)

#sigmahat <- result[[1]]
a_hat <- result[[2]] #Estimation of the AR coefficient
sigma_eta <-result[[3]] #Estimation of the sqrt of the variance 

#Fitting a curve with a data and calculating Noise to Signal Ratio for this curves
#in order to do simulations based on this ratio
noise_to_signal_vector <- c()
for (h in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25)){
  smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr_normalised, grid_points, h))
  # Plotting the estimated curve on top of the data:
  plot(grid_points, yearly_tempr_normalised, type = "l"); lines(grid_points, smoothed_curve, type = "l", col = "red")
  
  #Calculating the integral using the Simpsons method (with fitting a parabola to data)
  weights <- c(1, rep(c(4, 2), times = 178), 4, 1)
  
  integral <- sum(weights * smoothed_curve) / (3 * T_tempr)
  integral_squared_function <- sum(weights * smoothed_curve * smoothed_curve) / (3 * T_tempr)
  variance <- integral_squared_function - integral * integral
  noise_to_signal <- sigma_eta / sqrt(variance)
  noise_to_signal_vector <- c(noise_to_signal_vector, noise_to_signal)
  cat("noise_to_signal:", noise_to_signal, "\n")
}

##################################
#Calculating the size of the test#
##################################

#We don't have the trend function, the null hypothesis is true and the errors are from AR(1)
size_ar1 <- c()
for (T in different_T){
  L1 <- floor(sqrt(T))
  L2 <- floor(2 * sqrt(T))
  g_t_set = creating_g_set(T)
  
  for (alpha in different_alpha){
    #Calculating gaussian quantiles for given T and alpha
    gaussian_quantile = calculating_gaussian_quantile(T, g_t_set, kernel_ind, sigma_eta, alpha)
    
    #Replicating test procedure N times
    size_of_the_test_ar_1 = replicate(N, {
      y_data_ar_1 = arima.sim(model = list(ar = a_hat), n = T, innov = rnorm(T, 0, sigma_eta))
      sigmahat = estimating_sigma_for_AR1(y_data_ar_1, L1, L2)
      result_notrend_ar1 = psihat_statistic(y_data_ar_1, g_t_set, kernel_ind, sigmahat)[[2]]
      if (result_notrend_ar1 > gaussian_quantile) {d = 1} else {d = 0}
      d
    })
    size_ar1 <- c(size_ar1, sum(size_of_the_test_ar_1)/N)
    cat("Ratio of rejection in AR(1) under H0 case is ", sum(size_of_the_test_ar_1)/N, "with T = ", T, "and alpha = ", alpha, "\n")
  }
}
#Creating a nice-looking matrix for size of the test and writing it to a tex file
creating_matrix_and_texing(size_ar1, different_T, different_alpha, "Output/sizetable.tex")

###################################
#Calculating the power of the test#
###################################

power_ar1 = c()
#This is for partly linear function + AR(1) case
for (nts in c(2, 3, 4)){
  a = sqrt((48 * sigma_eta * sigma_eta)/(5*nts))#This is the value of m(1), it depends on Var(e) and Noise to Signal Ratio
  for (T in different_T){
    L1 <- floor(sqrt(T))
    L2 <- floor(2 * sqrt(T))
    g_t_set = creating_g_set(T)
    
    b = 0.5 * T/a #Constant needed just for calculation
    m = numeric(T)
    for (i in 1:T) {if (i/T < 0.5) {m[i] = 0} else {m[i] = (i - 0.5*T)/b}}
    
    for (alpha in different_alpha){
      #Calculating gaussian quantiles for given T and alpha
      gaussian_quantile = calculating_gaussian_quantile(T, g_t_set, kernel_ind, sigma_eta, alpha)
      
      #Replicating test procedure N times
      size_of_the_test_with_trend = replicate(N, {
        #Adding a function that is 0 on the first half and linear on the second half
        y_data_ar_1_with_trend = m + arima.sim(model = list(ar = a_hat), innov = rnorm(T, 0, sigma_eta), n = T)
        sigmahat = estimating_sigma_for_AR1(y_data_ar_1_with_trend, L1, L2)
        result_with_trend = psihat_statistic(y_data_ar_1_with_trend, g_t_set, kernel_ind, sigmahat)[[2]]
        if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
        d
      })
      power_ar1 = c(power_ar1, sum(size_of_the_test_with_trend)/N)
      cat("Ratio of rejection in AR(1) under H1 case with partially linear trend is", sum(size_of_the_test_with_trend)/N, "with NTS ratio =", nts, ", T =", T, "and alpha =", alpha, "\n")
    }
  }
}

#Creating nice-looking matrices for power of the test (each for different Noise to Signal ratio) and writing them to tex files
creating_matrix_and_texing(power_ar1[1:(num_T * num_alpha)], different_T, different_alpha, "Output/powertable_nts2.tex")
creating_matrix_and_texing(power_ar1[(num_T * num_alpha + 1):(2 * num_T * num_alpha)], different_T, different_alpha, "Output/powertable_nts3.tex")
creating_matrix_and_texing(power_ar1[(2 * num_T * num_alpha + 1):(3 * num_T * num_alpha)], different_T, different_alpha, "Output/powertable_nts4.tex")
