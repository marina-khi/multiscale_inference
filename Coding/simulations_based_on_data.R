library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("functions.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
dyn.load("C_code/psihat_statistic.dll")
source("C_code/psihat_statistic.R")


###############################
#Defining necessary parameters#
###############################

N <- 1000 #Number of replications for calculating the size and the power of the test
different_T     <- c(250, 350, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

num_T     <- length(different_T)
num_alpha <- length(different_alpha)

h         <- c(0.05, 0.1, 0.15, 0.2) #Different bandwidth for plotting. Number must be <=6
colors    <- c("black", "red", "green", "blue", "yellow") #Different colors for plotting. Number must be >= number of bandwidths

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
#Plotting the normalized data
grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating
plot(grid_points, yearly_tempr_normalised, type = 'h')


#####################################
#Estimating parameters from the data#
#####################################

#Tuning parameters
L1 <- floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))
result <- estimating_sigma_for_AR1(yearly_tempr_normalised, L1, L2)

#sigmahat <- result[[1]]
a_hat <- result[[2]] #Estimation of the AR coefficient
sigma_eta <-result[[3]] #Estimation of the sqrt of the variance 


###################################################################
#Calculating smoothed curve for the data using Epanechnikov kernel#
###################################################################

#Fitting a curve with a data and calculating Noise to Signal Ratio for this curves
#in order to do simulations based on this ratio

end_point <- c()
plot(NA, xlim = c(0, 1), ylim = c(-1.5, 1.5))
for (i in 1:4){
  #This part plots kernel smoothers for different bandwidths (all on one plot).
  smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr_normalised, grid_points, h[i]))
  lines(grid_points, smoothed_curve, lty = i, col = colors[i])
  end_point <- c(end_point, smoothed_curve[T_tempr])
  cat("End point:", smoothed_curve[T_tempr], "\n")
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

plot(NA, xlim = c(0, 1), ylim = c(-1.5, 1.5))

#This is for partly linear function + AR(1) case
for (a in c(1.0, 0.75, 0.65, 0.5, 0.25)){
  power_ar1 = c()
  #a = sqrt((48 * sigma_eta * sigma_eta)/(5*nts))#This is the value of m(1), it depends on Var(e) and Noise to Signal Ratio
  for (T in different_T){
    L1 <- floor(sqrt(T))
    L2 <- floor(2 * sqrt(T))
    g_t_set = creating_g_set(T)
    grid_points <- seq(from = 1/T, to = 1, length.out = T)
    
    b = 0.4 * T/a #Constant needed just for calculation
    m = numeric(T)

    for (i in 1:T) {if (i/T < 0.6) {m[i] = 0} else {m[i] = (i - 0.6*T)/b}}
    lines(grid_points, m, lty = i)#, col = colors[i])
    
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
      #cat("Ratio of rejection in AR(1) under H1 case with partially linear trend is", sum(size_of_the_test_with_trend)/N, "with a =", a, ", T =", T, "and alpha =", alpha, "\n")
    }
  }
  filename = paste0("Output/powertable_", a*100, ".tex")
  creating_matrix_and_texing(power_ar1, different_T, different_alpha, filename)
}

#creating_matrix_and_texing(power_ar1[1:(num_T * num_alpha)], different_T, different_alpha, "Output/powertable_v1.tex")
#creating_matrix_and_texing(power_ar1[(num_T * num_alpha + 1):(2 * num_T * num_alpha)], different_T, different_alpha, "Output/powertable_v2.tex")
#creating_matrix_and_texing(power_ar1[(2 * num_T * num_alpha + 1):(3 * num_T * num_alpha)], different_T, different_alpha, "Output/powertable_v3.tex")