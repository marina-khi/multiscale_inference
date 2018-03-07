source("functions.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
dyn.load("C_code/psihat_statistic.dll")
source("C_code/psihat_statistic.R")

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
N <- 1000

#Tuning parameters
L1 <- floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))

grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr)
plot(grid_points, yearly_tempr_normalised, type = 'h')

result <- estimating_sigma_for_AR1(yearly_tempr_normalised, L1, L2)

#sigmahat <- result[[1]]
a_hat <- result[[2]]
sigma_eta <-result[[3]]

noise_to_signal_vector <- c()
for (h in c(0.01, 0.05, 0.1, 0.15, 0.2, 0.25)){
  nonpar.reg <- ksmooth(grid_points, yearly_tempr_normalised, kernel="normal", bandwidth=h)
  
  # Plotting the estimated curve on top of the data:
  plot(grid_points, yearly_tempr_normalised, type = "l"); lines(nonpar.reg)
  
  #Calculating the integral using the Simpsons method (with fitting a parabola to data)
  weights <- c(1, rep(c(4, 2), times = 178), 4, 1)
  
  integral <- sum(weights * nonpar.reg$y) / (3 * T_tempr)
  integral_squared_function <- sum(weights * nonpar.reg$y * nonpar.reg$y) / (3 * T_tempr)
  variance <- integral_squared_function - integral * integral
  noise_to_signal <- sigma_eta / sqrt(variance)
  noise_to_signal_vector <- c(noise_to_signal_vector, noise_to_signal)
  cat("noise_to_signal:", noise_to_signal, "\n")
}
T_tempr <-300
power_ar1 = c()
#This is for partly linear function + AR(1) case
for (nts in c(2, 4)){
  a = sqrt((48 * sigma_eta * sigma_eta)/(5*nts))#This is the value of m(1), it depends on Var(e) and Noise to Signal Ratio
  g_t_set = creating_g_set(T_tempr)
  b = 0.5 * T_tempr/a #Constant needed just for calculation
  m = numeric(T_tempr)
  for (i in 1:T_tempr) {if (i/T_tempr < 0.5) {m[i] = 0} else {m[i] = (i - 0.5*T_tempr)/b}}
  for (alpha in c(0.01, 0.05, 0.1)){
    #Calculating gaussian quantiles for given T and alpha
    gaussian_quantile = calculating_gaussian_quantile(T_tempr, g_t_set, kernel_ind, sigma_eta, alpha)

    #Replicating test procedure N times
    size_of_the_test_with_trend = replicate(N, {
      #Adding a function that is 0 on the first half and linear on the second half
      y_data_ar_1_with_trend = m + arima.sim(model = list(ar = a_hat), innov = rnorm(T_tempr, 0, sigma_eta), n = T_tempr)
      sigmahat = estimating_sigma_for_AR1(y_data_ar_1_with_trend, L1, L2)
      result_with_trend = psihat_statistic(y_data_ar_1_with_trend, g_t_set, kernel_ind, sigmahat)[[2]]
      if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
      d
    })
    power_ar1 = c(power_ar1, sum(size_of_the_test_with_trend)/N)
    cat("Ratio of rejection in AR(1) under H1 case with partially linear trend is", sum(size_of_the_test_with_trend)/N, "with NTS ratio =", nts, ", T =", T_tempr, "and alpha =", alpha, "\n")
  }
}

# weights_trapezoid <- c(1, rep(2, times = 357), 1)
# 
# integral_t <- sum(weights_trapezoid * nonpar.reg$y) / (3 * T_tempr)
# integral_squared_function_t <- sum(weights_trapezoid * nonpar.reg$y * nonpar.reg$y) / (3 * T_tempr)
# variance_t <- integral_squared_function_t - integral_t * integral_t
# noise_to_signal_t <- sigma_eta / sqrt(variance_t)