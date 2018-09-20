library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/C_code/estimating_sigma.R")
dyn.load("Shape/C_code/estimating_sigma.dll")
dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/data_analysis.R")


###############################
#Defining necessary parameters#
###############################

alpha <- 0.05 #alpha for calculating quantiles
h     <- c(0.05, 0.1, 0.15, 0.2) #Different bandwidth for plotting. Number must be <=6 in order for the plot to be readable

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported
test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 

#############################
#Loading the data for income#
#############################

income <- read.table("Shape/data/France_income_top10.csv", header = FALSE, dec = ",")
yearly_income <- income[, 1]
T_income      <- length(yearly_income)


#####################################
#Estimating parameters from the data#
#####################################

#Tuning parameters
L1_i     <- floor(sqrt(T_income))
L2_i     <- floor(2 * sqrt(T_income))
result_i <- estimating_sigma_for_AR1(yearly_income, L1_i, L2_i)

a_hat_i     <- result_i[[2]] #Estimation of the AR coefficient
sigma_eta_i <- result_i[[3]] #Estimation of the sqrt of the variance of the innovation 
sigmahat_i  <- result_i[[1]]

###################################################################
#Calculating smoothed curve for the data using Epanechnikov kernel#
###################################################################

grid_points <- seq(from = 1/T_income, to = 1, length.out = T_income) #grid points for plotting and estimating
grid_time   <- seq(from = 1908, to = 2010, length.out = T_income) #grid points for plotting 


plot(NA, xlim = c(1908, 2010), ylim = c(0.25, 0.5), mgp=c(2,0.5,0))

for (i in 1:length(h)){
  if (kernel_method == "nw"){
    smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_income, grid_points, h[i]))
  } else if (kernel_method == "ll"){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(yearly_income, grid_points, h[i]))
  } else {
    print('Given method is currently not supported')
  }
  lines(grid_time, smoothed_curve, lty = i) 
}

###############
#Data analysis#
###############

if (test_problem == "zero"){
  kernel_ind = 1
} else if (test_problem == "constant"){
  kernel_ind = 2
} else {
  print('Given testing problem is currently not supported')
}
  
if (kernel_method == "nw"){
  quantile_function = calculating_gaussian_quantile
  statistic_function = psihat_statistic
  smoothing = epanechnikov_smoothing
} else if (kernel_method == "ll"){
  quantile_function = calculating_gaussian_quantile_ll
  statistic_function = psihat_statistic_ll
  smoothing = local_linear_smoothing
} else {
  print('Given method is currently not supported')
}
  
#Calculating gaussian quantile for T_data#

g_t_set           <- creating_g_set(T_income, kernel_method)
gaussian_quantile <- quantile_function(T_income, g_t_set, test_problem, kernel_ind, alpha)
  

result                 <- statistic_function(yearly_income, g_t_set, kernel_ind, sigmahat_i)
g_t_set_with_values    <- result[[1]]
psihat_statistic_value <- result[[2]]
  
#And now the testing itself
if (psihat_statistic_value > gaussian_quantile) {
  cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
      "Gaussian quantile value = ", gaussian_quantile, "\n")
  
  #The collection of intervals where the corrected test statistic lies above the critical value (as in (2.6))
  a_t_set <- subset(g_t_set_with_values, values > gaussian_quantile, select = c(u, h, values))
  p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*T_income + 1908, 'endpoint' = (a_t_set$u + a_t_set$h)*T_income + 1908, 'values' = a_t_set$values)
  p_t_set <- subset(p_t_set, endpoint <= 2017, select = c(startpoint, endpoint, values)) 
  p_t_set <- choosing_minimal_intervals(p_t_set)
  
} else {
  cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
      "Gaussian quantile value = ", gaussian_quantile, "\n")
}
  
  #Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
  ymaxlim = max(p_t_set$values)
  yminlim = min(p_t_set$values)
  plot(NA, xlim=c(1908,2010), ylim = c(yminlim - 0.5, ymaxlim + 0.5), mgp=c(2,0.5,0))
  segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set$endpoint, p_t_set[['values']])
  abline(h = gaussian_quantile, lty = 2)