dyn.load("C_code/psihat_statistic_ij.dll")
source("C_code/psihat_statistic_ij.R")
dyn.load("C_code/psihat_statistic_ij_ll.dll")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
source("functions.R")
source("testing_different_time_trends.R")

##############################
#Defining necessary constants#
##############################

N        <- 34 #number of different time series 
N_rep    <- 1000 #number of repetitions for calculating gaussian statistic
alpha    <- 0.05 #alpha for calculating quantiles

different_T     <- c(250, 350, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported


###########################################
#Loading the real station data for England#
###########################################

for (i in 1:N){
  filename = paste("data/txt", i, ".txt", sep = "")
  temperature_tmp  <- read.table(filename, header = FALSE, skip = 7,
                                 col.names = c("year", "month", "tmax", "tmin", "af", "rain", "sun", "aux"), fill = TRUE,  na.strings = c("---"))
  monthly_temp_tmp <- data.frame('1' = temperature_tmp[['year']], '2' = temperature_tmp[['month']],
                                 '3' = (temperature_tmp[["tmax"]] + temperature_tmp[["tmin"]]) / 2)
                                 #'3' = temperature_tmp[["tmax"]])
  colnames(monthly_temp_tmp) <- c('year', 'month', paste0("tmean", i))
  if (i == 1){
    monthly_temp <- monthly_temp_tmp
  } else {
    monthly_temp <- merge(monthly_temp, monthly_temp_tmp, by = c("year", "month"))
  }
}

monthly_temp <- na.omit(monthly_temp)#Deleting the rows with ommitted variables
T_tempr <- nrow(monthly_temp)

######################
#Deseasonalizing data#
######################

TemperatureColumns <- setdiff(names(monthly_temp), c("year", "month"))
for (i in TemperatureColumns){
  monthly_temp[[i]] <- monthly_temp[[i]] - ave(monthly_temp[[i]], monthly_temp[['month']], FUN = mean)
  #plot(monthly_temp[[i]],type = 'l')
}

###################################################################
#Calculating smoothed curve for the data using Epanechnikov kernel#
###################################################################

h <- c(0.01, 0.05, 0.1, 0.15)
grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating

plot(NA, xlim = c(0, 1), ylim = c(-1.5, 1.5))
for (column in TemperatureColumns){
  smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, 0.05))
  lines(grid_points, smoothed_curve)#, lty = i), col = colors[i]) 
}

plot(NA, xlim = c(0, 1), ylim = c(-1.5, 1.5))
for (column in TemperatureColumns){
  smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, 0.1))
  lines(grid_points, smoothed_curve)#, lty = i), col = colors[i]) 
}


#####################
#Estimating variance#
#####################

#Tuning parameters
L1 <- floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))

#Calculating each sigma_i separately
sigmahat_vector_2 <- c()
for (i in TemperatureColumns){
  sigma_i <- estimating_sigma_for_AR1(monthly_temp[[i]], L1, L2)[[1]]
  sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_i * sigma_i)
}
#Calculating sqrt root of the average long-run variance
sigmahat_tempr <- sqrt(sum(sigmahat_vector_2)/N)


#################################
#Testing equality of time trends#
#################################

results <- testing_different_time_trends(N, monthly_temp, alpha, kernel_method, sigmahat_vector_2)