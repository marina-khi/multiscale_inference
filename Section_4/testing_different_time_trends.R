dyn.load("C_code/psihat_statistic_ij.dll")
source("C_code/psihat_statistic_ij.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
source("functions.R")

##############################
#Defining necessary constants#
##############################

p <- 1 #Order of AR(p) process of the error terms. Currently only p=1 is supported
N <- 3 #number of different time series 

alpha    <- 0.05
kernel_f <- "epanechnikov" #Only "epanechnikov" and "biweight" kernel functions are currently supported

if (kernel_f == "epanechnikov"){
  kernel_ind = 1
} else if (kernel_f == "biweight"){
  kernel_ind = 2
} else {
  print('Currently only Epanechnikov and Biweight kernel functions are supported')
}

###########################################
#Loading the real station data for England#
###########################################

for (i in 1:N){
  filename = paste("data/", i, "_station.txt", sep = "")
  temperature_tmp  <- read.table(filename, header = FALSE, skip = 7,
                                 col.names = c("year", "month", "tmax", "tmin", "af", "rain", "sun", "aux"), fill = TRUE,  na.strings = c("---"))
  monthly_temp_tmp <- data.frame('1' = temperature_tmp[['year']], '2' = temperature_tmp[['month']],
                                 '3' = (temperature_tmp[["tmax"]] + temperature_tmp[["tmin"]]) / 2)
  colnames(monthly_temp_tmp) <- c('year', 'month', paste0("tmean", i))
  if (i == 1){
    monthly_temp <- monthly_temp_tmp
  } else {
    monthly_temp <- merge(monthly_temp, monthly_temp_tmp, by = c("year", "month"))
  }
}

monthly_temp <- na.omit(monthly_temp)#Deleting the rows with ommitted variables
T_tempr <- nrow(monthly_temp)

plot(monthly_temp[["tmean1"]],type = 'l')
plot(monthly_temp[["tmean2"]],type = 'l')
plot(monthly_temp[["tmean3"]],type = 'l')


######################
#Deseasonalizing data#
######################

TemperatureColumns <- setdiff(names(monthly_temp), c("year", "month"))
for (i in TemperatureColumns){
  monthly_temp[[i]] <- monthly_temp[[i]] - ave(monthly_temp[[i]], monthly_temp[['month']], FUN = mean)
  plot(monthly_temp[[i]],type = 'l')
}

#####################
#Estimating variance#
#####################

#Tuning parameters
L1 <- floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))

#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in TemperatureColumns){
  sigma_i <- estimating_sigma_for_AR1(monthly_temp[[i]], L1, L2)[[1]]
  sigmahat_vector <- c(sigmahat_vector, sigma_i * sigma_i)
}

sigmahat_tempr <- sqrt(sum(sigmahat_vector)/N)

###########################################
#Calculating gaussian quantile for T_tempr#
###########################################

g_t_set_tempr <- creating_g_set(T_tempr)
gaussian_quantile <- calculating_gaussian_quantile(T_tempr, N, g_t_set_tempr, kernel_ind, sigmahat_tempr, alpha)


grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating

