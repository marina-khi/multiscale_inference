#library(fpp)
##############################
#Defining necessary constants#
##############################

p <- 1 #Order of AR(p) process of the error terms. Currently only p=1 is supported
alpha <-0.05
N <- 3 #number of different time series 
kernel_f = "biweight" #Only "epanechnikov" and "biweight" kernel functions are currently supported

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

temperature_1  <- read.table("data/1_station.txt", header = FALSE, skip = 7,
                            col.names = c("year", "month", "tmax", "tmin", "af", "rain", "sun", "aux"), fill = TRUE,  na.strings = c("---"))
monthly_temp_1 <- data.frame('year' = temperature_1[['year']], 'month' = temperature_1[['month']],
                      'tmean1' = (temperature_1[["tmax"]] + temperature_1[["tmin"]]) / 2)

temperature_2 <- read.table("data/2_station.txt", header = FALSE, skip = 7,
                            col.names = c("year", "month", "tmax", "tmin", "af", "rain", "sun", "aux"), fill = TRUE, na.strings = c("---"))
monthly_temp_2 <- data.frame('year' = temperature_2[['year']], 'month' = temperature_2[['month']],
                             'tmean2' = (temperature_2[["tmax"]] + temperature_2[["tmin"]]) / 2)

temperature_3 <- read.table("data/3_station.txt", header = FALSE, skip = 7,
                            col.names = c("year", "month", "tmax", "tmin", "af", "rain", "sun", "aux"), fill = TRUE,  na.strings = c("---"))
monthly_temp_3 <- data.frame('year' = temperature_3[['year']], 'month' = temperature_3[['month']],
                             'tmean3' = (temperature_3[["tmax"]] + temperature_3[["tmin"]]) / 2)

monthly_temp <- merge(monthly_temp_1, monthly_temp_2, by = c("year", "month"))
monthly_temp <- merge(monthly_temp, monthly_temp_3, by = c("year", "month"))

monthly_temp <- na.omit(monthly_temp)#Deleting the rows with ommitted variables

plot(monthly_temp[["tmean1"]],type = 'l')
plot(monthly_temp[["tmean2"]],type = 'l')
plot(monthly_temp[["tmean3"]],type = 'l')


######################
#Deseasonalizing data#
######################

for 


yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
yearly_tempr_normalised <- yearly_tempr - mean(yearly_tempr)

T_tempr <- length(yearly_tempr)

#Tuning parameters
L1 <- floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))
grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating

sigmahat_tempr <- estimating_sigma_for_AR1(yearly_tempr_normalised, L1, L2)[[1]]

