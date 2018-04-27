library(car)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("C_code/psihat_statistic.R")
dyn.load("C_code/psihat_statistic_ll.dll")
dyn.load("C_code/psihat_statistic_nw.dll")

source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
source("functions.R")
source("testing_different_time_trends.R")
source("simulations_based_on_data.R")

##############################
#Defining necessary constants#
##############################

N_ts     <- 34 #number of different time series for application
N_ts_sim <- 15 #number of different time series for simulation
N_rep    <- 1000 #number of repetitions for calculating size and power
alpha    <- 0.05 #alpha for calculating quantiles

different_T     <- c(250, 300, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported

a_hat <- 0.267 #Parameters that are used in simulations but were estimated beforehand in apllication (shape)
sigma <- 0.59

###########################################
#Loading the real station data for England#
###########################################

for (i in 1:N_ts){
  filename = paste("data/txt", i, ".txt", sep = "")
  temperature_tmp  <- read.table(filename, header = FALSE, skip = 7,
                                 col.names = c("year", "month", "tmax", "tmin", "af", "rain", "sun", "aux"), fill = TRUE,  na.strings = c("---"))
  monthly_temp_tmp <- data.frame('1' = as.numeric(temperature_tmp[['year']]), '2' = as.numeric(temperature_tmp[['month']]),
                                 '3' = (temperature_tmp[["tmax"]] + temperature_tmp[["tmin"]]) / 2)
  colnames(monthly_temp_tmp) <- c('year', 'month', paste0("tmean", i))


  if (i == 1){
    monthly_temp <- monthly_temp_tmp
  } else {
    monthly_temp <- merge(monthly_temp, monthly_temp_tmp, by = c("year", "month"), all.x = TRUE, all.y = TRUE)
  }
}

monthly_temp <- subset(monthly_temp, year >= 1986)
monthly_temp <- monthly_temp[,colSums(is.na(monthly_temp)) <= 2]
monthly_temp <- na.omit(monthly_temp)#Deleting the rows with ommitted variables

date         <- paste(sprintf("%02d", monthly_temp$month), monthly_temp$year,  sep='-')
monthly_temp <- cbind(date, monthly_temp)
T_tempr      <- nrow(monthly_temp)
N_ts         <- ncol(monthly_temp) - 3

######################
#Deseasonalizing data#
######################


TemperatureColumns <- setdiff(names(monthly_temp), c("year", "month", "date"))
monthly_temp[4:(N_ts + 3)] <- lapply(monthly_temp[4:(N_ts + 3)], function(x) x - ave(x, monthly_temp[['month']], FUN=mean))


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


#################################
#Testing equality of time trends#
#################################

results <- testing_different_time_trends(N_ts, monthly_temp[-c(1, 2, 3)], monthly_temp['month'], alpha, kernel_method, sigmahat_vector_2)


############################
#Calculating size and power#
############################

results_size     <- simulations_size(a_hat, sigma, N_ts_sim, N_rep, different_T, different_alpha, kernel_method)
results_power    <- simulations_power(a_hat, sigma, N_ts_sim, different_T, different_alpha, kernel_method)
simulations_clustering(a_hat, sigma, N_ts_sim, N_rep, different_T, different_alpha, kernel_method)
results_clusters <- clustering_analysis(N_ts_sim, different_T, different_alpha)