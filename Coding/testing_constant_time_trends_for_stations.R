source("functions.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
dyn.load("C_code/psihat_statistic.dll")
source("C_code/psihat_statistic.R")

##############################
#Defining necessary constants#
##############################

p <- 1 #Order of AR(p) process of the error terms. Currently only p=1 is supported
alpha <-0.05
N <- 34
kernel_f = "biweight" #Only "epanechnikov" and "biweight" kernel functions are currently supported

if (kernel_f == "epanechnikov"){
  kernel_ind = 1
} else if (kernel_f == "biweight"){
  kernel_ind = 2
} else {
  print('Currently only Epanechnikov and Biweight kernel functions are supported')
}

###############################################
#Loading the real data for stations in England#
###############################################

for (i in 1:N){
  filename = paste("data/txt", i, ".txt", sep = "")
  temperature_tmp  <- read.table(filename, header = FALSE, skip = 3,
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


######################
#Deseasonalizing data#
######################

TemperatureColumns <- setdiff(names(monthly_temp), c("year", "month"))
for (i in TemperatureColumns){
  monthly_temp[[i]] <- monthly_temp[[i]] - ave(monthly_temp[[i]], monthly_temp[['month']], FUN = mean)
  #plot(monthly_temp[[i]],type = 'l')
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


#####################################
#Calculating gaussian quantile for T#
#####################################

g_t_set_tempr <- creating_g_set(T_tempr)
gaussian_quantile <- calculating_gaussian_quantile(T_tempr, g_t_set_tempr, kernel_ind, sigmahat_tempr, alpha)


#########################################
#Calculating the statistic for real data#
#########################################

for (i in 1:N){
  psihat_statistic_value <- psihat_statistic(monthly_temp[[i+2]], g_t_set_tempr, kernel_ind, sqrt(sigmahat_vector_2[i]))[[2]]
  if (psihat_statistic_value > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  } else {
  cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
      "Gaussian quantile value = ", gaussian_quantile, "\n")
  #plot(monthly_temp[[i]],type = 'l')
  }
}