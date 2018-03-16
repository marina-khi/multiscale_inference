dyn.load("C_code/psihat_statistic_ij.dll")
source("C_code/psihat_statistic_ij.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
source("functions.R")

##############################
#Defining necessary constants#
##############################

p <- 1 #Order of AR(p) process of the error terms. Currently only p = 1 is supported
N <- 34 #number of different time series 

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
  filename = paste("data/txt", i, ".txt", sep = "")
  temperature_tmp  <- read.table(filename, header = FALSE, skip = 7,
                                 col.names = c("year", "month", "tmax", "tmin", "af", "rain", "sun", "aux"), fill = TRUE,  na.strings = c("---"))
  monthly_temp_tmp <- data.frame('1' = temperature_tmp[['year']], '2' = temperature_tmp[['month']],
                                 #'3' = (temperature_tmp[["tmax"]] + temperature_tmp[["tmin"]]) / 2)
                                 '3' = temperature_tmp[["tmax"]])
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

# for (column in TemperatureColumns){
#   plot(NA, xlim = c(0, 1), ylim = c(-3, 3))
#   for (i in 1:4){
#     #This part plots kernel smoothers for different bandwidths (all on one plot).
#     smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, h[i]))
#     lines(grid_points, smoothed_curve, lty = i)#, col = colors[i]) 
#   }
# }

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


###########################################
#Calculating gaussian quantile for T_tempr#
###########################################

g_t_set_tempr <- creating_g_set(T_tempr)
gaussian_quantile <- calculating_gaussian_quantile_ij(T_tempr, N, g_t_set_tempr, kernel_ind, alpha)
pairwise_gaussian_quantile <- calculating_gaussian_quantile_ij(T_tempr, 2, g_t_set_tempr, kernel_ind, alpha)

#########################################
#Calculating the statistic for real data#
#########################################

matrix_of_statistic <- matrix(, nrow = N, ncol = N)
matrix_of_statistic_new <- matrix(, nrow = N, ncol = N)
for (i in (1 : (N - 1))){
  for (j in ((i + 1):N)){
    sigmahat_new = sqrt(sigmahat_vector_2[i] + sigmahat_vector_2[j])
    result = psihat_statistic_ij(monthly_temp[[i + 2]], monthly_temp[[j + 2]], g_t_set_tempr, kernel_ind, sqrt(2) * sigmahat_tempr)
    result_new = psihat_statistic_ij(monthly_temp[[i + 2]], monthly_temp[[j + 2]], g_t_set_tempr, kernel_ind, sigmahat_new)
    matrix_of_statistic[i, j] = result[[2]]
    matrix_of_statistic_new[i, j] = result_new[[2]]
    cat(matrix_of_statistic[i, j], "with i =", i, "and j = ",j, "\n")
    g_t_set = result[[1]]
    #cat(head(monthly_temp[[i + 2]]), "\n")
    #cat(head(monthly_temp[[j + 2]]), "\n")
    }
}
statistic = max(as.vector(matrix_of_statistic), na.rm=TRUE)
statistic_new = max(as.vector(matrix_of_statistic_new), na.rm=TRUE)
