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

N_ts     <- 34 #number of different time series 
N_rep    <- 1000 #number of repetitions for calculating size and power
alpha    <- 0.05 #alpha for calculating quantiles

different_T     <- c(250, 300, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported


###########################################
#Loading the real station data for England#
###########################################

for (i in 1:N_ts){
  filename = paste("data/txt", i, ".txt", sep = "")
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
date         <- paste(sprintf("%02d", monthly_temp$month), monthly_temp$year,  sep='-')
monthly_temp <- cbind(date, monthly_temp)
T_tempr      <- nrow(monthly_temp)


######################
#Deseasonalizing data#
######################

TemperatureColumns <- setdiff(names(monthly_temp), c("year", "month", "date"))
for (i in TemperatureColumns){
  monthly_temp[[i]] <- monthly_temp[[i]] - ave(monthly_temp[[i]], monthly_temp[['month']], FUN = mean)
  #plot(monthly_temp[[i]],type = 'l')
}


###################################################################
#Calculating smoothed curve for the data using Epanechnikov kernel#
###################################################################

# grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for estimating
# 
# pdf("../Plots/stations_data.pdf", width=10, height=10, paper="special")
# par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
# par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
# par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins
# 
# plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(-1.5, 1.5), yaxp  = c(-1.0, 1.0, 2), xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
# for (column in TemperatureColumns){
#   smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, 0.05))
#   lines(grid_points, smoothed_curve)#, lty = i), col = colors[i]) 
# }
# axis(1, at = grid_points[seq(1, 300, by = 20)], labels = monthly_temp$date[seq(1, 300, by = 20)])
# legend(0, 1.0, legend=c("h = 0.05"), lty = 1, cex = 0.95, ncol=1)
# 
# plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(-1.5, 1.5), yaxp  = c(-1.0, 1.0, 2),xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
# for (column in TemperatureColumns){
#   smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, 0.1))
#   lines(grid_points, smoothed_curve)#, lty = i), col = colors[i]) 
# }
# axis(1, at = grid_points[seq(1, 300, by = 20)], labels = monthly_temp$date[seq(1, 300, by = 20)])
# legend(0, 1.0, legend=c("h = 0.10"), lty = 1, cex = 0.95, ncol=1)
# 
# 
# plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(-1.5, 1.5), yaxp  = c(-1.0, 1.0, 2),xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
# for (column in TemperatureColumns){
#   smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, 0.15))
#   lines(grid_points, smoothed_curve)#, lty = i), col = colors[i]) 
# }
# axis(1, at = grid_points[seq(1, 300, by = 20)], labels = monthly_temp$date[seq(1, 300, by = 20)])
# legend(0, 1.0, legend=c("h = 0.15"), lty = 1, cex = 0.95, ncol=1)
# 
# dev.off()


#####################
#Estimating variance#
#####################

#Tuning parameters
# L1 <- floor(sqrt(T_tempr))
# L2 <- floor(2 * sqrt(T_tempr))
# 
# #Calculating each sigma_i separately
# sigmahat_vector_2 <- c()
# for (i in TemperatureColumns){
#   sigma_i <- estimating_sigma_for_AR1(monthly_temp[[i]], L1, L2)[[1]]
#   sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_i * sigma_i)
# }

#Calculating sqrt root of the average long-run variance
#sigmahat_tempr <- sqrt(sum(sigmahat_vector_2)/N)


#################################
#Testing equality of time trends#
#################################

#results <- testing_different_time_trends(N_ts, monthly_temp[-c('year', 'date', 'month')], alpha, kernel_method, sigmahat_vector_2)


############################
#Calculating size and power#
############################

source("simulations_based_on_data.R")
#results_size     <- simulations_size(15, N_rep, different_T, different_alpha, kernel_method)
results_power    <- simulations_power(15, N_rep, different_T, different_alpha, kernel_method)
#results_clusters <- simulations_clustering(15, N_rep, different_T, different_alpha, kernel_method)

