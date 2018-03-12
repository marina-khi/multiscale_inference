source("functions.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
dyn.load("C_code/psihat_statistic.dll")
source("C_code/psihat_statistic.R")

##############################
#Defining necessary constants#
##############################

noise_to_signal <- 1
p <- 1 #Order of AR(p) process of the error terms. Currently only p=1 is supported
alpha <-0.05
kernel_f = "biweight" #Only "epanechnikov" and "biweight" kernel functions are currently supported

if (kernel_f == "epanechnikov"){
  kernel_ind = 1
} else if (kernel_f == "biweight"){
  kernel_ind = 2
} else {
  print('Currently only Epanechnikov and Biweight kernel functions are supported')
}

###################################
#Loading the real data for England#
###################################

temperature <- read.table("data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
yearly_tempr_normalised <- yearly_tempr - mean(yearly_tempr)

T_tempr <- length(yearly_tempr)

#Tuning parameters
L1 <- floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))
grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating

sigmahat_tempr <- estimating_sigma_for_AR1(yearly_tempr_normalised, L1, L2)[[1]]


#####################################
#Calculating gaussian quantile for T#
#####################################

g_t_set_tempr <- creating_g_set(T_tempr)
gaussian_quantile <- calculating_gaussian_quantile(T_tempr, g_t_set_tempr, kernel_ind, sigmahat_tempr, alpha)


#########################################
#Calculating the statistic for real data#
#########################################

results <- plotting_all_rejected_intervals(yearly_tempr, g_t_set_tempr, gaussian_quantile, kernel_ind, sigmahat_tempr, "Output/", "temperature.jpg")
p_t_set <- results[[1]]
p_t_set_plus <- results[[2]]
p_t_set_minus <- results[[3]]


jpeg(filename="Output/threegraphics.jpg")
par(mfrow = c(3,1)) #Setting the layout of the graphs
par(mar = c(0, 0.5, 0, 0.5)) #Margins for each plot
par(oma = c(0, 2, 2, 0)) #Outer margins


# Plotting the real data
plot(grid_points, yearly_tempr_normalised, type = "l") 


#This part plots kernel smoothers for different bandwidths (all on one plot).
smoothed_curve_1 <-mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr_normalised, grid_points, 0.01))
plot(grid_points, smoothed_curve_1, type = "l") 

smoothed_curve_2 <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr_normalised, grid_points, 0.05))
lines(grid_points, smoothed_curve_2, type = "l", col = "red")

smoothed_curve_3 <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr_normalised, grid_points, 0.1))
lines(grid_points, smoothed_curve_3, type = "l", col = "green")

smoothed_curve_4 <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr_normalised, grid_points, 0.15))
lines(grid_points, smoothed_curve_4, type = "l", col = "blue")

smoothed_curve_5 <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr_normalised, grid_points, 0.2))
lines(grid_points, smoothed_curve_5, type = "l", col = "yellow")

legend(-1/T_tempr, 1.3, legend=c("0.01","0.05", "0.10", "0.15", "0.20"), col=c("black", "red","green", "blue", "yellow"), lty = 1, ncol=1)


#Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
ymaxlim = max(p_t_set$values)
yminlim = min(p_t_set$values)
plot(NA, xlim=c(0,1), ylim = c(yminlim - 1, ymaxlim + 1), xlab="x", ylab="y")
segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set[['endpoint']], p_t_set[['values']])
#legend(1/T_tempr, 1, legend = c("Minimal (positive) intervals"), lty = 1, ncol = 1)


title(main = "Plots for normalised yearly temperature data for England",  outer = TRUE)
dev.off()