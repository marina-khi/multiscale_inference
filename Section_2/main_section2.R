library(xtable)
library(stargazer)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("functions.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
dyn.load("C_code/psihat_statistic.dll")
dyn.load("C_code/psihat_statistic_ll.dll")
source("C_code/psihat_statistic.R")
source("simulations_based_on_data.R")
source("data_analysis.R")


###############################
#Defining necessary parameters#
###############################

N <- 10000 #Number of replications for calculating the size and the power of the test
different_T     <- c(250, 350, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power
alpha <- 0.05 #alpha for calculating quantiles

h         <- c(0.05, 0.1, 0.15, 0.2) #Different bandwidth for plotting. Number must be <=6 in order for the plot to be readable

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported
test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 

#########################################################
#Loading the real data for yearly temperature in England#
#########################################################

temperature             <- read.table("data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr            <- temperature[temperature$YEAR > -99, 'YEAR']
yearly_tempr_normalised <- yearly_tempr - mean(yearly_tempr) #Normalization of the data

T_tempr <- length(yearly_tempr)

###################################################################
#Calculating smoothed curve for the data using Epanechnikov kernel#
###################################################################
grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating
end_point <- c()
plot(grid_points, yearly_tempr, ylim = c(-1.5, 1.5), type = 'l')

for (i in 1:length(h)){
  smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr_normalised, grid_points, h[i]))
  end_point <- c(end_point, smoothed_curve[T_tempr])
  cat("End point:", smoothed_curve[T_tempr], "\n")
}

pdf("../Plots/temperature_data.pdf", width=10, height=3, paper="special")
par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
data <- ts(yearly_tempr, start=1659, end=2017, frequency=1)
plot(data, ylab="", xlab = "", yaxp  = c(7, 11, 4), xaxp = c(1675, 2025, 7), type = 'l', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
dev.off()

############################################
#Calculating the power and size of the test#
############################################
simulations_based_on_data(N, different_T, different_alpha, yearly_tempr, test_problem, kernel_method)


###############
#Data analysis#
###############
data_analysis(alpha, yearly_tempr, "constant", kernel_method)
