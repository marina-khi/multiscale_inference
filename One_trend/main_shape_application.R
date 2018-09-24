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


#########################################################
#Loading the real data for yearly temperature in England#
#########################################################

temperature  <- read.table("Shape/data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
T_tempr      <- length(yearly_tempr)


#############################
#Checking the order of AR(p)#
#############################

L1 <- 30
L2 <- 40

FPE <- c()
AIC <- c()
SIC <- c()
HQ <- c()

different_orders <- (1:9)

for (order in different_orders){
  K1 <- order + 1
  K2 <- 10
  cat(order, "\n")
  a_hat_method1         <- AR_coefficients(yearly_tempr, L1, L2, rep(0,L2), order)
  sigma_eta_hat_method1 <- calculating_sigma_eta(yearly_tempr, a_hat_method1, order)
  sigma_hat_method1     <- sqrt(sigma_eta_hat_method1^2 / (1 - sum(a_hat_method1))^2) 

  corrections_value     <- corrections(a_hat_method1, sigma_eta_hat_method1, K2+1)
  a_hat_method2         <- AR_coefficients(yearly_tempr, K1, K2, corrections_value, order)
  sigma_eta_hat_method2 <- calculating_sigma_eta(yearly_tempr, a_hat_method2, order)
  sigma_hat_method2     <- sqrt(sigma_eta_hat_method2^2 / (1 - sum(a_hat_method2))^2) 
  FPE_p <- (sigma_eta_hat_method2^2 * (T_tempr + order)) / (T_tempr - order)
  FPE <- c(FPE, FPE_p)
  AIC <- c(AIC, T_tempr * log(sigma_eta_hat_method2^2) + 2 * order)
  SIC <- c(SIC, log(sigma_eta_hat_method2^2) + order * log(T_tempr) / T_tempr)
  HQ <- c(HQ, log(sigma_eta_hat_method2^2) + 2 * order * log(log(T_tempr)) / T_tempr)
}

###########################
#Setting tuning parameters#
###########################

#Tuning parameters
p    <- 2
#different_L1 <- c(10)
#tuning_parameter_grid <- (5:25)
K1 <- p+1
K2 <- 10

###################################################################
#Calculating smoothed curve for the data using Epanechnikov kernel#
###################################################################

# grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating
# 
# for (i in 1:length(h)){
#   if (kernel_method == "nw"){
#     smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(yearly_tempr, grid_points, h[i]))
#   } else if (kernel_method == "ll"){
#     smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(yearly_tempr, grid_points, h[i]))
#   } else {
#     print('Given method is currently not supported')
#   }
#   cat("End point:", smoothed_curve[T_tempr], "\n")
# }
#
# pdf("Paper/Plots/temperature_data.pdf", width=10, height=3, paper="special")
# par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
# par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
# data <- ts(yearly_tempr, start=1659, end=2017, frequency=1)
# plot(data, ylab="", xlab = "", yaxp  = c(7, 11, 4), xaxp = c(1675, 2025, 7), type = 'l', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
# dev.off()


###############
#Data analysis#
###############

#for (L1 in different_L1){
#  for (q in tuning_parameter_grid){
#    L2 <- L1 + q
    result    <- estimating_sigma_for_AR1(yearly_tempr, L1, L2)
    a_hat     <- result[[2]] #Estimation of the AR coefficient
    sigma_eta <- result[[3]] #Estimation of the sqrt of the variance of the innovation 


    a_hat_method1         <- AR_coefficients(yearly_tempr, L1, L2, rep(0,L2), p)
    sigma_eta_hat_method1 <- calculating_sigma_eta(yearly_tempr, a_hat_method1, p)
    sigma_hat_method1     <- sqrt(sigma_eta_hat_method1^2 / (1 - sum(a_hat_method1))^2) 

    corrections_value     <- corrections(a_hat_method1, sigma_eta_hat_method1, K2+1)
    a_hat_method2         <- AR_coefficients(yearly_tempr, K1, K2, corrections_value, p)
    sigma_eta_hat_method2 <- calculating_sigma_eta(yearly_tempr, a_hat_method2, p)
    sigma_hat_method2     <- sqrt(sigma_eta_hat_method2^2 / (1 - sum(a_hat_method2))^2) 

    cat("L1 = ", L1, ", L2 = ", L2, ", a_hall = ", a_hat, ", a_method1 = ", a_hat_method1, ", a_method2 = ", a_hat_method2, "\n")
    data_analysis(alpha, yearly_tempr, test_problem, kernel_method, sigma_hat_method2, L1, L2)
#  }
#}