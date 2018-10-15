library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/estimating_sigma_new.R")
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

L1 <- 20:35
L2 <- 20:45
criterion_matrix <- expand.grid(L1 = L1, L2 = L2)
criterion_matrix <- subset(criterion_matrix, L2 >= L1)

criterion_matrix$FPE <- numeric(length = nrow(criterion_matrix))
criterion_matrix$AIC <- numeric(length = nrow(criterion_matrix))
criterion_matrix$SIC <- numeric(length = nrow(criterion_matrix))
criterion_matrix$HQ  <- numeric(length = nrow(criterion_matrix))


for (i in 1:nrow(criterion_matrix)){
  FPE <- c()
  AIC <- c()
  SIC <- c()
  HQ <- c()
    
  different_orders <- (1:9)
  
  for (order in different_orders){
    K1 <- order + 1
    K2 <- 10
    sigma_eta_hat_method2 <- estimating_variance_new(yearly_tempr, criterion_matrix$L1[[i]], criterion_matrix$L2[[i]], order, K1, K2)[[3]]
  
    FPE <- c(FPE, (sigma_eta_hat_method2^2 * (T_tempr + order)) / (T_tempr - order))
    AIC <- c(AIC, T_tempr * log(sigma_eta_hat_method2^2) + 2 * order)
    SIC <- c(SIC, log(sigma_eta_hat_method2^2) + order * log(T_tempr) / T_tempr)
    HQ <- c(HQ, log(sigma_eta_hat_method2^2) + 2 * order * log(log(T_tempr)) / T_tempr)
  }
  criterion_matrix$FPE[[i]] <- which.min(FPE)
  criterion_matrix$AIC[[i]] <- which.min(AIC)
  criterion_matrix$SIC[[i]] <- which.min(SIC)
  criterion_matrix$HQ[[i]]  <- which.min(HQ)
}



###########################
#Setting tuning parameters#
###########################

#Tuning parameters
p  <- 2
L1 <- 25
L2 <- 25
K1 <- p + 1
K2 <- 10

different_L1          <- (11:40)
tuning_parameter_grid <- (5:25)


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
    sigma_hat <- estimating_variance_new(yearly_tempr, L1, L2, p, K1, K2)[[1]]
    data_analysis(alpha, yearly_tempr, test_problem, kernel_method, sigma_hat, L1, L2)
#  }
#}