library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/SiZer_simulations.R")

source("Shape/C_code/psihat_statistic.R")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")


###############################
#Defining necessary parameters#
###############################

N_min_intervals <- 100 #Number of replications for calculating the minimal intervals and producing plots
N_rep           <- 1000
sigma_eta       <- 1    #We keep this as a constant parameter

different_alpha   <- c(0.01, 0.05, 0.10) #Level of significance
different_heights <- c(32/15) #Different strength of the signal calculated as height * 15/16
different_widths  <- c(10) #Different support of the signal calculated by [0.5 - 1/width, 0.5 + 1/width]
different_T       <- c(250, 350, 500) #Different lengths of time series for which we compare SiZer and our method
different_a       <- c(-0.5, -0.25, 0.25, 0.5) #Different a_1 in AR(1) model


########################################
#Producing plots with minimal intervals#
########################################

# T_size <- 500
# 
# different_i <- seq(from = 1/T_size, to = 1, by = 1/T_size)
# different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
# 
# for (a_1 in different_a){
#   gamma = c()
#   for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
#     gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
#   }
# 
#   #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
#   true_var <- gamma[1] / T_size
#   for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}
# 
#   T_star   <- gamma[1]/true_var
# 
#   SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)
# 
#   #Gaussian statistic for our own method
#   gaussian_quantile <- calculating_gaussian_quantile_ll(T_size, SiZer_matrix, "comparison", kernel_ind = 2, alpha)
# 
#   sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))
# 
#   for (height in different_heights){
#     for (width in different_widths){
#       plotting_many_minimal_intervals(height, width, T_size, SiZer_matrix, N_min_intervals, kernel_ind = 2, sigmahat, gaussian_quantile, a_1, sigma_eta)
#     }
#   }
# }
# 

#############################################
#Calculating size and power for both methods#
#############################################
for (a_1 in different_a){
  if (a_1 > 0){
    slopes <- c(4.0, 4.5, 5.0, 5.5, 6.0)
  } else {
    slopes <- c(1.5, 2.0, 2.5)
  }
  size <- c()
  power35 <- c()
  power40 <- c()
  power45 <- c()
  power50 <- c()
  for (alpha in different_alpha){
    for (T_size in different_T){
      result <- SiZer_simulations(T_size, a_1, sigma_eta, alpha, N_rep, slopes)
      size <- c(size, result[[1]])
      power35 <- c(power35, result[[2]][1:2])
      power40 <- c(power40, result[[2]][3:4])
      power45 <- c(power45, result[[2]][5:6])
      power50 <- c(power50, result[[2]][7:8])
    }
  }
  PDFpartialPath = paste0("Paper/Plots/SiZer_comparison_a1_", a_1*100)
  creating_matrix_and_texing_for_SiZer(size, different_T, different_alpha, filename = paste0(PDFpartialPath, "_size.tex"))
  creating_matrix_and_texing_for_SiZer(power35, different_T, different_alpha, filename = paste0(PDFpartialPath, "_power35.tex"))
  creating_matrix_and_texing_for_SiZer(power40, different_T, different_alpha, filename = paste0(PDFpartialPath, "_power40.tex"))
  creating_matrix_and_texing_for_SiZer(power45, different_T, different_alpha, filename = paste0(PDFpartialPath, "_power45.tex"))
  creating_matrix_and_texing_for_SiZer(power50, different_T, different_alpha, filename = paste0(PDFpartialPath, "_power50.tex"))
}