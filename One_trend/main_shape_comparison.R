#library(xtable)
#options(xtable.floating = FALSE)
#options(xtable.timestamp = "")
source("Shape/functions.R")

source("Shape/C_code/psihat_statistic.R")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
#source("Shape/SiZer_simulations.R")


###############################
#Defining necessary parameters#
###############################

N_min_intervals <- 100 #Number of replications for calculating the minimal intervals and producing plots

sigma_eta   <- 1    #We keep this as a constant parameter
alpha       <- 0.05 #Level of significance

different_heights <- c(1, 2, 3) #Different strength of the signal calculated as height * 15/16
different_widths  <- c(8, 10, 12) #Different support of the signal calculated by [0.5 - 1/width, 0.5 + 1/width]
different_T       <- c(200, 250, 500) #Different lengths of time series for which we compare SiZer and our method
different_a       <- c(-0.5, -0.25, 0.25, 0.5) #Different a_1 in AR(1) model


########################################
#Producing plots with minimal intervals#
########################################

for (T_size in different_T){
  different_i <- seq(from = 1/T_size, to = 1, by = 1/T_size)
  different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
  for (a_1 in different_a){
    gamma = c()
    for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
      gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
    }
  
    #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
    true_var <- gamma[1] / T_size
    for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}
  
    T_star   <- gamma[1]/true_var
  
    SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)
  
    #Gaussian statistic for our own method
    gaussian_statistic_distribution <- replicate(1000, {
      z = rnorm(T_size, 0, 1)
      psistar_statistic_ll(z, SiZer_matrix, kernel_ind = 2, 1)
    })
    gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)
  
    sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))
  
    for (height in different_heights){
      for (width in different_widths){
        plotting_many_minimal_intervals(height, width, T_size, SiZer_matrix, N_min_intervals, kernel_ind = 2, sigmahat, gaussian_quantile, a_1, sigma_eta)
#        cat("a_1 = ", a_1, ", height coeffiecient = ", height, ", width coefficient = ", width, "\n")
      }
    }
  }
}

#for (T_size in different_T){
#  PDFPath = paste0("Paper/Plots/SiZer_comparison_", T_size, ".pdf")
#  SiZer_simulations(T_size, a_1, sigma_eta, alpha, N_rep, PDFpath)
#}