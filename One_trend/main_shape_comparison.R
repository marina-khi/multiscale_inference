library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
#source("Shape/functions.R")
source("Shape/functions_new_SiZer.R")

source("Shape/SiZer_simulations.R")
source("Shape/C_code/psihat_statistic.R")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")


###############################
#Defining necessary parameters#
###############################
N_rep           <- 1000 #Number of simulations
sigma_eta       <- 1    #Standard deviation of the innovation term

different_alpha     <- c(0.01, 0.05, 0.10) #Level of significance
different_T         <- c(250, 350,500) #Different lengths of time series for which we compare SiZer and our method
different_a1        <- c(-0.25, 0.25) #Different a_1 in AR(1) model
slopes_for_negative <- c(1.0, 1.25, 1.5) #Slopes for power calculations for negative a_1
slopes_for_positive <- c(2.5, 2.75, 3.0)#Slopes for power calculations for positive a_1

PDFname <- "Paper/Plots/SiZer_comparison_" #Path for tables


########################################
#Producing plots with minimal intervals#
########################################
N_min_intervals   <- 100 #Number of replications for calculating the minimal intervals and producing plots
different_heights <- c(32/15) #Different strength of the signal calculated as height * 15/16
different_widths  <- c(10) #Different support of the signal calculated by [0.5 - 1/width, 0.5 + 1/width]

T_size <- 500 #Sample size for plotting minimal intervals
alpha  <- 0.05 #Confidence level

different_i <- seq(from = 5/T_size, to = 1, by = 5/T_size) #Grid G_T from Section 5.1
different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)

for (a_1 in different_a1){
  set.seed(1) #For reproducibility
  gamma = c()
  for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
    gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
  }

  #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
  true_var <- gamma[1] / T_size
  for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k + 1]}

  T_star   <- gamma[1]/true_var

  #THIS PART IS ONLY FOR FAST CALCULATIONS OF THE MATRIX. IF YOU CHANGE different_i OR different_h, YOU CAN'T USE IT!!!
  filename_aux = paste0("Shape/distribution/SiZer_matrix_T_", T_size, "_a_1_", a_1, "_sigma_eta_", sigma_eta, "_alpha_", alpha*100, ".RData")
  if(!file.exists(filename_aux)) {
    SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma, a_1, sigma_eta)
    save(SiZer_matrix, file = filename_aux)
  } else {
    load(filename_aux)
  }
  #USE THIS INSTEAD
  #SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma, a_1, sigma_eta)

  #Gaussian statistic for our own method
  gaussian_quantile <- calculating_gaussian_quantile_ll(T_size, SiZer_matrix, "comparison", kernel_ind = 2, alpha)
  
  #True sqrt of long-run variance
  sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))

  #Actually producing the  graphs with minimal intervals
  for (height in different_heights){
    for (width in different_widths){
      plotting_many_minimal_intervals(height, width, T_size, SiZer_matrix, N_min_intervals, kernel_ind = 2, sigmahat, gaussian_quantile, a_1, sigma_eta)
    }
  }
}


#############################################
#Calculating size and power for both methods#
#############################################
matrix_size  <- matrix(NA, nrow = length(different_T), ncol = (2 * length(different_alpha) + 1) * length(different_a1), byrow = TRUE)
rownames(matrix_size)  <- different_T

matrix_power  <- matrix(NA, nrow = length(slopes_for_negative) * length(different_T), ncol = (2 * length(different_alpha) + 1) * length(different_a1), byrow = TRUE)
rownames(matrix_power)  <- replicate(length(slopes_for_positive), different_T)

i <- 0
for (a_1 in different_a1){
  if (a_1 > 0){
    slopes <- slopes_for_positive
  } else {
    slopes <- slopes_for_negative
  }
  result_power <- SiZer_simulations_power(a_1, sigma_eta, N_rep, slopes, different_alpha, different_T)
  tmp_power    <- matrix(result_power, nrow = length(slopes_for_negative) * length(different_T), ncol = 2 * length(different_alpha), byrow = TRUE)
  matrix_power[, (i * (2 * length(different_alpha) + 1) + 2):((i + 1) * (2 * length(different_alpha) + 1))]  <- tmp_power

  result_size <- SiZer_simulations_size(a_1, sigma_eta, N_rep, different_alpha, different_T)
  tmp_size    <- matrix(result_size,nrow = length(different_T), ncol = 2 * length(different_alpha), byrow = TRUE)
  matrix_size[, (i * (2 * length(different_alpha) + 1) + 2):((i + 1) * (2 * length(different_alpha) + 1))]  <- tmp_size

  i <- i + 1
}


##################
#Producing tables#
##################
print.xtable(xtable(matrix_size, digits = c(3), align = paste(replicate((2 * length(different_alpha) + 1) * length(different_a1) + 1, "c"), collapse = "")),
            type="latex", file=paste0(PDFname, "_size_025.tex"), include.colnames = FALSE)

j <- 1
for (slope in slopes_for_negative){
  print.xtable(xtable(matrix_power[(length(different_T) * j - (length(different_T) - 1)):(length(different_T) * j),], digits = c(3), align = paste(replicate((2 * length(different_alpha) + 1) * length(different_a1) + 1, "c"), collapse = "")),
               type="latex", file=paste0(PDFname, "_power_", slope*100, "_025.tex"), include.colnames = FALSE)
  j <- j + 1
}