source("Shape/functions.R")
#source("Shape/C_code/estimating_sigma.R")
#dyn.load("Shape/C_code/estimating_sigma.dll")
dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/SiZer_simulations.R")

###############################
#Defining necessary parameters#
###############################

N_rep           <- 100 #Number of replications for calculating the size and the power of the test

sigma_eta   <- 0.5
alpha       <- 0.05

kernel_ind  <- 2

T_size      <- 500
different_i <- seq(from = 1/T_size, to = 1, by = 1/T_size)
different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)

different_heights <- c(0.5, 1, 1.5)
different_widths  <- c(8, 10, 12)


############################################################
#Calculating everything that does not depend on Y for SiZer#
############################################################
for (a_1 in c(-0.5, -0.25, 0.25, 0.5)){
  gamma = c()
  for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
    gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
  }

  #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
  true_var <- gamma[1] / T_size
  for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}

  T_star   <- gamma[1]/true_var

  SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)

  #################################################################
  #Calculating everything that does not depend on Y for our method#
  #################################################################

  #Gaussian statistic for our own method
  gaussian_statistic_distribution <- replicate(1000, {
    z = rnorm(T_size, 0, 1)
    psistar_statistic_ll(z, SiZer_matrix, kernel_ind, 1)
  })
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)

  sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))

  for (height in different_heights){
    for (width in different_widths){
      plotting_many_minimal_intervals(height, width, T_size, SiZer_matrix, N_rep, kernel_ind, sigmahat, gaussian_quantile, a_1, sigma_eta)
      cat("a_1 = ", a_1, ", height coeffiecient = ", height, ", width coefficient = ", width, "\n")
    }
  }
}