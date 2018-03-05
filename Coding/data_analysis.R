source("functions.R")
source("estimating_sigma.R")
dyn.load("psihat_statistic.dll")
source("psihat_statistic.R")

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

temperature <- read.table("cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']

T_tempr <- length(yearly_tempr)

#Tuning parameters
L1 <- floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))

sigmahat_tempr2 <- estimating_sigma_for_AR1(yearly_tempr, L1, L2)


#####################################
#Calculating gaussian quantile for T#
#####################################

g_t_set_tempr <- creating_g_set(T_tempr)
gaussian_quantile <- calculating_gaussian_quantile(T_tempr, g_t_set_tempr, kernel_ind, sqrt(sigmahat_tempr2), alpha)


#########################################
#Calculating the statistic for real data#
#########################################

results <- plotting_all_rejected_intervals(yearly_tempr, g_t_set_tempr, gaussian_quantile, kernel_ind, sqrt(sigmahat_tempr2), "Output/", "temperature.jpg")
p_t_set <- results[[1]]
p_t_set_plus <- results[[2]]
p_t_set_minus <- results[[3]]