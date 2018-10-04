library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/estimating_sigma_new.R")

dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/simulations_based_on_data.R")


###############################
#Defining necessary parameters#
###############################
N_rep            <- 1000 #Number of replications for calculating the size and the power of the test
different_T      <- c(250, 350, 500) #Different lengths of time series for which we calculate size and power
different_alpha  <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power
different_a1     <- c(-0.5, -0.25, 0.25, 0.5)
different_slopes <- c(1.25, 1.875, 2.5)
sigma_eta        <- 1

kernel_method <- "ll" #Only "nw" (Nadaraya-Watson) and "ll" (local linear) methods are currently supported
test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 

L1 <- 20

PDFname <- "Paper/Plots/finite_sample_properties_"

set.seed(1) #For reproducibility

################################################
#Loading data and estimating parameters from it#
################################################
temperature             <- read.table("Shape/data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr            <- temperature[temperature$YEAR > -99, 'YEAR']

results            <- estimating_variance_new(yearly_tempr, L1, L1, order = 2, K1 = 2 + 1, K2 = 10)
sigma_hat_data     <- results[[1]]
a_hat_data         <- results[[2]]
sigma_eta_hat_data <- results[[3]]


######################################################
#Calculating the power and size of the test for AR(1)#
######################################################
matrix_size  <- matrix(NA, nrow = length(different_T), ncol = (length(different_alpha) + 1) * (length(different_a1) + 1), byrow = TRUE)
matrix_power <- matrix(NA, nrow = length(different_slopes) * length(different_T), ncol = (length(different_alpha) + 1) * (length(different_a1) + 1), byrow = TRUE)

rownames(matrix_size)  <- different_T
rownames(matrix_power) <- replicate(length(different_slopes), different_T)

power <- c()

i <- 1
for (a_hat in different_a1){
 result      <- simulations_general(N_rep, different_T, different_alpha, different_slopes, a_hat, sigma_eta, order = 1, test_problem, kernel_method, L1, K1 = 1 + 1, K2 = 10)
 matrix_size[, (i * 4 - 2):(i * 4)]  <- result[[1]]
 matrix_power[, (i * 4 - 2):(i * 4)] <- result[[2]] 
 i <- i + 1
}


################################################################################
#Calculating the power and size of the test for one specification based on data#
################################################################################
result_data <- simulations_general(N_rep, different_T, different_alpha, different_slopes, a_hat_data, sigma_eta_hat_data, order = 2, test_problem, kernel_method, L1, K1 = 2 + 1, K2 = 10)
matrix_size[, ((length(different_alpha) + 1) * (length(different_a1) + 1) - 2):((length(different_alpha) + 1) * (length(different_a1) + 1))]  <- result_data[[1]]
matrix_power[, ((length(different_alpha) + 1) * (length(different_a1) + 1) - 2):((length(different_alpha) + 1) * (length(different_a1) + 1))] <- result_data[[2]] 


##############################
#Writing the results in files#
##############################
print.xtable(xtable(matrix_size, digits = c(3), align = paste(replicate((length(different_alpha) + 1) * (length(different_a1) + 1) + 1, "c"), collapse = "")),
             type="latex", file=paste0(PDFname, "size.tex"), include.colnames = FALSE)

j <- 1
for (slope in different_slopes){
  print.xtable(xtable(matrix_power[(3 * j - 2):(3 * j),], digits = c(3), align = paste(replicate((length(different_alpha) + 1) * (length(different_a1) + 1) + 1, "c"), collapse = "")),
               type="latex", file=paste0(PDFname, "power_", slope, ".tex"), include.colnames = FALSE)
  j <- j + 1
}