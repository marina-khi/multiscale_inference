library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions_for_referees.R")
source("Shape/estimating_sigma_new.R")
#dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_without_lambda.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/simulations_based_on_data.R")


###############################
#Defining necessary parameters#
###############################
N_rep            <- 1000 #Number of replications for calculating the size and the power of the test
different_T      <- c(250, 350, 500) #Different lengths of time series for which we calculate size and power
different_alpha  <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power
different_a1     <- c(-0.5, -0.25, 0.25, 0.5) #Different AR(1) parameters
different_slopes <- c(1.5, 2.0, 2.5) #Slopes for power simulations
sigma_eta        <- 1 #Sqrt(variance) for the innovation \eta_t

test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 

q <- 25 #Parameter q for the pilot estimator \widehat{a}_q

PDFname <- "JRSSB_submission/Plots/finite_sample_properties_" #Where to store the tables


######################################################
#Calculating the power and size of the test for AR(1)#
######################################################
matrix_size  <- matrix(NA, nrow = length(different_T), ncol = (length(different_alpha) + 1) * (length(different_a1) + 1), byrow = TRUE)
#matrix_power <- matrix(NA, nrow = length(different_slopes) * length(different_T), ncol = (length(different_alpha) + 1) * (length(different_a1) + 1), byrow = TRUE)

rownames(matrix_size)  <- different_T
#rownames(matrix_power) <- replicate(length(different_slopes), different_T)

#power <- c()

i <- 1
for (a_hat in different_a1){
 result <- simulations_general_without_lambda(N_rep, different_T, different_alpha, different_slopes, a_hat, sigma_eta, order = 1, test_problem, q, K2 = 10)
 
 matrix_size[, (i * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(i * (length(different_alpha) + 1))]  <- result[[1]]
 #matrix_power[, (i * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(i * (length(different_alpha) + 1))] <- result[[2]] 
 i <- i + 1
}


############################################################
#Loading temperature data and estimating parameters from it#
############################################################
temperature             <- read.table("Shape/data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr            <- temperature[temperature$YEAR > -99, 'YEAR']

results            <- estimating_variance_new(yearly_tempr, q, order = 2, K2 = 10)
sigma_hat_data     <- results[[1]]
a_hat_data         <- results[[2]]
sigma_eta_hat_data <- results[[3]]


################################################################################
#Calculating the power and size of the test for one specification based on data#
################################################################################
result_data <- simulations_general_without_lambda(N_rep, different_T, different_alpha, different_slopes, a_hat_data, sigma_eta_hat_data, order = 2, test_problem, q, K2 = 10)

matrix_size[, ((length(different_alpha) + 1) * (length(different_a1) + 1) - 2):((length(different_alpha) + 1) * (length(different_a1) + 1))]  <- result_data[[1]]
#matrix_power[, ((length(different_alpha) + 1) * (length(different_a1) + 1) - 2):((length(different_alpha) + 1) * (length(different_a1) + 1))] <- result_data[[2]] 


##############################
#Writing the results in files#
##############################
print.xtable(xtable(matrix_size, digits = c(3), align = paste(replicate((length(different_alpha) + 1) * (length(different_a1) + 1) + 1, "c"), collapse = "")),
             type="latex", file=paste0(PDFname, "size_without_lambda.tex"), include.colnames = FALSE)

#j <- 1
#for (slope in different_slopes){
#  print.xtable(xtable(matrix_power[(length(different_T) * j - (length(different_T) - 1)):(length(different_T) * j),], digits = c(3), align = paste(replicate((length(different_alpha) + 1) * (length(different_a1) + 1) + 1, "c"), collapse = "")),
#               type="latex", file=paste0(PDFname, "power_", slope, "_without_lambda.tex"), include.colnames = FALSE)
#  j <- j + 1
#}