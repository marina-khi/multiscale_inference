rm(list=ls())

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/grid_construction.r")
dyn.load("functions/kernel_weights.dll")
source("functions/kernel_weights.r")
source("functions/multiscale_statistics.r")
source("functions/multiscale_quantiles.r")
source("functions/multiscale_testing.r")

source("functions/long_run_variance.r")
source("functions/sim.r")
source("functions/calculating_size.r")
source("functions/SiZer_functions.r")
dyn.load("functions/SiZer_functions.dll")

##################################################################
#Calculating size for the estimated sigma for small parameters a1#
##################################################################

# Parameters
different_T     <- c(250, 500, 1000)
different_a1    <- c(-0.5, -0.25, 0.25, 0.5)
different_alpha <- c(0.01, 0.05, 0.1)

Nsim       <- 1000                  # number of simulation runs for size/power calculations
sigma_eta  <- 1                     # standard deviation of the innovation term in the AR model
SimRuns    <- 1000                  # number of simulation runs to produce critical values
kappa      <- 0.1                   # parameter to determine order statistic for the version
                                    # of the multiscale statistic from Section ?? 

number_of_cols            <- length(different_a1) * (length(different_alpha) + 1)
size_matrix_ms            <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix_ms)  <- different_T

k <- 1
for (a1 in different_a1){
  size <- calculating_size(a1, different_T, different_alpha, sigma_eta, Nsim = Nsim, kappa =kappa, SimRuns =SimRuns, type_of_sigma = 'estimated', q_ = 25, remove.small.ess = 'false')
  size_matrix_ms[, (k * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(k * (length(different_alpha) + 1))] <- size
  k <- k + 1
}

print.xtable(xtable(size_matrix_ms, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
  type="latex", file=paste0("plots/size_estimated_sigma.tex"), include.colnames = FALSE)


#############################################################
#Calculating size for the true sigma for small parameters a1#
#############################################################

size_matrix_ms            <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix_ms)  <- different_T

k <- 1
for (a1 in different_a1){
  size <- calculating_size(a1, different_T, different_alpha, sigma_eta, Nsim = Nsim, kappa =kappa, SimRuns =SimRuns, type_of_sigma = 'true', q_ = 25, remove.small.ess = 'false')[[1]]
  size_matrix_ms[, (k * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(k * (length(different_alpha) + 1))] <- size
  k <- k + 1
}

print.xtable(xtable(size_matrix_ms, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/size_true_sigma.tex"), include.colnames = FALSE)





#Calculating size for the estimated sigma for a_1 = +-0.9
different_T     <- c(250, 500, 1000, 2000, 3000, 4000, 5000)
different_a1    <- c(-0.9, 0.9)

number_of_cols              <- length(different_a1) * (length(different_T) + 1)
size_ms_est_sigma           <- matrix(NA, nrow = length(different_alpha), ncol = number_of_cols)
rownames(size_ms_est_sigma) <- different_alpha

k <- 1
for (a1 in different_a1){
  result <- calculating_size(a1, different_T, different_alpha, sigma_eta, Nsim = Nsim, kappa =kappa, SimRuns =SimRuns, type_of_sigma = 'estimated', q_ = 50, remove.small.ess = 'false')[[1]]
  size_ms_est_sigma[, (k * (length(different_T) + 1) - (length(different_T) - 1)):(k * (length(different_T) + 1))] <- t(result)
  k <- k + 1
}

print.xtable(xtable(size_ms_est_sigma, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/size_persistent_estimated_sigma.tex"), include.colnames = FALSE)


#Calculating size for the true sigma for for a_1 = +-0.9
size_ms_true_sigma           <- matrix(NA, nrow = length(different_alpha), ncol = number_of_cols)
rownames(size_ms_true_sigma) <- different_alpha

k <- 1
for (a1 in different_a1){
  result <- calculating_size(a1, different_T, different_alpha, sigma_eta, Nsim = Nsim, kappa =0.1, SimRuns =SimRuns, type_of_sigma = 'true', q_ = 50, remove.small.ess = 'false')[[1]]
  size_ms_true_sigma[, (k * (length(different_T) + 1) - (length(different_T) - 1)):(k * (length(different_T) + 1))] <- t(result)
  k <- k + 1
}

print.xtable(xtable(size_ms_true_sigma, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/size_persistent_true_sigma.tex"), include.colnames = FALSE)


#########################
#Calculating global size#
#########################

different_T     <- c(250, 500, 1000)
different_a1    <- c(-0.5, 0.5)
different_alpha <- c(0.05)

Nsim       <- 1000                  # number of simulation runs for size/power calculations
sigma_eta  <- 1                     # standard deviation of the innovation term in the AR model
SimRuns    <- 1000                  # number of simulation runs to produce critical values

number_of_cols         <- length(different_a1) * 5
size_matrix            <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix)  <- different_T

k <- 1
for (a1 in different_a1){
  size <- calculating_size(a1, different_T, different_alpha, sigma_eta, Nsim = Nsim, kappa =0.1, SimRuns =SimRuns, type_of_sigma = 'true', q_ = 25, remove.small.ess = 'true')
  size_matrix[, (k-1)*5 + 2] <- size[[1]]
  size_matrix[, (k-1)*5 + 3] <- size[[2]]
  size_matrix[, (k-1)*5 + 4] <- size[[3]]
  size.SiZer <- calculating_size_for_SiZer(a1, different_T, different_alpha, sigma_eta, Nsim = 1000, SimRuns =1000)
  size_matrix[, (k-1)*5 + 5] <- size.SiZer
  k <- k + 1
}

print.xtable(xtable(size_matrix, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/comparing_size.tex"), include.colnames = FALSE)

#################################################################
#Calculating rowwise size and creating parallel coordinate plots#
#################################################################

T            <- 500
different_a1 <- c(-0.5, 0.5)
alpha        <- 0.05

Nsim       <- 1000                  # number of simulation runs for size/power calculations
sigma_eta  <- 1                     # standard deviation of the innovation term in the AR model
SimRuns    <- 1000                  # number of simulation runs to produce critical values

for (a1 in different_a1){
  result <- calculating_size_rowwise(a1, T, alpha, sigma_eta, Nsim = Nsim, kappa =0.1, SimRuns =SimRuns)
  h.grid <- result$h.grid
  
  pdffilename <- paste0("plots/pcp_rowwise_T_", T, "_a1_", a1*100, ".pdf")
  pdf(pdffilename, height = 7, width = 10, paper = 'special')

  par(mfrow = c(1,1), cex = 1.0, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
  par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins
  
  plot(h.grid, result$size_ms*100, ylim = c(0, 15), xlab = 'bandwidth',
       ylab = "Rowwise % sig.", type = 'l')
  points(h.grid, result$size_uncor*100, type = 'l', lty = 2)
  points(h.grid, result$size_rows*100, type = 'l', lty = 3)
  points(h.grid, result$size_SiZer*100, type = 'l', lty = 4)
  
  legend(0.03, 15, legend=c("MS method", "Uncorrected version", "Rowwise version", "SiZer method"), lty=1:4, cex=1.1)
  dev.off()
}