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



# Parameters
Nsim       <- 5000                  # number of simulation runs for size/power calculations
sigma_eta  <- 1                     # standard deviation of the innovation term in the AR model
SimRuns    <- 5000                  # number of simulation runs to produce critical values


##########################################
#Calculating size for small parameters a1#
##########################################

different_T     <- c(250, 500, 1000)
different_a1    <- c(-0.5, -0.25, 0.25, 0.5)
different_alpha <- c(0.01, 0.05, 0.1)

number_of_cols            <- length(different_a1) * (length(different_alpha) + 1)
size_matrix_ms            <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix_ms)  <- different_T

#estimated lrv
k <- 1
for (a1 in different_a1){
  size <- calculating_size(a1, different_T, different_alpha, sigma_eta, Nsim = Nsim, kappa =0.1, SimRuns =SimRuns, type_of_sigma = 'estimated', q_ = 25, remove.small.ess = 'false')[[1]]
  size_matrix_ms[, (k * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(k * (length(different_alpha) + 1))] <- size
  k <- k + 1
}

print.xtable(xtable(size_matrix_ms, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
  type="latex", file=paste0("plots/size_estimated_sigma.tex"), include.colnames = FALSE)


#true lrv
size_matrix_ms            <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix_ms)  <- different_T

k <- 1
for (a1 in different_a1){
  size <- calculating_size(a1, different_T, different_alpha, sigma_eta, Nsim = Nsim, kappa =0.1, SimRuns =SimRuns, type_of_sigma = 'true', q_ = 25, remove.small.ess = 'false')[[1]]
  size_matrix_ms[, (k * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(k * (length(different_alpha) + 1))] <- size
  k <- k + 1
}

print.xtable(xtable(size_matrix_ms, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/size_true_sigma.tex"), include.colnames = FALSE)



######################################
#Calculating size for for a_1 = +-0.9#
######################################

different_T     <- c(250, 500, 1000, 2000, 3000, 4000, 5000)
different_a1    <- c(-0.9, 0.9)

number_of_cols              <- length(different_a1) * (length(different_T) + 1)
size_ms_est_sigma           <- matrix(NA, nrow = length(different_alpha), ncol = number_of_cols)
rownames(size_ms_est_sigma) <- different_alpha

#estimated lrv
k <- 1
for (a1 in different_a1){
  result <- calculating_size(a1, different_T, different_alpha, sigma_eta, Nsim = Nsim, kappa =0.1, SimRuns =SimRuns, type_of_sigma = 'estimated', q_ = 50, remove.small.ess = 'false')[[1]]
  size_ms_est_sigma[, (k * (length(different_T) + 1) - (length(different_T) - 1)):(k * (length(different_T) + 1))] <- t(result)
  k <- k + 1
}

print.xtable(xtable(size_ms_est_sigma, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/size_persistent_estimated_sigma.tex"), include.colnames = FALSE)


#true lrv
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

for (a1 in different_a1){
  result <- calculating_size_rowwise(a1, T, alpha, sigma_eta, Nsim = Nsim, SimRuns =SimRuns)
  h.grid <- result$h.grid
  
  pdffilename <- paste0("plots/pcp_rowwise_T_", T, "_a1_", a1*100, "_new_look.pdf")
  pdf(pdffilename, width=5.5, height=4.16, paper="special")
  
  par(mar = c(4, 4, 1.5, 0)) #Margins for each plot
  par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins 
  
  plot(x = h.grid, y = result$size_ms*100, ylim=c(0, 17), yaxp = c(0, 15, 3), type="l", lty=1, xaxt='n', ylab = "Rowwise % sig.", xlab='bandwidth')
  points(x = h.grid, y = result$size_ms*100, pch=19, cex = 0.6)
  
  lines(x = h.grid, y = result$size_uncor*100, lwd=1.5)#lty="dashed", 
  points(x = h.grid, y = result$size_uncor*100, pch=15, cex = 0.6)
  
  lines(x = h.grid, y = result$size_rows*100, lwd=1.5)#lty="dotted", 
  points(x = h.grid, y = result$size_rows*100, pch=17, cex = 0.6)

  lines(x = h.grid, y = result$size_SiZer*100, lwd=1.5)#lty="dotdash", 
  points(x = h.grid, y = result$size_SiZer*100, pch=8, cex = 0.6)
  
  abline(h = alpha*100, lty = 'dashed')
  
  axis(1, at=h.grid)
  legend('topleft', cex = 0.65, bty = "n", legend = c(expression(italic(T)[MS]), expression(italic(T)[UC]), expression(italic(T)[RW]), expression(italic(T)[SiZer])),
         #lty=c("solid","dashed","dotted", 'dotdash'),
         pch = c(19, 15, 17, 8), lwd=1.5, y.intersp=1.25)
  dev.off()
}