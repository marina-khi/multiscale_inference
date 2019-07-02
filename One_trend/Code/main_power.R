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
source("functions/calculating_power.r")
source("functions/SiZer_functions.r")
dyn.load("functions/SiZer_functions.dll")



# Parameters
Nsim       <- 1000                  # number of simulation runs for size/power calculations
sigma_eta  <- 0.6                     # standard deviation of the innovation term in the AR model
SimRuns    <- 1000                  # number of simulation runs to produce critical values
sim.design <- 'bump'

######################################
#Calculating global and rowwise power#
######################################

T            <- 500
different_a1 <- c(-0.5, 0.5)
alpha        <- 0.05

power_matrix           <- matrix(NA, nrow = 2, ncol = 5*length(different_a1))
rownames(power_matrix) <- c('Power', 'Spurious power')

k <- 1
for (a1 in different_a1){
  power <- calculating_power(a1, T, alpha, sigma_eta, Nsim = Nsim, SimRuns =SimRuns, type_of_sigma = 'true', remove.small.ess = 'true', sim.design = sim.design)
  power_matrix[1, ((k - 1) * 5 + 2):(k * 5)] <- power$overall
  power_matrix[2, ((k - 1) * 5 + 2):(k * 5)] <- power$spurious
  k <- k + 1
}

print.xtable(xtable(power_matrix, digits = c(3), align = paste(replicate(5*length(different_a1) + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/power_table.tex"), include.colnames = FALSE)