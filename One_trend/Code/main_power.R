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
sigma_eta  <- 1                     # standard deviation of the innovation term in the AR model
SimRuns    <- 1000                  # number of simulation runs to produce critical values
sim.design <- 'bump'

##########################
#Calculating global power#
##########################

T            <- 500
different_a1 <- c(-0.5, 0.5)
alpha        <- 0.05

power_matrix           <- matrix(NA, nrow = 2, ncol = 5*length(different_a1))
rownames(power_matrix) <- c('Power', 'Spurious power')

k <- 1
for (a1 in different_a1){
  if (a1 < 0) {height = 1.2}
  else {height = 3.8}
  power_overall <- calculating_power(a1, T, alpha, sigma_eta, Nsim = Nsim, SimRuns =SimRuns, 
                             type_of_sigma = 'true', remove.small.ess = 'true', 
                             sim.design = 'bump', bump.height = height, region = 'increase',
                             type_of_power = 'global')
  power_spurious <- calculating_power(a1, T, alpha, sigma_eta, Nsim = Nsim, SimRuns =SimRuns, 
                             type_of_sigma = 'true', bump.height = height, remove.small.ess = 'true', 
                             sim.design = 'bump', region = 'increase',
                             type_of_power = 'spurious')
  
  power_matrix[1, ((k - 1) * 5 + 2):(k * 5)] <- power_overall
  power_matrix[2, ((k - 1) * 5 + 2):(k * 5)] <- power_spurious
  k <- k + 1
}

print.xtable(xtable(power_matrix, digits = c(3), align = paste(replicate(5*length(different_a1) + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/power_table.tex"), include.colnames = FALSE)

############################
#Calculating rowwise power#
############################

for (a1 in different_a1){
  if (a1 < 0) {height = 1.2}
  else {height = 3.8}
  
  result <- calculating_power_rowwise(a1, T, alpha, sigma_eta, Nsim = Nsim, SimRuns =SimRuns,
                                      type_of_sigma = 'true', remove.small.ess = 'true',
                                      sim.design = 'bump', bump.height = height,
                                      region = 'increase', type_of_power = 'global')

  result_spurious <- calculating_power_rowwise(a1, T, alpha, sigma_eta, Nsim = Nsim, SimRuns =SimRuns,
                                      type_of_sigma = 'true', remove.small.ess = 'true',
                                      sim.design = 'bump', bump.height = height,
                                      region = 'increase', type_of_power = 'spurious')
  
  h.grid <- result$h.grid
  
  pdffilename <- paste0("plots/pcp_power_T_", T, "_a1_", a1*100, ".pdf")
  pdf(pdffilename, width = 13, height = 3, paper = 'special')

  par(mfrow = c(1, 3), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
  
  par(mar = c(4, 4, 1.5, 0)) #Margins for each plot
  par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins 
  
  data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = sim.design, slope.fac = height)
  data           <- data.simulated$data
  trend          <- data.simulated$trend
  
  plot(seq(1/T, 1, by = 1/T), data, ylab = "", xlab = "")
  lines(seq(1/T, 1, by = 1/T), trend)
  
  plot(x = h.grid, y = result$power_ms*100, ylim=c(0, 65), yaxp = c(0, 60, 4), type="l", lty=1, xaxt='n', ylab = "Rowwise % sig.", xlab='bandwidth')
  points(x = h.grid, y = result$power_ms*100, pch=19, cex = 0.6)
  
  lines(x = h.grid, y = result$power_uncor*100, lwd=1.5)#lty="dashed", 
  points(x = h.grid, y = result$power_uncor*100, pch=15, cex = 0.6)
  
  lines(x = h.grid, y = result$power_rows*100, lwd=1.5)#lty="dotted", 
  points(x = h.grid, y = result$power_rows*100, pch=17, cex = 0.6)
  
  lines(x = h.grid, y = result$power_SiZer*100, lwd=1.5)#lty="dotdash", 
  points(x = h.grid, y = result$power_SiZer*100, pch=8, cex = 0.6)
  
  axis(1, at=h.grid)
  legend('topleft', cex = 0.65, bty = "n", legend = c(expression(italic(T)[MS]), expression(italic(T)[UC]), expression(italic(T)[RW]), expression(italic(T)[SiZer])),
         #lty=c("solid","dashed","dotted", 'dotdash'),
         pch = c(19, 15, 17, 8), lwd=1.5, y.intersp=1.25)
  
  plot(x = h.grid, y = result_spurious$power_ms*100, ylim=c(0, 8), type="l", lty=1, xaxt='n', ylab = "Rowwise % sig.", xlab='bandwidth')
  points(x = h.grid, y = result_spurious$power_ms*100, pch=19, cex = 0.6)
  
  lines(x = h.grid, y = result_spurious$power_uncor*100, lwd=1.5)#lty="dashed", 
  points(x = h.grid, y = result_spurious$power_uncor*100, pch=15, cex = 0.6)
  
  lines(x = h.grid, y = result_spurious$power_rows*100, lwd=1.5)#lty="dotted", 
  points(x = h.grid, y = result_spurious$power_rows*100, pch=17, cex = 0.6)
  
  lines(x = h.grid, y = result_spurious$power_SiZer*100, lwd=1.5)#lty="dotdash", 
  points(x = h.grid, y = result_spurious$power_SiZer*100, pch=8, cex = 0.6)
  
  axis(1, at=h.grid)
  legend('topleft', cex = 0.65, bty = "n", legend = c(expression(italic(T)[MS]), expression(italic(T)[UC]), expression(italic(T)[RW]), expression(italic(T)[SiZer])),
         #lty=c("solid","dashed","dotted", 'dotdash'),
         pch = c(19, 15, 17, 8), lwd=1.5, y.intersp=1.25)
  
  dev.off()
}