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
Nsim       <- 1000    # number of simulation runs for size/power calculations
sigma_eta  <- 1       # standard deviation of the innovation term in the AR model
SimRuns    <- 5000    # number of simulation runs to produce critical values

T            <- 1000         #Sample size
different_a1 <- c(-0.5, 0.5) #different a_1 parameters
alpha        <- 0.05         #Significance level
different_T  <- c(250, 500, 1000)

##########################
#Calculating global power#
##########################
sim.design <- 'bump'  # type of trend function m()
height.neg <- 0.85    # height of the bump in case of sim.design = 'bump'. This height is used to calculate power for negative a1     
height.pos <- 2.65    # height of the bump in case of sim.design = 'bump'. This height is used to calculate power for positive a1 

power_matrix           <- matrix(NA, nrow = 2 * length(different_T), ncol = 5*length(different_a1) + 1)
rownames(power_matrix) <- c(250, 250, 500, 500, 1000, 1000)
power_matrix[, 1] <- rep(c("Power", "Spurious power"), length(different_T)) 

k <- 1
for (a1 in different_a1){
  if (a1 < 0) {height = height.neg}
  else {height = height.pos}
  i <- 1
  for (T in different_T){
    power_overall  <- calculating_power(a1, T, alpha, sigma_eta, Nsim = Nsim, SimRuns =SimRuns,
                             type_of_sigma = 'true', remove.small.ess = 'true',
                             sim.design = 'bump', bump.height = height, region = 'increase',
                             type_of_power = 'global')
    power_spurious <- calculating_power(a1, T, alpha, sigma_eta, Nsim = Nsim, SimRuns =SimRuns,
                             type_of_sigma = 'true', remove.small.ess = 'true',
                             sim.design = 'bump', bump.height = height, region = 'increase',
                             type_of_power = 'spurious')

    power_matrix[2 * i - 1, ((k - 1) * 5 + 3):(k * 5 + 1)] <- power_overall
    power_matrix[2 * i, ((k - 1) * 5 + 3):(k * 5 + 1)] <- power_spurious
    i <- i + 1 
  }
  k <- k + 1
}

print.xtable(xtable(power_matrix, digits = c(3), align = paste(replicate(5*length(different_a1) + 2, "c"), collapse = "")),
             type="latex", file=paste0("plots/power_table.tex"), include.colnames = FALSE)

############################
#Calculating rowwise power#
############################
sim.design <- 'bump'  # type of trend function m()
height.neg <- 0.85    # height of the bump in case of sim.design = 'bump'. This height is used to calculate power for negative a1     
height.pos <- 2.65    # height of the bump in case of sim.design = 'bump'. This height is used to calculate power for positive a1 

for (a1 in different_a1){
  if (a1 < 0) {height = height.neg}
  else {height = height.pos}
  
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
  pdf(pdffilename, width = 6, height = 10, paper = 'special')

  par(mfrow = c(3, 1), cex = 1,  tck = -0.025) #Setting the layout of the graphs
  
  par(mar = c(3.5, 3.5, 0, 0)) #Margins for each plot
  par(oma = c(0, 0.2, 0.2, 0.2)) #Outer margins 
  
  set.seed(T)
  data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = sim.design, slope.fac = height)
  data           <- data.simulated$data
  trend          <- data.simulated$trend

  plot(seq(1/T, 1, by = 1/T), data, ylab = "", xlab = "", col= 'grey', type = 'l', lty = 1,mgp=c(1.8,0.5,0))
  lines(seq(1/T, 1, by = 1/T), trend, lwd = 2)
  
  plot(x = h.grid, y = result$power_ms*100, ylim=c(0, 70), yaxp = c(0, 60, 4),mgp=c(1.8,0.5,0), type="l", lty=1, xaxt='n', ylab = "Percentage (%)", xlab='bandwidth h')
  points(x = h.grid, y = result$power_ms*100, pch=19, cex = 0.8)
  
  lines(x = h.grid, y = result$power_uncor*100, lwd=1.5, lty = 'dashed') 
  lines(x = h.grid, y = result$power_rows*100, lwd=1.5, lty = 'dotted') 
  lines(x = h.grid, y = result$power_SiZer*100, lwd=1.5) 
  title(main = "Rowwise power", line = -1.5)
  
  
  axis(1, at=h.grid, mgp=c(1.8,0.5,0))
  legend('topleft', cex = 0.8, bty = "n", legend = c(expression(italic(T)[MS]), expression(italic(T)[UC]), expression(italic(T)[RW]), expression(italic(T)[SiZer])),
         pch = c(19, NA, NA, NA), lty = c('solid', 'dashed', 'dotted', 'solid'), y.intersp=1.25)
  
  plot(x = h.grid, y = result_spurious$power_ms*100, ylim=c(0, 9), yaxp = c(0, 8, 4), type="l",mgp=c(1.8,0.5,0), lty=1, xaxt='n', ylab = "Percentage (%)", xlab='bandwidth h', lwd = 2)
  points(x = h.grid, y = result_spurious$power_ms*100, pch=19, cex = 0.8)
  
  lines(x = h.grid, y = result_spurious$power_uncor*100, lwd=2, lty = 'dashed') 
  lines(x = h.grid, y = result_spurious$power_rows*100, lwd=2, lty = 'dotted')
  lines(x = h.grid, y = result_spurious$power_SiZer*100, lwd=2)
  title(main = "Rowwise spurious power", line = -1.5)
  
  axis(1, at=h.grid,mgp=c(1.8,0.5,0))
  legend('topleft', cex = 0.8, bty = "n", legend = c(expression(italic(T)[MS]), expression(italic(T)[UC]), expression(italic(T)[RW]), expression(italic(T)[SiZer])),
         pch = c(19, NA, NA, NA), lty = c('solid', 'dashed', 'dotted', 'solid'), y.intersp=1.25)
  
  dev.off()
}

##########################
#Blocks and sine examples#
##########################

sim.design <- c('blocks', 'sine') # type of trend function m()

for (a1 in different_a1){
  for (sim.design_ in sim.design){
    set.seed(1)
    if (sim.design_ == 'sine') {
      sigma_eta = sqrt((1 -a1^2) )
    } else {
      sigma_eta = sqrt((1 - a1^2) * 0.01)
    }

    #Simulating the time series
    data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = sim.design_)
    data           <- data.simulated$data
    trend          <- data.simulated$trend
    sigma_true     <- data.simulated$sigma

    #Construct grid
    grid      <- grid_construction(T)
    gset      <- grid$gset
    u.grid    <- sort(unique(gset[,1]))
    h.grid    <- sort(unique(gset[,2]))
    autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))

    ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
    deletions <- ess$del
    grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)

    sizer.wghts   <- SiZer_weights(T=T, grid=grid)
    sizer.std     <- SiZer_std(weights=sizer.wghts, autocov=autocov, T)
    sizer.quants  <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
    SiZer.values  <- sizer.wghts %*% data
    sizer.vals    <- as.vector(SiZer.values)

    # Compute kernel weights and critical value for multiscale test
    wghts  <- kernel_weights(T=T, grid=grid)
    quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
    stats  <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_true, grid=grid)
    vals   <- stats$values

    test.res      <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
    SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)

    pdffilename <- paste0("plots/SiZer_map_T_", T, "_", sim.design_, "_a1_", a1*100, ".pdf")
    pdf(pdffilename, width = 6, height = 16, paper = 'special')

    par(mfrow = c(5, 1), cex = 1,  tck = -0.025) #Setting the layout of the graphs

    par(mar = c(1, 4, 2, 0)) #Margins for each plot
    par(oma = c(0.5, 0.2, 0, 0.2)) #Outer margins

    plot(seq(1/T, 1, by = 1/T), data, xlim = c(0, 1), ylab = "", xlab = "", col= 'grey', type = 'l', lty = 1,mgp=c(1.8,0.5,0), xaxt = "n")
    lines(seq(1/T, 1, by = 1/T), trend, lwd = 2)

    SiZermap(u.grid, h.grid, test.res$test_ms, plot.title = expression(italic(T)[MS]))
    axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = FALSE, mgp=c(1,0.5,0))

    SiZermap(u.grid, h.grid, test.res$test_uncor, plot.title = expression(italic(T)[UC]))
    axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = FALSE, mgp=c(1,0.5,0))

    SiZermap(u.grid, h.grid, test.res$test_rows, plot.title = expression(italic(T)[RW]))
    axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = FALSE, mgp=c(1,0.5,0))

    SiZermap(u.grid, h.grid, SiZer_results$test, plot.title = expression(italic(T)[SiZer]))

    axis(1, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(NA, 0.2, 0.4, 0.6, 0.8, NA), mgp=c(1,0.5,0))
    dev.off()
  }
}