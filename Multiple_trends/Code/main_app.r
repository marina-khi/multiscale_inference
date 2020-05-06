# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

library(Rcpp)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

# The following file contains functions that are used to estimate the AR parameters, in particular, to compute
# the estimators of the coefficients, a_1, ... ,a_p, the estimator of the variance of the innovation term, sigma_eta^2 
# and the estimator of the long-run variance sigma^2.
source("functions/long_run_variance.r")

#Load necessary functions  
source("functions/ConstructGrid.r")
source("functions/functions.r")
source("functions/multiscale_statistics.r")
source("functions/multiscale_quantiles.r")
source("functions/multiscale_testing.r")
source("functions/minimal_intervals.r")
sourceCpp("functions/kernel_weights.cpp")


###############################################
#Analysis of UK annual temperature time series#
###############################################

alpha   <- 0.05 # alpha for calculating quantiles
SimRuns <- 5000 # number of simulation runs to produce critical values

#Loading the real data for yearly temperature in England
temperature  <- read.table("data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
T_tempr      <- length(yearly_tempr)

#Order selection
q <- 10:35
r <- 10:15
criterion_matrix <- expand.grid(q = q, r = r)

criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))

for (i in 1:nrow(criterion_matrix)){
  FPE <- c()
  AIC <- c()
  AICC <- c()
  SIC <- c()
  HQ <- c()
  
  different_orders <- (1:9)
  
  for (order in different_orders){
    AR.struc      <- AR_lrv(data=yearly_tempr, q=criterion_matrix$q[[i]], r.bar=criterion_matrix$r[[i]], p=order)
    sigma_eta_hat <- sqrt(AR.struc$vareta)
    FPE <- c(FPE, (sigma_eta_hat^2 * (T_tempr + order)) / (T_tempr - order))
    AIC <- c(AIC, T_tempr * log(sigma_eta_hat^2) + 2 * order)
    AICC <- c(AICC, T_tempr * log(sigma_eta_hat^2) + T_tempr* (1 + order / T_tempr)/(1 - (order +2)/T_tempr))
    SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(T_tempr) / T_tempr)
    HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(T_tempr)) / T_tempr)
  }
  criterion_matrix$FPE[[i]]  <- which.min(FPE)
  criterion_matrix$AIC[[i]]  <- which.min(AIC)
  criterion_matrix$AICC[[i]] <- which.min(AICC)
  criterion_matrix$SIC[[i]]  <- which.min(SIC)
  criterion_matrix$HQ[[i]]   <- which.min(HQ)
}

write.csv(criterion_matrix,"plots/criterion_matrix_UK.csv", row.names = FALSE)

#Setting tuning parameters for testing
order <- 2
q     <- 25
r.bar <- 10

#And finally testing itself!
data = yearly_tempr
ts_start= 1659
filename_table = "plots/minimal_intervals_UK.tex"
filename_plot = "plots/temperature_UK.pdf"
axis_at = seq(17/T_tempr, 367/T_tempr, by = 50/T_tempr)
axis_labels = seq(1675, 2025, by = 50)
xaxp = c(1675, 2025, 7)
yaxp = c(1.9, 2.2, 3)
sigma_supplied = 'no'
sigma = NA
  # Function that performs different multiscale tests on the data to determine where the underlying trend
  # is increasing or decreasing. It produces all the plots from the paper such as SiZer map or the plot of minimal intervals
  #
  # Arguments:
  #   data:           Time series needed to be analysed.
  #   ts_start:       First point in time of the data. Needed for the plots.
  #   filename_table: Name (and path) of the .tex file which will contain the table with the minimal intervals.
  #   filename_plot:  Name (and path) of the .pdf file which will contain the plots with minimal intervals and SiZer maps.
  #   axis_at, axis_labels, xaxp, yaxp:
  #                   Parameters needed for plotting.
  #   alpha:          Significance level. Default is 0.05.     
  #   SimRuns:        Number of simulations for calculating gaussian quantile for the multiscale tests. Default is 1000.
  #   order:          Order of the AR process. Default is 1.
  #   q, r.bar:       Tuning parameters for estimating the long-run variance from Section 4. Default are 25 and 10. 
  #   sigma_supplied: If 'yes', then the estimator of the square root of the long-run variance for the data needs to be supplied.
  #                   Default is 'no'.
  #   sigma:          Value of the square root of the long-run variance estimated by other methods. Needed only if sigma_supplied = 'yes'.
  #
  # Outputs:
  #   physical files  
  
Tlen <- length(data)         #length of the data
ts_end <- ts_start + Tlen -1 #the last point of time series
  
#Estimate the long-run variance
AR.struc  <- AR_lrv(data=data, q=q, r.bar=r.bar, p=order)
sigma_hat <- sqrt(AR.struc$lrv)
  
#Construct grid
grid    <- grid_construction(Tlen)
gset    <- grid$gset
u.grid  <- sort(unique(gset[,1]))
h.grid  <- sort(unique(gset[,2]))
correct <- sqrt(2*log(1/(2*gset[,2])))
  
# Compute kernel weights and critical value for multiscale test
Tlen                   <- as.integer(Tlen) 
N                      <- as.integer(dim(grid$gset)[1])
gset_cpp               <- as.matrix(grid$gset)
gset_cpp               <- as.vector(gset_cpp) 
storage.mode(gset_cpp) <- "double"
  
wghts <- matrix(kernel_weights(Tlen, gset_cpp, N), ncol = Tlen, byrow = TRUE)
  
  if (sigma_supplied == 'yes'){
    stats <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma, grid=grid)
  } else {
    stats <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid)
  }
  vals  <- stats$values
  
  #Compute the quantile for the multiscale method
  filename = paste0("quantiles/distr_T_", Tlen,".RData")
  if(!file.exists(filename)) {
    quants <- multiscale_quantiles(T=Tlen, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
    save(quants, file = filename)
  } else {
    load(filename)
  }
  
  #Compute test results for the multiscale method
  test.res   <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
  gset$test  <- as.vector(t(test.res$test_ms))
  gset$vals2 <- as.vector(abs(vals) - correct) 
  quant.ms   <- test.res$quant.ms
  
  #Produce minimal intervals (Here - only the increases! But this is because we do not have decreases for our applications)
  a_t_set <- subset(gset, test == 1, select = c(u, h, vals2))
  p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*Tlen + ts_start, 'endpoint' = (a_t_set$u + a_t_set$h)*Tlen + ts_start, 'values' = a_t_set$vals2)
  p_t_set <- subset(p_t_set, endpoint <= ts_end, select = c(startpoint, endpoint, values)) 
  p_t_set <- choosing_minimal_intervals(p_t_set)
  
  #Writing down thetable with minimal intervals for T_ms
  print.xtable(xtable(subset(p_t_set, select = c(startpoint, endpoint)), digits = c(0)), type="latex", file=filename_table)
  
  #postscript(filename_plot, paper="special", width=7, height=6.6, horizontal=FALSE)
  pdf(filename_plot, width=7, height=6.6, paper="special")
  layout(matrix(c(1,2,3),ncol=1), widths=c(3,3,3), heights=c(1,0.8,1), TRUE) # Setting the layout of the graphs

  #Parameters for plotting
  grid_time <- seq(from = ts_start, to = ts_end, length.out = Tlen) #grid points for plotting 
  p_t_set$plottingindex <- (1:nrow(p_t_set))/nrow(p_t_set)
  
  par(cex = 1.1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  
  # Plotting the real data
  plot(grid_time, xlim = c(ts_start, ts_end), data, type = "l", mgp=c(1,0.5,0), xaxp = xaxp, xaxs='i')
  title(main = expression((a) ~ observed ~ temperature ~ time ~ series), line = 1)
  
  par(mar = c(0.5, 0.5, 3.5, 0)) #Margins for each plot
  
  #Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
  #ymaxlim = max(p_t_set$values)
  #yminlim = min(min(p_t_set$values), quant.ms)
  plot(NA, xlim=c(ts_start,ts_end), xaxt = "n",  ylim = c(0, 1 + 1/nrow(p_t_set)), yaxp  = yaxp, mgp=c(2,0.5,0))
  title(main = expression((b) ~ minimal ~ intervals ~ produced ~ by ~ italic(T)[MS]), line = 1)
  segments(p_t_set[['startpoint']], p_t_set$plottingindex, p_t_set$endpoint, p_t_set$plottingindex, lwd = 2)
  #abline(h = quant.ms, lty = 2)
  
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  SiZermap(u.grid, h.grid, test.res$test_ms, plot.title = expression((c) ~ SiZer ~ map ~ 'for' ~ italic(T)[MS]))
  
  axis(1, at=axis_at, labels = axis_labels, mgp=c(1,0.5,0))
  
  dev.off()