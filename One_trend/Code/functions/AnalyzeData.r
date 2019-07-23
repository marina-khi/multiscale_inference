AnalyzeData <- function(data, ts_start, filename_table, filename_plot, axis_at, axis_labels, xaxp, yaxp, 
                          alpha = 0.05, SimRuns = 1000, order = 1, q = 25, r.bar = 10, plot_SiZer = 'yes',
                          sigma_supplied = 'no', sigma = NA){
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
  #   plot_SiZer:     If "yes", then the function produces SiZer map for SiZer method. Since SiZer needs an estimate 
  #                   for the whole autocovariance function of the error process and sometimes this estimate is not availabe,
  #                   in this cases SiZer map for SiZer can not be plotted. Then set the value of plot_SiZer to 'no'. Default is "yes".
  #   sigma_supplied: If 'yes', then the estimator of the square root of the long-run variance for the data needs to be supplied.
  #                   Default is 'no'.
  #   sigma:          Value of the square root of the long-run variance estimated by other methods. Needed only if sigma_supplied = 'yes'.
  #
  # Outputs:
  #   physical files  

  #Load necessary functions  
  source("functions/ConstructGrid.r")
  source("functions/multiscale_statistics.r")
  source("functions/multiscale_quantiles.r")
  source("functions/multiscale_testing.r")
  source("functions/minimal_intervals.r")
  source("functions/SiZer_functions.r")
  sourceCpp("functions/kernel_weights.cpp")
  sourceCpp("functions/SiZer_functions.cpp")
  
    
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
  
  #In some cases we are plotting the SiZer map for SiZer method
  if (plot_SiZer == 'yes'){
    gamma <- AR_acf(AR.struc$ahat, AR.struc$vareta, Tlen)
    
    ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=Tlen, autocov=gamma)
    deletions <- ess$del
    grid      <- grid_construction(T=Tlen, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
    
  
    sizer.wghts   <- SiZer_weights(T=Tlen, grid=grid)
    sizer.std     <- SiZer_std(weights=sizer.wghts, autocov=gamma, Tlen)
    sizer.quants  <- SiZer_quantiles(alpha=alpha, T=Tlen, grid=grid, autocov=gamma)
    SiZer.values  <- sizer.wghts %*% data
    sizer.vals    <- as.vector(SiZer.values)
    SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)
    
    test.res      <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
  
    pdf(filename_plot, width=7, height=10, paper="special")
    par(mfrow = c(4,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  } else {
    pdf(filename_plot, width=7, height=7.5, paper="special")
    par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  }
  
  #Parameters for plotting
  grid_time <- seq(from = ts_start, to = ts_end, length.out = Tlen) #grid points for plotting 
  p_t_set$plottingindex <- (1:nrow(p_t_set))/nrow(p_t_set)
  
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
  segments(p_t_set[['startpoint']], p_t_set$plottingindex, p_t_set$endpoint, p_t_set$plottingindex)
  #abline(h = quant.ms, lty = 2)

  #SiZer 
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  SiZermap(u.grid, h.grid, test.res$test_ms, plot.title = expression((c) ~ SiZer ~ map ~ 'for' ~ italic(T)[MS]))
  
  #SiZer map for SiZer method
  if (plot_SiZer == 'yes'){
    axis(1, at=axis_at, labels = FALSE, mgp=c(1,0.5,0))
    SiZermap(u.grid, h.grid, SiZer_results$test, plot.title = expression((d) ~ SiZer ~ map ~ 'for' ~ italic(T)[SiZer]))
  }
  axis(1, at=axis_at, labels = axis_labels, mgp=c(1,0.5,0))
  
  dev.off()
}