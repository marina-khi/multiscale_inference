data_analysis <- function(data, ts_start, ts_end, filename_table, filename_plot, axis_at, axis_labels, xaxp, yaxp, 
                          alpha = 0.05, SimRuns = 1000, order = 1, q = 25, r.bar = 10, plot_SiZer = 'yes',
                          sigma_supplied = 'no', sigma = NA){
  Tlen <- length(data)
  
  #estimate the long-run variance
  AR.struc  <- AR_lrv(data=data, q=q, r.bar=r.bar, p=order)
  sigma_hat <- sqrt(AR.struc$lrv)
  
  #Construct grid
  grid    <- grid_construction(Tlen)
  gset    <- grid$gset
  u.grid  <- sort(unique(gset[,1]))
  h.grid  <- sort(unique(gset[,2]))
  correct <- sqrt(2*log(1/(2*gset[,2])))
  
  # Compute kernel weights and critical value for multiscale test
  wghts <- kernel_weights(T=Tlen, grid=grid)
  if (sigma_supplied == 'yes'){
    stats <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma, grid=grid)
  } else {
    stats <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid)
  }
  vals  <- stats$values
  
  #Compute the quantile
  filename = paste0("quantiles/distr_T_", Tlen,".RData")
  if(!file.exists(filename)) {
    quants <- multiscale_quantiles(T=Tlen, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
    save(quants, file = filename)
  } else {
    load(filename)
  }
  
  # Select (1-alpha) quantile of the multiscale statistic under the null
  probs.ms <- as.vector(quants$quant_ms[1,])
  quant.ms <- as.vector(quants$quant_ms[2,])
  
  if(sum(probs.ms == (1-alpha)) == 0){
    pos.ms <- which.min(abs(probs.ms-(1-alpha)))
  } else {
    pos.ms <- which.max(probs.ms == (1-alpha)) 
  }
  quant.ms    <- quant.ms[pos.ms]
  
  #Compute test results
  vals2            <- abs(vals) - correct
  gset_result      <- cbind(gset, vals, vals2)
  gset_result$test <- (gset_result$vals2 > quant.ms) * sign(gset_result$vals)
  
  #Produce minimal intervals
  a_t_set <- subset(gset_result, test == 1, select = c(u, h, vals2))
  p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*Tlen + ts_start, 'endpoint' = (a_t_set$u + a_t_set$h)*Tlen + ts_start, 'values' = a_t_set$vals2)
  p_t_set <- subset(p_t_set, endpoint <= ts_end, select = c(startpoint, endpoint, values)) 
  p_t_set <- choosing_minimal_intervals(p_t_set)
  
  print.xtable(xtable(subset(p_t_set, select = c(startpoint, endpoint)), digits = c(0)), type="latex", file=filename_table)
  
  #Parameters for plotting
  grid_points <- seq(from = 1/Tlen, to = 1, length.out = Tlen) #grid points for estimating
  grid_time <- seq(from = ts_start, to = ts_end, length.out = Tlen) #grid points for plotting 
  
  
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
    test.res      <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
    
    pdf(filename_plot, width=7, height=10, paper="special")
    par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  }
  
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  
  # Plotting the real data
  plot(grid_time, xlim = c(ts_start, ts_end), data, type = "l", mgp=c(1,0.5,0), xaxp = xaxp, xaxs='i') 
  
  
  #Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
  ymaxlim = max(p_t_set$values)
  yminlim = min(p_t_set$values)
  plot(NA, xlim=c(ts_start,ts_end), xaxt = "n",  ylim = c(yminlim - 0.3, ymaxlim + 0.3), yaxp  = yaxp, mgp=c(2,0.5,0))
  segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set$endpoint, p_t_set[['values']])
  abline(h = quant.ms, lty = 2)
  
  SiZermap(u.grid, h.grid, test.res$test_ms, plot.title = expression(italic(T)[MS]))

  if (plot_SiZer == 'yes'){
    axis(1, at=axis_at, labels = FALSE, mgp=c(1,0.5,0))
    SiZermap(u.grid, h.grid, SiZer_results$test, plot.title = expression(italic(T)[SiZer]))
  }
  axis(1, at=axis_at, labels = axis_labels, mgp=c(1,0.5,0))
  
  dev.off()
}