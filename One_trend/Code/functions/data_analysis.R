data_analysis <- function(data, ts_start, filename_table, filename_plot, axis_at, axis_labels, xaxp, yaxp, 
                          alpha = 0.05, SimRuns = 1000, order = 1, q = 25, r.bar = 10, plot_SiZer = 'yes',
                          sigma_supplied = 'no', sigma = NA){
  # Function that carries out the multiscale test for the given time series including all the necessary calculations of   
  # multiscale statistics and critical values. 
  #
  # Arguments:
  # data        time series needed to be analysed
  # ts_start    first point in time of the data. Needed for the plots
  # alpha       significance level     
  #  
  # Outputs: physical files  
  # test_ms     matrix of test results for our multiscale test defined in Section 3
  #             test_ms[i,j] = -1: test rejects the null for the j-th location u and the 
  #                                i-th bandwidth h and indicates a decrease in the trend
  #             test_ms[i,j] = 0:  test does not reject the null for the j-th location u  
  #                                and the i-th bandwidth h 
  #             test_ms[i,j] = 1:  test rejects the null for the j-th location u and the 
  #                                i-th bandwidth h and indicates an increase in the trend
  #             test_ms[i,j] = 2:  no test is carried out at j-th location u and i-th 
  #                                bandwidth h (because the point (u,h) is excluded from  
  #                                the grid as specified by the 'deletions'-option in the
  #                                function 'grid_construction')  
  # test_uncor  matrix of test results for the uncorrected version of the multiscale test 
  # test_order  matrix of test results for the version of the multiscale test from Section ??
  # test_rows   matrix of test results for the rowwise version of the multiscale test 
  
  
  Tlen <- length(data)         #length of the data
  ts_end <- ts_start + Tlen -1 #the last point of time series
  
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
  
  #Compute test results
  test.res   <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
  gset$test  <- as.vector(t(test.res$test_ms))
  gset$vals2 <- as.vector(abs(vals) - correct) 
  quant.ms   <- test.res$quant.ms
  
  #Produce minimal intervals
  a_t_set <- subset(gset, test == 1, select = c(u, h, vals2))
  p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*Tlen + ts_start, 'endpoint' = (a_t_set$u + a_t_set$h)*Tlen + ts_start, 'values' = a_t_set$vals2)
  p_t_set <- subset(p_t_set, endpoint <= ts_end, select = c(startpoint, endpoint, values)) 
  p_t_set <- choosing_minimal_intervals(p_t_set)
  
  print.xtable(xtable(subset(p_t_set, select = c(startpoint, endpoint)), digits = c(0)), type="latex", file=filename_table)
  
  #Parameters for plotting
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
    pdf(filename_plot, width=7, height=7.5, paper="special")
    par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  }
  
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  
  # Plotting the real data
  plot(grid_time, xlim = c(ts_start, ts_end), data, type = "l", mgp=c(1,0.5,0), xaxp = xaxp, xaxs='i')
  title(main = expression((a) ~ observed ~ temperature ~ time ~ series), line = 1)
  
  par(mar = c(0.5, 0.5, 3.5, 0)) #Margins for each plot
  #Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
  ymaxlim = max(p_t_set$values)
  yminlim = min(min(p_t_set$values), quant.ms)
  plot(NA, xlim=c(ts_start,ts_end), xaxt = "n",  ylim = c(yminlim - 0.1, ymaxlim + 0.1), yaxp  = yaxp, mgp=c(2,0.5,0))
  title(main = expression((b) ~ minimal ~ intervals ~ produced ~ by ~ italic(T)[MS]), line = 1)
  segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set$endpoint, p_t_set[['values']])
  abline(h = quant.ms, lty = 2)

  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  SiZermap(u.grid, h.grid, test.res$test_ms, plot.title = expression((c) ~ SiZer ~ map ~ 'for' ~ italic(T)[MS]))
  #title(main =, line = 1)
    
  if (plot_SiZer == 'yes'){
    axis(1, at=axis_at, labels = FALSE, mgp=c(1,0.5,0))
    SiZermap(u.grid, h.grid, SiZer_results$test, plot.title = expression((d) ~ SiZer ~ map ~ 'for' ~ italic(T)[SiZer]))
  }
  axis(1, at=axis_at, labels = axis_labels, mgp=c(1,0.5,0))
  
  dev.off()
  #return(list(gset, quant.ms))
}