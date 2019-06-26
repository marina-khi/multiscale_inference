data_analysis <- function(alpha, y_data, test_problem, sigmahat, pdffilename){
  #Recoding testing problem and type of kernel estimator 
  if (test_problem == "zero"){
    kernel_ind = 1
  } else if (test_problem == "constant"){
    kernel_ind = 2
  } else {
    print('Given testing problem is currently not supported')
  }

  T_data = length(y_data)
  
  ##########################################
  #Calculating gaussian quantile for T_data#
  ##########################################
  
  g_t_set           <- creating_g_set(T_data)
  gaussian_quantile <- calculating_gaussian_quantile_ll(T_data, g_t_set, test_problem, kernel_ind, alpha)
  
  
  #########################################
  #Calculating the statistic for real data#
  #########################################
  
  result                 <- psihat_statistic_ll(y_data, g_t_set, kernel_ind, sigmahat)
  g_t_set_with_values    <- result[[1]]
  psihat_statistic_value <- result[[2]]
  
  #And now the testing itself
  if (psihat_statistic_value > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    
    #The collection of intervals where the corrected test statistic lies above the critical value (as in (2.6))
    a_t_set <- subset(g_t_set_with_values, values > gaussian_quantile, select = c(u, h, values))
    p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*T_data + 1659, 'endpoint' = (a_t_set$u + a_t_set$h)*T_data + 1659, 'values' = a_t_set$values)
    p_t_set <- subset(p_t_set, endpoint <= 2017, select = c(startpoint, endpoint, values)) 
    p_t_set <- choosing_minimal_intervals(p_t_set)

    print.xtable(xtable(subset(p_t_set, select = c(startpoint, endpoint)), digits = c(0)), type="latex", file="Paper/Plots/minimal_intervals.tex")

    #The collection of intervals where the values_with_sign > gaussian_quantile + lambda (as in (2.6))
    #a_t_set_plus <- subset(g_t_set_with_values, values_with_sign > gaussian_quantile + lambda, select = c(u, h, values_with_sign))
    #p_t_set_plus <- data.frame('startpoint' = a_t_set_plus$u - a_t_set_plus$h, 'endpoint' = a_t_set_plus$u + a_t_set_plus$h, 'values' = a_t_set_plus$values_with_sign)
    #p_t_set_plus <- choosing_minimal_intervals(p_t_set_plus)
    
    #The collection of intervals where the -values_with_sign > gaussian_quantile + lambda (as in (2.6))
    #a_t_set_minus <- subset(g_t_set_with_values, -values_with_sign > gaussian_quantile + lambda, select = c(u, h, values_with_sign))
    #p_t_set_minus <- data.frame('startpoint' = a_t_set_minus$u - a_t_set_minus$h, 'endpoint' = a_t_set_minus$u + a_t_set_minus$h, 'values' = a_t_set_minus$values_with_sign)
    #p_t_set_minus <- choosing_minimal_intervals(p_t_set_minus)
    
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  }
  
  #Parameters for plotting
  grid_points <- seq(from = 1/T_data, to = 1, length.out = T_data) #grid points for estimating
  grid_time <- seq(from = 1659, to = 2017, length.out = T_data) #grid points for plotting 
    
  pdf(pdffilename, width=10, height=6, paper="special")

  par(mfrow = c(2,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  
  # Plotting the real data
  plot(grid_time, xlim = c(1659, 2019), y_data, type = "l", mgp=c(2,0.5,0), xaxp = c(1675, 2025,7)) 
  
  #Calculating smoothed curve for the data using local linear estimator#
  h <- c(0.01, 0.05, 0.1, 0.15, 0.2)
  plot(NA, xlim = c(1659, 2019), ylim = c(7.5, 11), yaxp  = c(8, 10, 2), xaxp = c(1675, 2025,7), mgp=c(2,0.5,0))
  for (i in 1:5){
   #This part plots kernel smoothers for different bandwidths (all on one plot).
   smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(y_data, grid_points, h[i]))
   lines(grid_time, smoothed_curve, lty = i) 
  }
  legend(1950, 9, legend=c("h = 0.01", "h = 0.05", "h = 0.10", "h = 0.15", "h = 0.2"), lty = 1:5, cex = 0.95, ncol=1)
  
  #Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
  ymaxlim = max(p_t_set$values)
  yminlim = min(p_t_set$values)
  plot(NA, xlim=c(1659,2019), xaxp = c(1675, 2025,7),  ylim = c(yminlim - 0.2, ymaxlim + 0.2), yaxp  = c(1.75, 2.5, 3), mgp=c(2,0.5,0))
  segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set$endpoint, p_t_set[['values']])
  abline(h = gaussian_quantile, lty = 2)
  
  dev.off()
}


data_analysis_global <- function(alpha, y_data, test_problem, sigmahat, pdffilename){
  #Recoding testing problem and type of kernel estimator 
  if (test_problem == "zero"){
    kernel_ind = 1
  } else if (test_problem == "constant"){
    kernel_ind = 2
  } else {
    print('Given testing problem is currently not supported')
  }
  
  T_data = length(y_data)
  
  ##########################################
  #Calculating gaussian quantile for T_data#
  ##########################################
  
  g_t_set           <- creating_g_set(T_data)
  gaussian_quantile <- calculating_gaussian_quantile_ll(T_data, g_t_set, test_problem, kernel_ind, alpha)
  
  
  #########################################
  #Calculating the statistic for real data#
  #########################################
  
  result                 <- psihat_statistic_ll(y_data, g_t_set, kernel_ind, sigmahat)
  g_t_set_with_values    <- result[[1]]
  psihat_statistic_value <- result[[2]]
  
  #And now the testing itself
  if (psihat_statistic_value > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    
    #The collection of intervals where the corrected test statistic lies above the critical value (as in (2.6))
    a_t_set <- subset(g_t_set_with_values, values > gaussian_quantile, select = c(u, h, values))
    p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*T_data + 1850, 'endpoint' = (a_t_set$u + a_t_set$h)*T_data + 1850, 'values' = a_t_set$values)
    #p_t_set <- subset(p_t_set, endpoint <= 2017, select = c(startpoint, endpoint, values)) 
    p_t_set <- choosing_minimal_intervals(p_t_set)
    
    print.xtable(xtable(subset(p_t_set, select = c(startpoint, endpoint)), digits = c(0)), type="latex", file="Paper/Plots/minimal_intervals.tex")
    
    #The collection of intervals where the values_with_sign > gaussian_quantile + lambda (as in (2.6))
    #a_t_set_plus <- subset(g_t_set_with_values, values_with_sign > gaussian_quantile + lambda, select = c(u, h, values_with_sign))
    #p_t_set_plus <- data.frame('startpoint' = a_t_set_plus$u - a_t_set_plus$h, 'endpoint' = a_t_set_plus$u + a_t_set_plus$h, 'values' = a_t_set_plus$values_with_sign)
    #p_t_set_plus <- choosing_minimal_intervals(p_t_set_plus)
    
    #The collection of intervals where the -values_with_sign > gaussian_quantile + lambda (as in (2.6))
    #a_t_set_minus <- subset(g_t_set_with_values, -values_with_sign > gaussian_quantile + lambda, select = c(u, h, values_with_sign))
    #p_t_set_minus <- data.frame('startpoint' = a_t_set_minus$u - a_t_set_minus$h, 'endpoint' = a_t_set_minus$u + a_t_set_minus$h, 'values' = a_t_set_minus$values_with_sign)
    #p_t_set_minus <- choosing_minimal_intervals(p_t_set_minus)
    
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  }
  
  #Parameters for plotting
  grid_points <- seq(from = 1/T_data, to = 1, length.out = T_data) #grid points for estimating
  grid_time <- seq(from = 1850, to = 2015, length.out = T_data) #grid points for plotting 
  
  pdf(pdffilename, width=10, height=10, paper="special")
  
  par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  
  # Plotting the real data
  plot(grid_time, xlim = c(1850, 2020), y_data, type = "l", mgp=c(2,0.5,0), xaxp = c(1850, 2020,5)) 
  
  #Calculating smoothed curve for the data using local linear estimator#
  h <- c(0.01, 0.05, 0.1, 0.15, 0.2)
  plot(NA, xlim = c(1850, 2020), ylim = c(-1, 1), yaxp  = c(-1, 1, 3), xaxp = c(1850, 2020,5), mgp=c(2,0.5,0))
  for (i in 1:5){
    #This part plots kernel smoothers for different bandwidths (all on one plot).
    smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(y_data, grid_points, h[i]))
    lines(grid_time, smoothed_curve, lty = i) 
  }
  legend(1950, 0.5, legend=c("h = 0.01", "h = 0.05", "h = 0.10", "h = 0.15", "h = 0.2"), lty = 1:5, cex = 0.95, ncol=1)
  
  #Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
  ymaxlim = max(p_t_set$values)
  yminlim = min(p_t_set$values)
  plot(NA, xlim=c(1850,2020), xaxp = c(1850, 2020,5),  ylim = c(yminlim - 0.2, ymaxlim + 0.2), yaxp  = c(-1, 1, 3), mgp=c(2,0.5,0))
  segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set$endpoint, p_t_set[['values']])
  abline(h = gaussian_quantile, lty = 2)
  
  dev.off()
}
