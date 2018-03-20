data_analysis <- function(alpha, y_data, test_problem = "zero", kernel_t = "nw"){
  #Recoding testing problem and type of kernel estimator 
  if (test_problem == "zero"){
    kernel_ind = 1
  } else if (test_problem == "constant"){
    kernel_ind = 2
  } else {
    print('Given testing problem is currently not supported')
  }
  
  if (kernel_t == "nw"){
    quantile_function = calculating_gaussian_quantile
    statistic_function = psihat_statistic
    defining_set = creating_g_set
  } else if (kernel_t == "ll"){
    quantile_function = calculating_gaussian_quantile_ll
    statistic_function = psihat_statistic_ll
    defining_set = creating_g_set_ll
  } else {
    print('Given method is currently not supported')
  }
  
  T_data = length(data)
  
  #Tuning parameters
  L1 <- floor(sqrt(T_data))
  L2 <- floor(2 * sqrt(T_data))
  sigmahat_data <- estimating_sigma_for_AR1(data, L1, L2)[[1]]
  
  #Parameters for plotting
  grid_points <- seq(from = 1/T_data, to = 1, length.out = T_data) #grid points for plotting and estimating

  ##########################################
  #Calculating gaussian quantile for T_data#
  ##########################################
  
  g_t_set           <- defining_set(T_data)
  gaussian_quantile <- quantile_function(T_data, g_t_set, kernel_ind, alpha)
  
  
  #########################################
  #Calculating the statistic for real data#
  #########################################
  
  result                 <- statistic_function(data, g_t_set, kernel_ind, sigmahat_data)
  g_t_set_with_values    <- result[[1]]
  psihat_statistic_value <- result[[2]]
  
  #And now the testing itself
  if (psihat_statistic_value > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    
    #The collection of intervals where the corrected test statistic lies above the critical value (as in (2.6))
    a_t_set <- subset(g_t_set_with_values, values > gaussian_quantile, select = c(u, h, values))
    p_t_set <- data.frame('startpoint' = a_t_set$u - a_t_set$h, 'endpoint' = a_t_set$u + a_t_set$h, 'values' = a_t_set$values)
    p_t_set <- choosing_minimal_intervals(p_t_set)
    
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
    
  pdffilename = paste0("Output/threegraphics_testproblem_", test_problem, "_method_", kernel_t, ".pdf")
  pdf(pdffilename,width=10,height=10,paper="special")

  par(mfrow = c(3,1)) #Setting the layout of the graphs
  par(mar = c(0, 0.5, 0, 0.5)) #Margins for each plot
  par(oma = c(0, 2, 0, 0)) #Outer margins
  
  # Plotting the real data
  plot(grid_points, data, type = "l") 
  
  ###################################################################
  #Calculating smoothed curve for the data using Epanechnikov kernel#
  ###################################################################
  h <- c(0.01, 0.05, 0.1, 0.15, 0.2)
  plot(NA, xlim = c(0, 1), ylim = c(-1.5, 1.5))
  for (i in 1:5){
    #This part plots kernel smoothers for different bandwidths (all on one plot).
    smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(data, grid_points, h[i]))
    lines(grid_points, smoothed_curve, lty = i) 
  }
  legend(-1/T_tempr, 1.3, legend=h, lty = i, ncol=1)
  
  #Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
  ymaxlim = max(p_t_set$values)
  yminlim = min(p_t_set$values)
  plot(NA, xlim=c(0,1), ylim = c(yminlim - 1, ymaxlim + 1), xlab="x", ylab="y")
  segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set[['endpoint']], p_t_set[['values']])
  
  #title(main = "Plots for normalized yearly temperature data for England",  outer = TRUE)
  dev.off()
}