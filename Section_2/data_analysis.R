data_analysis <- function(alpha, y_data, test_problem, kernel_t){
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
  
  T_data = length(y_data)
  
  #Tuning parameters
  L1 <- floor(sqrt(T_data))
  L2 <- floor(2 * sqrt(T_data))
  sigmahat_data <- estimating_sigma_for_AR1(y_data, L1, L2)[[1]]
  

  ##########################################
  #Calculating gaussian quantile for T_data#
  ##########################################
  
  g_t_set           <- defining_set(T_data)
  gaussian_quantile <- quantile_function(T_data, g_t_set, test_problem, kernel_ind, alpha)
  
  
  #########################################
  #Calculating the statistic for real data#
  #########################################
  
  result                 <- statistic_function(y_data, g_t_set, kernel_ind, sigmahat_data)
  g_t_set_with_values    <- result[[1]]
  psihat_statistic_value <- result[[2]]
  
  #And now the testing itself
  if (psihat_statistic_value > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    
    #The collection of intervals where the corrected test statistic lies above the critical value (as in (2.6))
    a_t_set <- subset(g_t_set_with_values, values > gaussian_quantile, select = c(u, h, values))
    p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*T_data + 1659, 'endpoint' = (a_t_set$u + a_t_set$h)*T_data + 1659, 'values' = a_t_set$values)
    p_t_set$endpoint <- pmin(p_t_set$endpoint, 2017) 
    p_t_set <- choosing_minimal_intervals(p_t_set)

    print.xtable(xtable(subset(p_t_set, select = c(startpoint, endpoint)), digits = c(0)), type="latex", file="../Plots/mimimal_intervals.tex")

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
    
  pdffilename = paste0("../Plots/threegraphics_testing_", test_problem, "_method_", kernel_t, ".pdf")
  pdf(pdffilename, width=10, height=10, paper="special")

  par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  
  # Plotting the real data
  plot(grid_time, xlim = c(1659, 2019), y_data, type = "l", mgp=c(2,0.5,0), xaxp = c(1675, 2025,7)) 
  
  ###################################################################
  #Calculating smoothed curve for the data using Epanechnikov kernel#
  ###################################################################
  h <- c(0.01, 0.05, 0.1, 0.15, 0.2)
  plot(NA, xlim = c(1659, 2019), ylim = c(7.5, 11), yaxp  = c(8, 10, 2), xaxp = c(1675, 2025,7), mgp=c(2,0.5,0))
  for (i in 1:5){
    #This part plots kernel smoothers for different bandwidths (all on one plot).
    smoothed_curve <- mapply(epanechnikov_smoothing, grid_points, MoreArgs = list(y_data, grid_points, h[i]))
    lines(grid_time, smoothed_curve, lty = i) 
  }
  legend(1950, 9, legend=c("h = 0.01", "h = 0.05", "h = 0.10", "h = 0.15", "h = 0.2"), lty = 1:5, cex = 0.95, ncol=1)
  
  #Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
  ymaxlim = max(p_t_set$values)
  yminlim = min(p_t_set$values)
  plot(NA, xlim=c(1659,2019), xaxp = c(1675, 2025,7),  ylim = c(yminlim - 0.5, ymaxlim + 0.5), yaxp  = c(1.5, 3, 3), mgp=c(2,0.5,0))
  segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set$endpoint, p_t_set[['values']])
  abline(h = gaussian_quantile, lty = 2)
  
  #title(main = "Plots for normalized yearly temperature data for England",  outer = TRUE)
  dev.off()
}