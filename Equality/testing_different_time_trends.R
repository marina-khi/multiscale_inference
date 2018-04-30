testing_different_time_trends <- function(N_ts, y_data, month_column, alpha, kernel_method, sigmahat_vector_2){

  T_data <- nrow(y_data)

  ######################################################################
  #Calculating smoothed curve for the data using local linear estimator#
  ######################################################################
  
  grid_points <- seq(from = 1/T_data, to = 1, length.out = T_data) #grid points for estimating
  
  pdf("../Plots/stations_data.pdf", width=10, height=10, paper="special")
  par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
  par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins
  
  plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(-3.0, 2.5), yaxp  = c(-2.0, 2.0, 4), xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
  for (column in TemperatureColumns){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, 0.05))
    lines(grid_points, smoothed_curve)
  }
  axis(1, at = grid_points[seq(1, 380, by = 20)], labels = monthly_temp$date[seq(1, 380, by = 20)])
  legend(0, 1.5, legend=c("h = 0.05"), lty = 1, cex = 0.95, ncol=1)
  
  plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(-2.5, 2.5), yaxp  = c(-2.0, 2.0, 4),xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
  for (column in TemperatureColumns){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, 0.1))
    lines(grid_points, smoothed_curve)
  }
  axis(1, at = grid_points[seq(1, 380, by = 20)], labels = monthly_temp$date[seq(1, 380, by = 20)])
  legend(0, 1.5, legend=c("h = 0.10"), lty = 1, cex = 0.95, ncol=1)
  
  
  plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(-2.5, 2.5), yaxp  = c(-2.0, 2.0, 4),xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
  for (column in TemperatureColumns){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(monthly_temp[[column]], grid_points, 0.15))
    lines(grid_points, smoothed_curve)
  }
  axis(1, at = grid_points[seq(1, 380, by = 20)], labels = monthly_temp$date[seq(1, 380, by = 20)])
  legend(0, 1.5, legend=c("h = 0.15"), lty = 1, cex = 0.95, ncol=1)
  
  dev.off()
  
  
  ##########################################
  #Calculating gaussian quantile for T_data#
  ##########################################
  
  g_t_set <- creating_g_set(T_data, kernel_method)
  
  filename = paste("distribution/distr_for_application_T_", T_data, "_and_N_ts_", N_ts, "_and_method_", kernel_method, ".RData", sep = "")
  if(!file.exists(filename)) {
    gaussian_statistic_distribution <- replicate(1000, {
      z_matrix <- matrix(0, nrow = T_data, ncol = N_ts + 1)
      for (i in 1:N_ts){
        z_matrix[, i] <- rnorm(T_data, 0, sqrt(sigmahat_vector_2[i]))
      }
      z_matrix <- cbind(z_matrix, month_column) #Adding auxiliary month column to "deseasonalize" simulated data
      z_matrix[1:N_ts] <- lapply(z_matrix[1:N_ts], function(x) x - ave(x, z_matrix[N_ts + 1], FUN=mean))
      z_matrix <- z_matrix[-(N_ts + 1)] #Dropping auxiliary month
      psistar_statistic(z_matrix, N_ts, g_t_set, sigmahat_vector_2, kernel_method)
    })
    save(gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)  
  cat("Gaussian quantile is", gaussian_quantile, "\n")
  
  #########################################
  #Calculating the statistic for real data#
  #########################################
  
  statistic <- psihat_statistic(y_data, N_ts, g_t_set, sigmahat_vector_2, kernel_method)

  #And now the testing itself
  if (statistic[[2]] > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", statistic[[2]],
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", statistic[[2]],
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  }
  
  return(list(statistic[[2]], gaussian_quantile))
}