#Calculate autocovariance function for AR(1) model: \varepsilon_t = a_1 \varepsilon_{t-1} + \eta_t based on true coefficients of the model
autocovariance_function_AR1 <- function(k, a_1, sigma_eta){
  #if (k%%1==0)
  #{
    result = sigma_eta * sigma_eta * a_1^(abs(k)) / (1 - a_1 * a_1)
  #} else {
  #  print('Check the input: k is not integer')
  #}
  return(result)
}

#Function that compares the minimal intervals for both SiZer and our method.
#For each run of this function the result is a pdf file with three plots:
#1. Possible time series together with underlying time trend
#2. Union of minimal intervals for our method for N_rep number of simulations
#3. Union of minimal intervals for SiZer for N_rep number of simulations
#The output file is hardcoded!
plotting_many_minimal_intervals <- function(trend_height, trend_width, T_size, SiZer_matrix, N_rep, kernel_ind, sigmahat, gaussian_quantile, a_1, sigma_eta){
  matrix_our_results   <- data.frame('startpoint' = SiZer_matrix$u - SiZer_matrix$h, 'endpoint' = SiZer_matrix$u + SiZer_matrix$h)
  matrix_their_results <- data.frame('startpoint' = SiZer_matrix$u - SiZer_matrix$h, 'endpoint' = SiZer_matrix$u + SiZer_matrix$h)
  
  biweight_trend  <- numeric(T_size)
  for (i in 1:T_size) {biweight_trend[i] = trend_height * biweight_kernel(trend_width *(i/T_size - 0.5))}
  
  for (col in 1:N_rep){
    y_data_ar_1_with_trend <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta)) + biweight_trend
    
    g_t_set     <- psihat_statistic_ll(y_data_ar_1_with_trend, SiZer_matrix, kernel_ind, sigmahat)[[1]]
    
    results_our   <- c()
    results_their <- c()
    
    for (row in 1:nrow(g_t_set)){
      i                = g_t_set[row, 'u']
      h                = g_t_set[row, 'h']
      q_h              = g_t_set[row, 'q_h']
      sd_m_hat_prime   = g_t_set[row, 'sd']
      XtWX_inverse_XtW = g_t_set$XtWX_inv_XtW[[row]]
      
      if (!is.null(XtWX_inverse_XtW)){
        m_hat_prime <- (XtWX_inverse_XtW %*% y_data_ar_1_with_trend)[2]
        if (m_hat_prime - q_h * sd_m_hat_prime > 0){
          results_their = c(results_their, 1)
        } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
          results_their = c(results_their, -1)
        } else {
          results_their = c(results_their, 0)
        }
      } else {results_their = c(results_their, 0)}
      
      if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
        results_our = c(results_our, 1)
      } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
        results_our = c(results_our, -1)
      } else {
        results_our = c(results_our, 0)
      }
    }
    matrix_our_results <- cbind(matrix_our_results, results_our)
    matrix_their_results <- cbind(matrix_their_results, results_their)
  }
  
  grid_points <- seq(from = 1/T_size, to = 1, length.out = T_size) #grid points for estimating
  
  pdffilename = paste0("Paper/Plots/min_int_with_T_", T_size, "_a1_", a_1*100, ".pdf")
  pdf(pdffilename, width=8, height=10, paper="special")
  
  par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(1.5, 0.5, 0, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  
  plot(grid_points, y_data_ar_1_with_trend, ylim = c(min(y_data_ar_1_with_trend) - 0.2, max(y_data_ar_1_with_trend)+0.2), type = "l")
  lines(grid_points, biweight_trend)

  par(mar = c(1.5, 0.5, 2, 0)) #Margins for each plot
  
  plot(NA, xlim=c(0,1), ylim = c(-1, N_rep +1), main = "Our test")
  for (col in 3:(N_rep+2)){
    a_t_set <- subset(matrix_our_results, matrix_our_results[,col] != 0, select = c(startpoint, endpoint, col))
    colnames(a_t_set) <- c('startpoint', 'endpoint', 'values')
    p_t_set <- choosing_minimal_intervals(a_t_set)
    if (nrow(p_t_set) > 0) {segments(pmax(0, p_t_set$startpoint), col-2, pmin(1, p_t_set$endpoint), col-2)}
    #cat("Number of repetition:", col, "\n")
  }
  
  plot(NA, xlim=c(0,1), ylim = c(-1, N_rep +1), main = "SiZer")
  for (col in 3:(N_rep+2)){
    a_t_set <- subset(matrix_their_results, matrix_their_results[,col] != 0, select = c(startpoint, endpoint, col))
    colnames(a_t_set) <- c('startpoint', 'endpoint', 'values')
    p_t_set <- choosing_minimal_intervals(a_t_set)
    if (nrow(p_t_set) > 0){segments(pmax(0, p_t_set$startpoint), col-2, pmin(1, p_t_set$endpoint), col-2)}
  }
  dev.off()
}




plotting_MSE_graphs <- function(data1, data2, data3, pdfname, margin_, legend_position,
                                ylab_, legend_, title_, different_a, zoomed){
  data_min <- min(c(data1, data2, data3))
  data_max <- max(c(data1, data2, data3))
  
  pdf(pdfname, width=5.5, height=4.16, paper="special")
  
  par(mar = c(4, 4, margin_, 0)) #Margins for each plot
  par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins 
  
  if (zoomed == 'yes'){
    plot(data1, type="l", lty=1, xaxt='n', ylim = c(0.000, 0.020), ylab = ylab_, xlab=expression(a[1]))
  } else {
    plot(data1,ylim=c(data_min,data_max),type="l",lty=1,xaxt='n', ylab = ylab_, xlab=expression(a[1]))
  }
  points(data1,pch=19)
  lines(data2,lty="dashed",lwd=1.5)
  points(data2,pch=19)
  lines(data3,lty="dotted",lwd=1.5)
  points(data3,pch=19)
  axis(1, at=1:length(different_a), labels=different_a)
  legend(legend_position, cex = 1, bty = "n", legend = legend_, lty=c("solid","dashed","dotted"),lwd=1.5,y.intersp=1.25)
  title(main = title_, line = 1.5)
  dev.off()
}

PlotHistograms <- function(data1, data2, data3, star_value, pdfname, text1, text2, text3, cut.at.end = FALSE){
  smallest_value <- min(c(data1, data2, data3))
  biggest_value  <- max(c(data1, data2, data3))
  step_len       <- (biggest_value-smallest_value)/50
  
  breaks_grid <- seq(smallest_value, biggest_value, by = step_len)
  breaks_grid[length(breaks_grid)] <- biggest_value
  
  if (cut.at.end){
    if (smallest_value < -1.4){
      steps <- ceiling(50*(biggest_value-smallest_value)/(biggest_value+1.4))
      step_len <- (biggest_value-smallest_value)/steps
    }
    smallest_value <- max(c(-1.4, smallest_value))
  }
  
  hist1 <- hist(data1, breaks = breaks_grid, plot = FALSE)
  hist2 <- hist(data2, breaks = breaks_grid, plot = FALSE)
  hist3 <- hist(data3, breaks = breaks_grid, plot = FALSE)
  
  highestCount <- max(hist1$counts, hist2$counts, hist3$counts)
  
  pdf(pdfname, width=8, height=2.9, paper="special")
  par(mfrow = c(1,3))
  par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.5, 0.2)) #Outer margins
  
  hist(data1, main = NULL, breaks = breaks_grid, freq=TRUE, xlim=c(smallest_value,biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab = 1.1)
  mtext(side=1,text= text1,line=2.75)
  segments(x0=star_value,y0=0,x1=star_value,y1=highestCount,col="red",lwd=1.5)
  
  hist(data2, main = NULL, breaks = breaks_grid, freq=TRUE, xlim=c(smallest_value,biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab = 1.1)
  mtext(side=1,text= text2,line=2.75)
  segments(x0=star_value,y0=0,x1=star_value,y1=highestCount,col="red",lwd=1.5)
  
  hist(data3, main = NULL, breaks = breaks_grid, freq=TRUE, xlim=c(smallest_value,biggest_value), ylim=c(0,highestCount), xlab = "", mgp=c(2,0.5,0), cex.lab = 1.1)
  mtext(side=1,text= text3,line=2.75)
  segments(x0=star_value,y0=0,x1=star_value,y1=highestCount,col="red",lwd=1.5)
  dev.off()
}

plot.SiZer <- function(results, different_i, different_h, ylab=expression(log[10](h)), 
                       colorlist=c('red', 'purple', 'blue', 'grey'), title = "Results", ...){
  temp <- factor(results);
  final.colorlist <- NULL;
  if( is.element( '-1', levels(temp) ) )
    final.colorlist <- c(final.colorlist, colorlist[1]);
  if( is.element( '0', levels(temp) ) )
    final.colorlist <- c(final.colorlist, colorlist[2]);
  if( is.element( '1', levels(temp) ) )
    final.colorlist <- c(final.colorlist, colorlist[3]);
  if( is.element( '2', levels(temp) ) )
    final.colorlist <- c(final.colorlist, colorlist[4]);
  
  # Convert the slopes to a factor list to match up with final.colorlist
  temp <- matrix( as.integer(factor(results)), nrow=length(different_h), byrow = TRUE) 
  
  graphics::image( different_i, log(different_h,10), t(temp), 
                   col=final.colorlist, ylab=ylab, main = title, ...)
}

calculating_sigma_eta <- function(y_data, coefficients, p){
  y_diff <- c(0, diff(y_data))
  T_size <- length(y_diff)
  
  res    <- rep(0, T_size)
  
  for(i in (p+1):T_size){
    res[i] <- y_diff[i] - sum(coefficients * y_diff[(i-1):(i-p)])
  }
  
  var_eta <- mean(res^2)/2
  return(sqrt(var_eta))
}
