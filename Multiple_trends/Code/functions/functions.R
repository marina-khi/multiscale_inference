add.quarters <- function(n, date_) {
  seq(date_, by = paste (3 * n, "months"), length = 2)[2]
}

#' Epanechnikov kernel function.
#' @param x A number.
#' @return 3/4(1-x^2) for |x|<=1 and 0 elsewhere.
#' @example 
#' epanechnikov_kernel(1)
epanechnikov_kernel <- function(x)
{
  if (abs(x)<=1)
  {
    result = 3/4 * (1 - x*x)
  } else {
    result = 0
  }
  return(result)
}

#' Function needed for local linear smoothing
#' @param x      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
#' @param x_vec  Vector of values for the X variables, length should be T_size
s_t_0 <- function(x, h, T_size, x_vec) {
  result = 0
  for (i in 1:T_size) {
    u = (x_vec[i] - x) / h
    result = result + epanechnikov_kernel(u)
  }
  return(result / (T_size * h));
}

#' Function needed for local linear smoothing
#' @param x      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
#' @param x_vec  Vector of values for the X variables, length should be T_size
s_t_1 <- function(x, h, T_size, x_vec) {
  result = 0
  for (i in 1:T_size) {
    u = (x_vec[i] - x) / h
    result = result + epanechnikov_kernel(u) * u
  }
  return(result / (T_size * h));
}

#' Function needed for local linear smoothing
#' @param x      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
#' @param x_vec  Vector of values for the X variables, length should be T_size
s_t_2 <- function(x, h, T_size, x_vec) {
  result = 0
  for (i in 1:T_size) {
    u = (x_vec[i] - x) / h
    result = result + epanechnikov_kernel(u) * u * u
  }
  return(result / (T_size * h));
}

#Local Linear estimator using the Epanechnikov kernel. 
local_linear_smoothing <- function(x_, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    t_len       = length(data_p)
    s_t_2_value = s_t_2(x = x_, h = bw, T_size = t_len, x_vec = grid_p)
    s_t_1_value = s_t_1(x = x_, h = bw, T_size = t_len, x_vec = grid_p)
    for (i in 1:t_len){
      u = (grid_p[i] - x_) / bw
      k = (s_t_2_value - s_t_1_value * u) * epanechnikov_kernel(u)
      result = result + k * data_p[i]
      norm = norm + k
    }
    return(result/norm)
  }
}

#Function used for simulated the three groups of time series and calculate 
#the corresponding test statistics, needed for parallel computations
repl <- function(rep, t_len_, n_ts_, sigma_, a_hat_, q_, r_, grid_, m1_, m2_){
  library(multiscale)
  
  simulated_data           <- matrix(NA, nrow = t_len_, ncol = n_ts_)
  colnames(simulated_data) <- 1:n_ts_
  
  sigmahat_vector <- c()
  for (i in 1:(floor(n_ts_ / 3))){
    simulated_data[, i] <- arima.sim(model = list(ar = a_hat_), innov = rnorm(t_len_, 0, sigma_), n = t_len_)
    simulated_data[, i] <- simulated_data[, i] - mean(simulated_data[, i])
    AR.struc            <- estimate_lrv(data = simulated_data[, i], q = q_, r_bar = r_, p = 1)
    sigma_hat_i         <- sqrt(AR.struc$lrv)
    sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
  }
  for (i in (floor(n_ts_ / 3) + 1):(floor(2 * n_ts_ / 3))){
    simulated_data[, i] <- m1_ + arima.sim(model = list(ar = a_hat_), innov = rnorm(t_len_, 0, sigma_), n = t_len_)
    simulated_data[, i] <- simulated_data[, i] - mean(simulated_data[, i])
    AR.struc            <- estimate_lrv(data = simulated_data[, i], q = q_, r_bar = r_, p = 1)
    sigma_hat_i         <- sqrt(AR.struc$lrv)
    sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
  }
  for (i in (floor(2 * n_ts_ / 3) + 1):n_ts_){
    simulated_data[, i] <- m2_ + arima.sim(model = list(ar = a_hat_), innov = rnorm(t_len_, 0, sigma_), n = t_len_)
    simulated_data[, i] <- simulated_data[, i] - mean(simulated_data[, i])
    AR.struc            <- estimate_lrv(data = simulated_data[, i], q = q_, r_bar = r_, p = 1)
    sigma_hat_i         <- sqrt(AR.struc$lrv)
    sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
  }
  psi     <- compute_statistics(data = simulated_data, sigma_vec = sigmahat_vector,
                                n_ts = n_ts_, grid = grid_)
  results <- as.vector(psi$stat_pairwise)
  return(results)
}


produce_plots_gdp <- function(results, data_i, data_j, ticks_, labels_,
                              name_i, name_j){
  filename    <- paste0("output/plots/gdp/", name_i, "_vs_", name_j, ".pdf")
  t_len       <- length(data_i)
  grid_points <- seq(1/t_len, 1, by = 1/t_len)
  
  pdf(filename, width = 5.5, height = 10.5, paper = "special")
  layout(matrix(c(1, 2, 3),ncol = 1), widths = c(2.2, 2.2, 2.2),
         heights = c(1.5, 1.5, 1.8), TRUE)
    
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 2, 0.2)) #Outer margins
  
  #if ((name_i %in% c("FRA", "USA")) & (name_j %in% c("FRA", "USA")))
  plot(data_i, ylim = c(-0.07, 0.06), type="l", col = "black", ylab = "",
       xlab = "", xaxt = "n", yaxt = "n", mgp = c(1, 0.5, 0))

  lines(data_j, col = "red")
  axis(side = 1, at = ticks_, cex.axis = 0.95, mgp = c(1, 0.5, 0), 
       labels = labels_)
  axis(side = 2, at = c(-0.05, 0, 0.05), mgp = c(1, 0.5, 0), cex.axis = 0.9)    
  title(main = "(a) adjusted GDP growth rate", font.main = 1, line = 0.5)
  legend("topright", inset = 0.01, legend=c(name_i, name_j),
         col = c("black", "red"), lty = 1, cex = 0.95, ncol = 2)
    
  par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot
    
  #Plotting the smoothed version of the time series that we have
  smoothed_i  <- mapply(local_linear_smoothing, grid_points,
                        MoreArgs = list(data_i, grid_points, bw = 0.1))
  smoothed_j  <- mapply(local_linear_smoothing, grid_points,
                        MoreArgs = list(data_j, grid_points, bw = 0.1))

  plot(smoothed_i, ylim = c(-0.027, 0.015), type = "l", col = "black", ylab = "",
       xlab = "", xaxt = "n", mgp = c(1,0.5,0), cex.axis = 0.9)
  axis(side = 1, at = ticks_, labels = labels_, cex.axis = 0.95,
       mgp = c(1, 0.5, 0))
  title(main = "(b) smoothed curves from (a)", font.main = 1, line = 0.5)
  lines(smoothed_j, col = "red")
    
  par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot
  gset    <- results$gset_with_values[[l]]
  a_t_set <- subset(gset, test == TRUE, select = c(u, h))
  if (nrow(a_t_set) > 0){
    p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * t_len,
                          'endpoint' = (a_t_set$u + a_t_set$h) * t_len,
                          'values' = 0)
    p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
      
    #Produce minimal intervals
    p_t_set2  <- compute_minimal_intervals(p_t_set)
      
    plot(NA, xlim=c(0, t_len),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab = "",
         xaxt = "n", mgp=c(2, 0.5, 0), yaxt = "n")
    axis(side = 1, at = ticks_, labels = labels_, cex.axis = 0.95,
         mgp = c(1, 0.5, 0))
    title(main = "(c) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
      #title(xlab = "quarter", line = 1.7, cex.lab = 0.9)
    segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint,
             p_t_set2$values, lwd = 2)
    segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint,
             p_t_set$values, col = "gray")
    mtext(paste0("Comparison of ", name_i, " and ", name_j), side = 3,
          line = 0, outer = TRUE, font = 1, cex = 1.2)
    dev.off()
    
    p_t_set_tex <- data.frame("from" = as.character(as.Date(sapply((p_t_set$startpoint + 0.5), add.quarters,
                                                      date_ = as.Date('01-10-1975', format = "%d-%m-%Y")))),
                              "to" = as.character(as.Date(sapply((p_t_set$endpoint - 0.5), add.quarters,
                                            date_ = as.Date('01-10-1975', format = "%d-%m-%Y")))))
    p_t_set2_tex <- data.frame("from" = as.character(as.Date(sapply((p_t_set2$startpoint + 0.5), add.quarters,
                                              date_ = as.Date('01-10-1975', format = "%d-%m-%Y")))),
                              "to" = as.character(as.Date(sapply((p_t_set2$endpoint - 0.5), add.quarters,
                                            date_ = as.Date('01-10-1975', format = "%d-%m-%Y")))))
    print.xtable(xtable(p_t_set_tex[order(p_t_set_tex$from), ], digits = c(0),
                        align = paste(replicate(3, "c"), collapse = "")),
                 file = paste0("output/tables/gdp/", name_i, "_vs_", name_j, ".tex"),
                 type = "latex", include.colnames = FALSE)
    print.xtable(xtable(p_t_set2_tex[order(p_t_set2_tex$from), ], digits = c(0),
                        align = paste(replicate(3, "c"), collapse = "")),
                 file = paste0("output/tables/gdp/", name_i, "_vs_", name_j, "_min_intervals.tex"),
                 type = "latex", include.colnames = FALSE)
    
  } else {
    #If there are no intervals where the test rejects, we produce empty plots
    plot(NA, xlim=c(0, t_len),  ylim = c(0, 1), xlab="", ylab = "", xaxt = "n",
         mgp=c(2,0.5,0), yaxt = "n")
    axis(side = 1, at = ticks_, labels = labels_, cex.axis = 0.95,
         mgp = c(1, 0.5, 0))
    title(main = "(c) (minimal) intervals produced by our test", font.main = 1,
          line = 0.5)
    mtext(paste0("Comparison of ", name_i, " and ", name_j), side = 3,
          line = 0, outer = TRUE, font = 1, cex = 1.2)
    dev.off()
  }
}


produce_plots_talk <- function(results, l, data_i, data_j, at_, labels_, dates_,
                               dates, name_i, name_j, filename){
  filename    <- paste0("output/plots/talk/VOC/", name_i, "_vs_", name_j, ".pdf")
  grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
  
  Tlen <- length(data_i)
  gset <- results$gset_with_values[[l]]
  
  pdf(filename, width = 5, height = 6.5, paper="special")
  layout(matrix(c(1, 2), ncol=1), widths=c(2.4, 2.4),
         heights=c(1.5, 1.8), TRUE)
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins

  plot(x = dates_, y = data_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
       col = "#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
  lines(x = dates_, y = data_j, col="#604c38")
  title(main = "(a) adjusted logarithm of house prices", font.main = 1, line = 0.5)
  legend("topright", inset = 0.02, legend=c(name_i, name_j),
         col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
  
  par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot
  
  a_t_set <- subset(gset, test == TRUE, select = c(u, h))
  if (nrow(a_t_set) > 0){
    p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * Tlen + 0.5,
                          'endpoint' = (a_t_set$u + a_t_set$h) * Tlen - 0.5, 'values' = 0)
    p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
    
    #Produce minimal intervals
    p_t_set2  <- compute_minimal_intervals(p_t_set)
    
    plot(NA, xlim=c(0, Tlen),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab="", mgp=c(2, 0.5, 0), yaxt = "n")
    title(main = "(b) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
    segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint, p_t_set2$values, lwd = 2)
    segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint, p_t_set$values, col = "gray")
  } else {
    #If there are no intervals where the test rejects, we produce empty plots
    plot(NA, xlim=c(0, Tlen),  ylim = c(0, 1), xlab="", ylab = "", mgp=c(2,0.5,0), yaxt = "n")
    title(main = "(b) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
  }
  dev.off()
}


produce_plots_hp <- function(results, data_i, data_j,
                             at_, labels_, name_i, name_j, l){
  filename    <- paste0("output/plots/hp/", name_i, "_vs_", name_j, ".pdf")
  t_len       <- length(data_i)
  grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
  
  pdf(filename, width = 5.5, height = 10.5, paper="special")
  layout(matrix(c(1, 2, 3),ncol=1), widths=c(2.2, 2.2, 2.2),
         heights = c(1.5, 1.5, 1.8), TRUE)
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 2, 0.2)) #Outer margins
  
  if ((name_i %in% c("AUS", "NLD")) & (name_j %in% c("AUS", "NLD"))) {
    shift <- 0.4
  } else {
    shift <- 0
  }
  
  plot(data_i, ylim = c(-0.9, 1.5),
       #ylim = c(min(data_i, data_j), max(data_i, data_j) + shift),
       type = "l", col = "black", ylab = "", xlab = "", xaxt = "n",
       mgp = c(1, 0.5, 0))
  lines(data_j, col="red")
  axis(side = 1, at = at_, labels = labels_,
       cex.axis = 0.95, mgp=c(1, 0.5, 0))
  
  title(main = "(a) adjusted log of housing prices", font.main = 1, line = 0.5)
  legend("topright", inset = 0.02, legend = c(name_i, name_j),
         col = c("black", "red"), lty = 1, cex = 0.95, ncol = 2)
  
  par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot
  
  #Plotting the smoothed version of the time series that we have
  smoothed_i  <- mapply(local_linear_smoothing, grid_points,
                        MoreArgs = list(data_i, grid_points, bw = 7/t_len))
  smoothed_j  <- mapply(local_linear_smoothing, grid_points,
                        MoreArgs = list(data_j, grid_points, bw = 7/t_len))
  
  plot(smoothed_i, ylim = c(-0.9, 1.5),
       #ylim = c(min(data_i, data_j), max(data_i, data_j)),
       type = "l",
       col="black", ylab = "", xlab = "", xaxt = "n", mgp = c(1,0.5,0))
  axis(side = 1, at = at_, labels = labels_, cex.axis = 0.95,
       mgp = c(1, 0.5, 0))
  title(main = "(b) smoothed curves from (a)", font.main = 1, line = 0.5)
  lines(smoothed_j, col="red")
  
  par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot
  gset    <- results$gset_with_values[[l]]
  a_t_set <- subset(gset, test == TRUE, select = c(u, h))
  if (nrow(a_t_set) > 0){
    p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * t_len + 0.5,
                          'endpoint' = (a_t_set$u + a_t_set$h) * t_len - 0.5,
                          'values' = 0)
    p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
    
    #Produce minimal intervals
    p_t_set2  <- compute_minimal_intervals(p_t_set)
    
    plot(NA, xlim=c(0, t_len),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab = "",
         xaxt = "n", mgp = c(2, 0.5, 0), yaxt = "n")
    axis(side = 1, at = at_, labels = labels_, cex.axis = 0.95,
         mgp = c(1, 0.5, 0))
    title(main = "(c) (minimal) intervals produced by our test", font.main = 1,
          line = 0.5)
    segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint,
             p_t_set2$values, lwd = 2)
    segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint,
             p_t_set$values, col = "gray")
    
    p_t_set_tex <- data.frame("from" = as.character((p_t_set$startpoint - 0.5) + dates[1] - 1),
                              "to" = as.character((p_t_set$endpoint + 0.5) + dates[1] - 1))
    p_t_set2_tex <- data.frame("from" = as.character((p_t_set2$startpoint - 0.5) + dates[1] - 1),
                              "to" = as.character((p_t_set2$endpoint + 0.5) + dates[1] - 1))
    print.xtable(xtable(p_t_set_tex[order(p_t_set_tex$from), ], digits = c(0),
                        align = paste(replicate(3, "c"), collapse = "")),
                 file = paste0("output/tables/hp/", name_i, "_vs_", name_j, ".tex"),
                 type = "latex", include.colnames = FALSE)
    print.xtable(xtable(p_t_set2_tex[order(p_t_set2_tex$from), ], digits = c(0),
                        align = paste(replicate(3, "c"), collapse = "")),
                 file = paste0("output/tables/hp/", name_i, "_vs_", name_j, "_min_intervals.tex"),
                 type = "latex", include.colnames = FALSE)
    mtext(paste0("Comparison of ", name_i, " and ", name_j), side = 3, line = 0,
          outer = TRUE, font = 1, cex = 1.2)
    dev.off()
  } else {
    #If there are no intervals where the test rejects, we produce empty plots
    plot(NA, xlim = c(0, t_len),  ylim = c(0, 1), xlab = "", ylab = "",
         xaxt = "n", mgp = c(2,0.5,0), yaxt = "n")
    axis(side = 1, at = at_, labels = labels_, cex.axis = 0.95,
         mgp = c(1, 0.5, 0))
    title(main = "(c) (minimal) intervals produced by our test", font.main = 1,
          line = 0.5)
    mtext(paste0("Comparison of ", name_i, " and ", name_j), side = 3, line = 0,
          outer = TRUE, font = 1, cex = 1.2)
    dev.off()
  }
}


#Create a matrix (for size and power table for example) and write them in the tex file
output_matrix <- function(matrix_, filename){
  addtorow     <- list()
  addtorow$pos <- list(0, 0)
  addtorow$command <- c("& \\multicolumn{3}{c}{nominal size $\\alpha$} \\\\\n",
                        "$T$ & 0.01 & 0.05 & 0.1 \\\\\n") 
  print.xtable(xtable(matrix_, digits = c(3), align = "cccc"), type = "latex",
               file = filename, add.to.row = addtorow, include.colnames = FALSE)
}

#Function that simulates the covariates as AR(1), the error terms as AR(1)
#the time series as y = beta_ %*% covariates + m_matrix_ + errors,
#estimates the parameters, and then computes the test statistics
repl2 <- function(rep, t_len_, n_ts_, grid_, gaussian_sim = FALSE,
                  a_ = 0, sigma_ = 1, beta_ = 0, a_x_ = 0, sigma_x_ = 1,
                  m_matrix_ = NULL, q_ = 25, r_ = 10){
  library(multiscale)
  library(dplyr)
  
  if (gaussian_sim){
    z_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    z_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    sigma_vector  <- rep(1, n_ts_)
    
    for (i in 1:n_ts_){
      z_matrix[, i]      <- rnorm(t_len_, 0, sigma_)
      z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    }
    
    psi <- compute_statistics(data = z_augm_matrix,
                              sigma_vec = sigma_vector,
                              n_ts = n_ts_, grid = grid_)
  } else {
    y_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    y_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    
    sigmahat_vector <- c()
    
    if (beta_ != 0){
      x_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      for (i in 1:n_ts_){
        error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                       innov = rnorm(t_len_, 0, sigma_),
                                       n = t_len_)
        x_matrix[, i]     <- arima.sim(model = list(ar = a_x_),
                                       innov = rnorm(t_len_, 0, sigma_x_),
                                       n = t_len_)
        y_matrix[, i]     <- m_matrix_[, i] + beta_ * x_matrix[, i] + error_matrix[, i]
        
        #First differences
        x_diff    <- x_matrix[, i]- dplyr::lag(x_matrix[, i], n = 1, default = NA)
        y_diff    <- y_matrix[, i]- dplyr::lag(y_matrix[, i], n = 1, default = NA)
        
        #Estimating beta
        x_diff_tmp <- as.matrix(x_diff)[-1, ]
        y_diff_tmp <- as.matrix(y_diff)[-1, ]
        
        beta_hat_tmp  <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
        alpha_hat_tmp <- mean(y_matrix[, i] - x_matrix[, i] * as.vector(beta_hat_tmp))

        y_augm_matrix[, i] <- y_matrix[, i] - x_matrix[, i] * as.vector(beta_hat_tmp) - alpha_hat_tmp
        AR.struc           <- estimate_lrv(data = y_augm_matrix[, i], q = q_,
                                           r_bar = r_, p = 1)
        sigma_hat_i        <- sqrt(AR.struc$lrv)
        sigmahat_vector    <- c(sigmahat_vector, sigma_hat_i) 
      }
    } else {
      for (i in 1:n_ts_){
        error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                       innov = rnorm(t_len_, 0, sigma_),
                                       n = t_len_)
        y_matrix[, i]     <- m_matrix_[, i] + error_matrix[, i]
        
        #Estimating alpha_i
        alpha_hat_tmp     <- mean(y_matrix[, i])

        y_augm_matrix[, i]  <- y_matrix[, i] - alpha_hat_tmp
        AR.struc            <- estimate_lrv(data = y_augm_matrix[, i], q = q_,
                                            r_bar = r_, p = 1)
        sigma_hat_i         <- sqrt(AR.struc$lrv)
        sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)   
      }     
    }
    psi <- compute_statistics(data = y_augm_matrix,
                              sigma_vec = sigmahat_vector,
                              n_ts = n_ts_, grid = grid_)    
  }
  results <- as.vector(psi$stat_pairwise)
  return(results)
}