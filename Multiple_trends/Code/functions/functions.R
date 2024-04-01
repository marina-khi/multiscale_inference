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

#Simulating a bump function
bump  <- function(u){
  u.lower <- 0.3
  u.upper <- 0.7
  arg <- (u - 0.5)/(u.upper - 0.5)
  return(as.double(u >= u.lower & u <= u.upper) * (1 - arg^2)^2)
}

#Simulating a kind of a bump function for clustering
b_function <- function(u, x_0, h){
  arg <- (u - x_0)/h
  return(as.double((abs(arg) <= 1) * (1 - arg^2)^2))
}


#Function used for simulated the three groups of time series and calculate 
#the corresponding test statistics, needed for parallel computations
repl <- function(rep, t_len_, n_ts_, sigma_, a_hat_, q_, r_, grid_, m1_, m2_){
  library(MSinference)
  
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

#Function used for simulated the three groups of time series and calculate 
#the corresponding test statistics, needed for parallel computations
repl_clustering_revision <- function(rep, t_len_, n_ts_, grid_,
                                     m1_ = NULL, m2_ = NULL, 
                                     a_hat_ = 0, sigma_ = 1,
                                     q_ = 25, r_ = 10, h_ = 0.05,
                                     gaussian_sim = FALSE){
  library(MSinference)

  if (gaussian_sim){
    z_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    z_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    sigma_vector  <- rep(sigma_, n_ts_)
    
    for (i in 1:n_ts_){
      z_matrix[, i]      <- rnorm(t_len_, 0, sigma_)
      z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    }
    
    psi <- compute_statistics(data = z_augm_matrix,
                              sigma_vec = sigma_vector,
                              n_ts = n_ts_, grid = grid_)
    results <- c(as.vector(psi$stat_pairwise))
  } else {
    simulated_data           <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    #colnames(simulated_data) <- 1:n_ts_
  
    sigmahat_vector <- c()
    for (i in 1:(floor(n_ts_ / 3))){
      simulated_data[, i] <- arima.sim(model = list(ar = a_hat_), innov = rnorm(t_len_, 0, sigma_), n = t_len_)
      simulated_data[, i] <- simulated_data[, i] - mean(simulated_data[, i])
      #AR.struc            <- estimate_lrv(data = simulated_data[, i], q = q_, r_bar = r_, p = 1)
      #sigma_hat_i         <- sqrt(AR.struc$lrv)
      #sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
    }
    for (i in (floor(n_ts_ / 3) + 1):(floor(2 * n_ts_ / 3))){
      simulated_data[, i] <- m1_ + arima.sim(model = list(ar = a_hat_), innov = rnorm(t_len_, 0, sigma_), n = t_len_)
      simulated_data[, i] <- simulated_data[, i] - mean(simulated_data[, i])
      #AR.struc            <- estimate_lrv(data = simulated_data[, i], q = q_, r_bar = r_, p = 1)
      #sigma_hat_i         <- sqrt(AR.struc$lrv)
      #sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
    }
    for (i in (floor(2 * n_ts_ / 3) + 1):n_ts_){
      simulated_data[, i] <- m2_ + arima.sim(model = list(ar = a_hat_), innov = rnorm(t_len_, 0, sigma_), n = t_len_)
      simulated_data[, i] <- simulated_data[, i] - mean(simulated_data[, i])
      #AR.struc            <- estimate_lrv(data = simulated_data[, i], q = q_, r_bar = r_, p = 1)
      #sigma_hat_i         <- sqrt(AR.struc$lrv)
      #sigmahat_vector     <- c(sigmahat_vector, sigma_hat_i)
    }
    
    sigmahat_vector <- rep(sqrt(sigma_^2/((1 - a_hat_)^2)), n_ts_)
    psi <- compute_statistics(data = simulated_data, sigma_vec = sigmahat_vector,
                              n_ts = n_ts_, grid = grid_)
    
    grid_points    <- seq(1/t_len_, 1, by = 1/t_len_)
    grid_points_2  <- seq(1/(2 * t_len_), 1, by = 1/t_len_)
    smoothed_data  <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    smoothed_data2 <- matrix(NA, nrow = t_len_, ncol = n_ts_)

    for (i in 1:n_ts_){
      smoothed_data[, i] <- mapply(local_linear_smoothing, grid_points_2, MoreArgs = list(simulated_data[, i], grid_points, h_))
      smoothed_data2[, i] <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(simulated_data[, i], grid_points, h_))
    }
    
    #Calculating the benchmark model
    benchmark_results <- matrix(0, nrow = n_ts_, ncol = n_ts_)
    for (i in 1:(n_ts_ - 1)){
      for (j in (i + 1):n_ts_){
        benchmark_results[i, j] <- sum((smoothed_data[, i] - smoothed_data[, j])^2) * (1/t_len_)
      }
    }
    
    #Calculating the benchmark model
    benchmark_results2 <- matrix(0, nrow = n_ts_, ncol = n_ts_)
    for (i in 1:(n_ts_ - 1)){
      for (j in (i + 1):n_ts_){
        benchmark_results2[i, j] <- max(abs(smoothed_data2[, i] - smoothed_data2[, j]))
      }
    }
    results <- c(as.vector(psi$stat_pairwise), as.vector(benchmark_results), as.vector(benchmark_results2))
  }
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

produce_plots_gdp2 <- function(results, data_i, data_j, ticks_, labels_,
                              name_i, name_j){
  filename    <- paste0("output/revision/gdp/", name_i, "_vs_", name_j, ".pdf")
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

produce_plots_hp2 <- function(results, data_i, data_j,
                             at_, labels_, name_i, name_j, l){
  filename    <- paste0("output/revision/hp/", name_i, "_vs_", name_j, ".pdf")
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

output_matrix2 <- function(matrix_, filename){
  addtorow     <- list()
  addtorow$pos <- list(0, 0)
  addtorow$command <- c("& \\multicolumn{6}{c}{nominal size $\\alpha$} \\\\\n",
                        "$T$ & 0.01 & 0.05 & 0.1 & 0.01 & 0.05 & 0.1 \\\\\n") 
  print.xtable(xtable(matrix_, digits = c(3), align = "ccccccc"), type = "latex",
               file = filename, add.to.row = addtorow, include.colnames = FALSE)
}


#Function that simulates 3 covariates as VAR(3) process with the given
#coefficients (a_x_mat_ and sigma_x_mat_),
#the error terms as AR(1) also with the given coefficient a_ and sigma_,
#the fixed effect term alpha_ as a normally distributed random vector,
#N(0, c_ * Sigma_a_mat_), the time series as
#y = alpha_ + beta_ %*% covariates + m_matrix_ + errors,
#estimates the parameters, and then computes the test statistics
repl_revision <- function(rep_, n_ts_, t_len_, grid_, a_ = 0, sigma_ = 1,
                          beta_ = NULL,
                          a_x_vec_ = c(0, 0, 0), phi_ = 0,
                          rho_ = 0, m_matrix_ = NULL,
                          q_ = 25, r_ = 10,
                          gaussian_sim = FALSE){

  library(MSinference)
  library(dplyr)
  
  if (gaussian_sim){
    z_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    z_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    sigma_vector  <- rep(sigma_, n_ts_)

    for (i in 1:n_ts_){
      z_matrix[, i]      <- rnorm(t_len_, 0, sigma_)
      z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    }

    psi <- compute_statistics(data = z_augm_matrix,
                              sigma_vec = sigma_vector,
                              n_ts = n_ts_, grid = grid_)
    results <- c(as.vector(psi$stat_pairwise))
  } else {
    
    if (is.null(m_matrix_)){
      m_matrix_ <- matrix(0, nrow = t_len_, ncol = n_ts_)
    }
    
    library(mvtnorm)
    y_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    y_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)

    big_sigma_matrix       <- matrix(rho_, nrow = n_ts_, ncol = n_ts_)
    diag(big_sigma_matrix) <- 1
    alpha_vec              <- rmvnorm(1, mean = rep(0, n_ts_), sigma = big_sigma_matrix)

    sigmahat_vector <- c()

    if (!is.null(beta_)){
      phi_matrix       <- matrix(phi_, nrow = 3, ncol = 3)
      diag(phi_matrix) <- 1
      
      for (i in 1:n_ts_){
        error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                       innov = rnorm(t_len_, 0, sigma_),
                                       n = t_len_)
        
        nu       <- rmvnorm(t_len_ + 10, mean = c(0, 0, 0), sigma = phi_matrix)
        a_matrix <- diag(a_x_vec_)
        x_matrix <- matrix(0, 3, t_len_ + 10)

        for (t in 2:(t_len_ + 10)){
          x_matrix[, t] <- a_matrix %*% x_matrix[, t - 1] + nu[t, ]
        }
        x_matrix <- t(x_matrix[, -(1:10)])
        
        y_matrix[, i] <- alpha_vec[i] + m_matrix_[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
        
        #First differences
        x_diff_1  <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
        x_diff_2  <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
        x_diff_3  <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)
        y_diff    <- y_matrix[, i] - dplyr::lag(y_matrix[, i], n = 1, default = NA)
        
        #Estimating beta
        x_diff_tmp <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]
        y_diff_tmp <- as.matrix(y_diff)[-1, ]
        
        beta_hat_tmp       <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
        alpha_hat_tmp      <- mean(y_matrix[, i] - x_matrix %*% as.vector(beta_hat_tmp))
        
        y_augm_matrix[, i] <- y_matrix[, i] - x_matrix %*% as.vector(beta_hat_tmp) - alpha_hat_tmp
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
        
        y_matrix[, i]     <- alpha_vec[i] + m_matrix_[, i] + error_matrix[, i]
        
        #Estimating the fixed effects
        alpha_hat_tmp      <- mean(y_matrix[, i])
        y_augm_matrix[, i] <- y_matrix[, i] - alpha_hat_tmp
        AR.struc           <- estimate_lrv(data = y_augm_matrix[, i], q = q_,
                                           r_bar = r_, p = 1)
        sigma_hat_i        <- sqrt(AR.struc$lrv)
        sigmahat_vector    <- c(sigmahat_vector, sigma_hat_i) 
      }     
    }
    psi <- compute_statistics(data = y_augm_matrix,
                              sigma_vec = sigmahat_vector,
                              n_ts = n_ts_, grid = grid_)    
    results <- c(as.vector(psi$stat_pairwise))
  }
  return(results)
}


# pdf(paste0("output/revision/bump_function.pdf"),
#     width = 12, height = 8, paper="special")
# par(mfrow = c(2, 2))
# par(mar = c(4, 3, 0.5, 0)) #Margins for each plot
# par(oma = c(0.5, 0.5, 0.5, 0.2)) #Outer margins
# 
# for (b in different_b){
#   errors <- arima.sim(model = list(ar = a),
#                       innov = rnorm(t_len, 0, sigma),
#                       n = t_len)
#   plot(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
#        y = (bump((1:t_len)/t_len) * b), ylim = c(-1, 1.6),
#        xlab = "", ylab = "", main = NULL,
#        type = 'l', cex = 0.8)
#   lines(x = seq(from = 1 / t_len, to = 1, by = 1 / t_len),
#         y = (bump((1:t_len)/t_len) * b) + errors, type = "l",
#         col = "red")
#   mtext(side = 1, text = paste0("b = ", b), line = 2.3, cex = 1)
# }
# dev.off()


#################
#SIZER FUNCTIONS#
#################


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


ESS.star <- function(u.grid, h.grid, T, autocov)
  
{ # compute the effective sample size ESS.star
  #
  # Arguments:
  # u.grid       grid of locations
  # h.grid       grid of bandwidths
  # T            time series length
  # autocov      vector of (estimated) error autocovariances 
  # 
  # Outputs:
  # ess          matrix with length(u.grid) columns and length(h.grid) rows
  #              specifying ESS for each point (u,h)
  # ess.star     matrix with length(u.grid) columns and length(h.grid) rows
  #              specifying ESS.star for each point (u,h)
  # deletions    vector of length length(u.grid)*length(h.grid) with NA elements
  #              in places where ESS.star<5
  
  N.u <- length(u.grid)
  N.h <- length(h.grid)
  
  ess      <- matrix(NA,ncol=N.u,nrow=N.h)
  ess.star <- matrix(NA,ncol=N.u,nrow=N.h)
  
  for(i in 1:N.h){
    bw <- h.grid[i]
    
    pos.int <- 1:N.u     
    temp    <- ( u.grid - bw >= 0 & u.grid + bw <= 1 )
    if(sum(temp) > 0){
      pos.int <- pos.int[temp]
      u       <- u.grid[pos.int[1]]
      arg     <- ((1:T)/T - u)/bw
      ess[i,pos.int] <- sum(epanechnikov_kernel(arg)) / 0.75
    }
    
    pos.bnd <- 1:N.u
    temp    <- ( u.grid - bw < 0 | u.grid + bw > 1 )
    if(sum(temp) > 0){
      pos.bnd <- pos.bnd[temp]
      for(j in 1:length(pos.bnd)){
        u   <- u.grid[pos.bnd[j]]
        arg <- ((1:T)/T - u)/bw
        ess[i,pos.bnd[j]] <- sum(epanechnikov_kernel(arg)) / 0.75
      }
    }
  }
  
  cov.wghts <- 1 - (1:(T-1))/T
  V.bar     <- autocov[1]/T + (2/T) * sum(cov.wghts * autocov[2:T])
  T.star    <- autocov[1] / V.bar
  ess.star  <- (T.star/T) * ess
  
  deletions <- 1:(N.u*N.h)
  temp      <- as.vector(t(ess.star))
  deletions[temp < 5] <- NA
  
  return(list(ess = ess, ess.star=ess.star,del=deletions))
}

SiZer_weights <- function(t_len, grid){
  # calculate the kernel weights for SiZer 
  #
  # Arguments:
  # t_len        sample size 
  # grid         grid of location-bandwidth points as produced by the function 'grid_construction',
  #              list with the element 'gset' (and possibly others)
  #
  # Outputs: 
  # weights      matrix of kernel weights
  #              w_1(u_1,h_1), ..., w_T(u_1,h_1)
  #              w_1(u_2,h_2), ..., w_T(u_2,h_2)
  #                          ...
  #              w_1(u_N,h_N), ..., w_T(u_N,h_N)
  
  t_len <- as.integer(t_len)
  gset  <- grid$gset
  N     <- as.integer(dim(gset)[1])
  gset  <- as.matrix(gset)
  gset  <- as.vector(gset) 
  
  storage.mode(gset) <- "double"
  
  wghts <- vector(mode = "double", length = N*T)
  
  result <- sizer_weights(t_len, gset, N)
  
  return(matrix(result, ncol=t_len, byrow=TRUE))
}


SiZer_std <- function(weights, autocov1, autocov2, t_len)

{ # compute local linear derivative estimator and its standard deviation on the
  # location-bandwidth grid.
  #
  # Arguments:
  # data      time series of length T
  # weights   kernel weights matrix produced by the function 'SiZer_weights'
  # autocov   vector of error autocovariances (gamma[0],...,gamma[T-1])
  #
  # Outputs:
  # std       vector of standard deviations (length = number of location-bandwidth
  #           points in the grid)

  autocov.mat1 <- matrix(NA, ncol=t_len, nrow=t_len)
  autocov.mat2 <- matrix(NA, ncol=t_len, nrow=t_len)

  for(ell in 1:(t_len-1)){
    autocov.mat1[ell,] <- c(autocov1[ell:1],autocov1[2:(t_len-ell+1)])
    autocov.mat2[ell,] <- c(autocov2[ell:1],autocov2[2:(t_len-ell+1)])
  }
  autocov.mat1[t_len,] <- autocov1[t_len:1]
  autocov.mat2[t_len,] <- autocov2[t_len:1]

  temp1     <- autocov.mat1 %*% t(weights)
  temp1     <- t(temp1)
  temp1     <- weights * temp1
  temp1     <- temp1 %*% rep(1,dim(temp1)[2])
  temp2     <- autocov.mat2 %*% t(weights)
  temp2     <- t(temp2)
  temp2     <- weights * temp2
  temp2     <- temp2 %*% rep(1,dim(temp2)[2])

  std.devs <- sqrt(temp1 + temp2)
  std.devs <- as.vector(std.devs)

  return(std=std.devs)
}

SiZer_quantiles <- function(alpha, t_len, grid, autocov1, autocov2)
  
{ # compute quantiles for SiZer as described in Park et al. (2009), 
  # 'Improved SiZer for time series' 
  
  gset  <- grid$gset
  u.vec <- sort(unique(gset[,1]))
  h.vec <- sort(unique(gset[,2]))
  
  Delta.tilde <- u.vec[2] - u.vec[1]
  quants      <- rep(0,length(h.vec))
  
  for(i in 1:length(h.vec)){
    gg        <- sum(gset[,2] == h.vec[i])
    
    integrand_1   <- function(s, h_, delta_, gamma_) {1000 * gamma_[floor(s * h_ / delta_) + 1] * exp(-s^2/4) * (2 - s^2)/8}
    I_gamma_num   <- 2 * integrate(integrand_1, lower = 0, upper = (t_len - 1) / (t_len * h.vec[i]),
                                   h_ = h.vec[i], delta_ = 1/t_len,
                                   gamma_ = autocov1 + autocov2,
                                   subdivisions = 500)[[1]]
    integrand_2   <- function(s, h_, delta_, gamma_) {1000 * gamma_[floor(s * h_ / delta_) + 1] * exp(-s^2/4)}
    I_gamma_denom <- 2 * integrate(integrand_2, lower = 0, upper = (t_len - 1) / (t_len * h.vec[i]),
                                   h_ = h.vec[i], delta_ = 1/t_len,
                                   gamma_ = autocov1 + autocov2,
                                   subdivisions = 500)[[1]]
    I.gamma <- I_gamma_num/I_gamma_denom
    
    # arg       <- seq(-(t_len-1), (t_len-2), by = 1)/(t_len * h.vec[i])
    # autocovs1 <- c(autocov1[t_len:2],autocov1[1:(t_len-1)])
    # autocovs2 <- c(autocov2[t_len:2],autocov2[1:(t_len-1)])
    # int1      <- sum((autocovs1 + autocovs2) * exp(-arg^2/4) * (2 - arg^2) / 8) 
    # int2      <- sum((autocovs1 + autocovs2) * exp(-arg^2/4))
    # 
    # arg       <- seq(-(t_len - 2), (t_len - 1), by = 1)/(t_len * h.vec[i])
    # autocovs1 <- c(autocov1[(t_len - 1):2], autocov1[1:t_len])
    # autocovs2 <- c(autocov2[(t_len - 1):2], autocov2[1:t_len])
    # 
    # int1      <- int1 + sum((autocovs1 + autocovs2) * exp(-arg^2/4) * (2 - arg^2) / 8) 
    # int2      <- int2 + sum((autocovs1 + autocovs2) * exp(-arg^2/4))
    #
    # I.gamma   <- int1/int2
    
    #Clustering coefficient
    theta     <- 2 * pnorm(sqrt(I.gamma) * sqrt(log(gg)) * Delta.tilde/h.vec[i]) - 1
    x         <- (1 - alpha/2)^(1/(theta * gg))
    quants[i] <- qnorm(x)
  }
  return(quants)
}


SiZer_test <- function(values1, values2, std.devs, quants, grid){ 
  
  # carry out row-wise SiZer test
  #
  # Arguments:
  # values1     vector of local linear derivative estimators of the first time series
  #             (length = number of location-bandwidth points in grid)
  # values2     vector of local linear derivative estimators of the second time series
  #             (length = number of location-bandwidth points in grid)
  # std.devs    vector of standard deviations of the local linear derivative estimators
  #             (length = number of location-bandwidth points in grid)
  # quants      vector of quantiles (length = number of bandwidth levels)
  # grid        grid of location-bandwidth points as produced by the function 'grid_construction'
  #
  # Outputs: 
  # test_sizer  matrix of SiZer test results 
  #             test_sizer[i,j] = -1: test rejects the null for the j-th location u and the 
  #                                   i-th bandwidth h and indicates a decrease in the trend
  #             test_sizer[i,j] = 0:  test does not reject the null for the j-th location u  
  #                                   and the i-th bandwidth h 
  #             test_sizer[i,j] = 1:  test rejects the null for the j-th location u and the 
  #                                   i-th bandwidth h and indicates an increase in the trend
  #             test_sizer[i,j] = 2:  no test is carried out at j-th location u and i-th 
  #                                   bandwidth h (because the point (u,h) is excluded from  
  #                                   the grid as specified by the 'deletions'-option in the
  #                                   function 'grid_construction')  
  
  gset    <- grid$gset
  h.vec   <- grid$bws   
  N       <- dim(gset)[1]
  N.h     <- length(h.vec)
  N.u     <- grid$lens
  
  quants   <- rep(quants,N.u)
  critvals <- std.devs * quants
  
  test.sizer <- rep(0,N)
  test.sizer[values1 - values2 > critvals]  <- 1
  test.sizer[values1 - values2 < -critvals] <- -1
  
  gset.full   <- grid$gset_full
  u.grid.full <- unique(gset.full[,1])
  h.grid.full <- unique(gset.full[,2])  
  pos.full    <- grid$pos_full
  
  test.full  <- rep(2,length(pos.full))  
  test.full[!is.na(pos.full)] <- test.sizer
  test.sizer <- matrix(test.full, ncol=length(u.grid.full), byrow=TRUE)
  
  return(list(ugrid=u.grid.full, hgrid=h.grid.full, test=test.sizer))
}

cluster_analysis <- function(t_len_, n_rep_, alpha_, results_matrix_){
  correct_number_of_groups   <- 0 #Starting the counter from zero
  correctly_specified_groups <- 0
  
  num_of_errors     <- c()
  
  for (i in 1:n_rep_){
    if (results_matrix_[1, i] == 3) {
      correct_number_of_groups = correct_number_of_groups + 1
    }
    if ((results_matrix_[1, i] == 2) | (results_matrix_[1, i] == 3)){
      groups123  <- results_matrix_[2:(n_ts + 1), i]
      groups132  <- recode(groups123, "2=3;3=2")
      groups213  <- recode(groups123, "1=2;2=1")
      groups231  <- recode(groups123, "1=2;2=3;3=1")
      groups312  <- recode(groups123, "1=3;2=1;3=2")
      groups321  <- recode(groups123, "1=3;3=1")
      difference <- min(sum(correct_specification != groups132),
                        sum(correct_specification != groups213),
                        sum(correct_specification != groups231),
                        sum(correct_specification != groups312),
                        sum(correct_specification != groups321),
                        sum(correct_specification != groups123))
    }
    if ((results_matrix_[1, i] == 1) | (results_matrix_[1, i] > 4)){
      difference <- 10
    }
    if (results_matrix_[1, i] == 4) {
      groups1234  <- results_matrix_[2:(n_ts + 1), i]
      groups1243  <- recode(groups1234, "4=3;3=4")
      groups1342  <- recode(groups1234, "2=3;3=4;4=2")
      groups1324  <- recode(groups1234, "2=3;3=2")
      groups1423  <- recode(groups1234, "2=4;3=2;4=3")
      groups1432  <- recode(groups1234, "2=4;4=2")
      
      groups2134  <- recode(groups1234, "1=2;2=1")
      groups2143  <- recode(groups1234, "1=2;2=1;3=4;4=3")
      groups2314  <- recode(groups1234, "1=2;2=3;3=1")
      groups2341  <- recode(groups1234, "1=2;2=3;3=4;4=1")
      groups2413  <- recode(groups1234, "1=2;2=4;3=1;4=3")
      groups2431  <- recode(groups1234, "1=2;2=4;4=1")
      
      groups3124  <- recode(groups1234, "1=3;2=1;3=2")
      groups3142  <- recode(groups1234, "1=3;2=1;3=4;4=2")
      groups3214  <- recode(groups1234, "1=3;3=1")
      groups3241  <- recode(groups1234, "1=3;3=4;4=1")
      groups3412  <- recode(groups1234, "1=3;2=4;3=1;4=2")
      groups3421  <- recode(groups1234, "1=3;2=4;3=2;4=1")
      
      groups4123  <- recode(groups1234, "1=4;2=1;3=2;4=3")
      groups4132  <- recode(groups1234, "1=4;2=1;4=2")
      groups4213  <- recode(groups1234, "1=4;3=1;4=3")
      groups4231  <- recode(groups1234, "1=4;4=1")
      groups4312  <- recode(groups1234, "1=4;2=3;3=1;4=2")
      groups4321  <- recode(groups1234, "1=4;2=3;3=2;4=1")
      
      difference <- min(sum(correct_specification != groups1234),
                        sum(correct_specification != groups1243),
                        sum(correct_specification != groups1342),
                        sum(correct_specification != groups1324),
                        sum(correct_specification != groups1423),
                        sum(correct_specification != groups1432),
                        
                        sum(correct_specification != groups2134),
                        sum(correct_specification != groups2143),
                        sum(correct_specification != groups2314),
                        sum(correct_specification != groups2341),
                        sum(correct_specification != groups2413),
                        sum(correct_specification != groups2431),
                        
                        sum(correct_specification != groups3124),
                        sum(correct_specification != groups3142),
                        sum(correct_specification != groups3214),
                        sum(correct_specification != groups3241),
                        sum(correct_specification != groups3412),
                        sum(correct_specification != groups3421),
                        
                        sum(correct_specification != groups4123),
                        sum(correct_specification != groups4132),
                        sum(correct_specification != groups4213),
                        sum(correct_specification != groups4231),
                        sum(correct_specification != groups4312),
                        sum(correct_specification != groups4321))
    }
    if (difference == 0){
      correctly_specified_groups = correctly_specified_groups + 1
    }
    num_of_errors <- c(num_of_errors, difference)
  }
  
  cat("Percentage of detecting true number of clusters",
      correct_number_of_groups/n_rep_, "with alpha = ", alpha_,
      ", T = ", t_len_, "\n")
  cat("Percentage of detecting true clustering",
      correctly_specified_groups/n_rep_, "with alpha = ", alpha_,
      ", T = ", t_len_, "\n")
  cat("Maximum number of errors is ", max(num_of_errors), "\n")
  return(list(num_of_errors = num_of_errors,
              correct_number_of_groups = correct_number_of_groups,
              correctly_specified_groups = correctly_specified_groups))
}

produce_hist_plots <- function(file_extension_, different_T_, n_rep_,
                               group_count_, error_count_){ 
  labels_ <- c("(a)", "(b)", "(c)")
  
  pdf(paste0("output/revision/hist_groups", file_extension_, ".pdf"), width = 8, height = 2.9, paper="special")
  par(mfrow = c(1,2))
  par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
  par(oma = c(1.4, 1.5, 0.5, 0.2)) #Outer margins
  
  for (j in 1:length(different_T_)){
    t_len <- different_T_[j]
    bp1 <- barplot(group_count_[[j]], ylim = c(0, 1.1 * n_rep_), xlab = "",
                   main = "", ylab = "", xaxt = 'n', space = 0)
    text(x = bp1, y = group_count_[[j]], label = group_count_[[j]], cex = 0.8, pos = 3)
    axis(1, at = bp1, labels = 1:5, tick = FALSE, line = -0.5, cex.axis = 1)
    mtext(side = 1, text= paste0(labels_[j], " T = ", t_len), line = 3.4)
    title(xlab="number of groups", mgp=c(1.5,1,0), cex.lab=1)
  }
  dev.off()
  
  pdf(paste0("output/revision/hist_errors", file_extension_, ".pdf"), width = 8, height = 2.9, paper="special")
  par(mfrow = c(1,2))
  par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
  par(oma = c(1.4, 1.5, 0.5, 0.2)) #Outer margins
  
  for (j in 1:length(different_T_)){
    t_len <- different_T_[j]
    bp2 <- barplot(error_count_[[j]], ylim = c(0, 1.1 *  n_rep_), xlab = "",
                   main = "", ylab = "", xaxt = 'n', space = 0)
    text(x = bp2, y = error_count_[[j]], label = error_count_[[j]], cex = 0.8, pos = 3)
    mtext(side = 1, text = paste0(labels_[j], " T = ", t_len), line = 3.4)
    axis(1, at = bp2, labels = 0:8, tick = FALSE, line = -0.5, cex.axis = 1)
    title(xlab = "number of errors", mgp = c(1.5,1,0), cex.lab = 1)
  }
  dev.off() 
}
