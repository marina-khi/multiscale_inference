#####################
#AUXILIARY FUNCTIONS#
#####################

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
#  if (length(data_p) != length(grid_p)){
#    cat("Dimensions of the grid and the data do not match, please check the arguments")
#    return(NULL)
#  } else {
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
#  }
}

#Simulating a bump function
bump  <- function(u){
  u.lower1 <- 0.2
  u.upper1 <- 0.4
  u.lower2 <- 0.6
  u.upper2 <- 0.8
  arg1 <- (u - 0.3)/(u.upper1 - 0.3)
  arg2 <- (u - 0.7)/(u.upper2 - 0.7)
  return(as.double(u >= u.lower1 & u <= u.upper1) * (1 - arg1^2)^2 - as.double(u >= u.lower2 & u <= u.upper2) * (1 - arg2^2)^2)
}

#Simulating a kind of a bump function for clustering
b_function <- function(u, x_0, h){
  arg <- (u - x_0)/h
  return(as.double((abs(arg) <= 1) * (1 - arg^2)^2))
}

################################
#FUNCTIONS FOR THE APPLICATIONS#
################################

produce_plots_gdp <- function(results, data_i, data_j, ticks_, labels_,
                              name_i, name_j, folder){
  filename    <- paste0(folder, name_i, "_vs_", name_j, ".pdf")
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

###############################
#FUNCTIONS FOR THE SIMULATIONS#
###############################

#Function that simulates 3 covariates as VAR(3) process with the given
#coefficients (a_x_mat_ and sigma_x_mat_),
#the error terms as AR(1) also with the given coefficient a_ and sigma_,
#the fixed effect term alpha_ as a normally distributed random vector,
#N(0, Sigma_a_mat_), the time series as
#y = alpha_ + beta_ %*% covariates + m_matrix_ + errors,
#estimates the parameters, and then computes the test statistics
repl <- function(rep_, n_ts_, t_len_, grid_, ijset_ = NULL, a_ = 0, sigma_ = 1,
                 beta_ = NULL, a_x_vec_ = c(0, 0, 0), phi_ = 0,
                 rho_ = 0, different_b_ = c(0),
                 q_ = 25, r_ = 10, type_of_m_ = "", gaussian_sim = FALSE){
  
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
                              n_ts = n_ts_, ijset = ijset_, grid = grid_)
    results <- c(as.vector(psi$stat_pairwise))
  } else {
    m_matrix        <- matrix(0, nrow = t_len_, ncol = n_ts_)    
    y_matrices      <- list()
    y_augm_matrices <- list()
    sigmahat_list   <- list()
    for (k in 1:length(different_b_)){
      y_matrices[[k]]      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      y_augm_matrices[[k]] <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      sigmahat_list[[k]]   <- rep(NA, n_ts_)  
    }
    
    error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)

    library(mvtnorm)    
    big_sigma_matrix       <- matrix(rho_, nrow = n_ts_, ncol = n_ts_)
    diag(big_sigma_matrix) <- 1
    alpha_vec              <- rmvnorm(1, mean = rep(0, n_ts_), sigma = big_sigma_matrix)

    if (!is.null(beta_)){
      phi_matrix       <- matrix(phi_, nrow = 3, ncol = 3)
      diag(phi_matrix) <- 1
      a_matrix         <- diag(a_x_vec_)      
      
      for (i in 1:n_ts_){
        error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                       innov = rnorm(t_len_, 0, sigma_),
                                       n = t_len_)
        
        nu       <- rmvnorm(t_len_ + 10, mean = c(0, 0, 0), sigma = phi_matrix)
        x_matrix <- matrix(0, 3, t_len_ + 10)
        
        for (t in 2:(t_len_ + 10)){
          x_matrix[, t] <- a_matrix %*% x_matrix[, t - 1] + nu[t, ]
        }
        x_matrix <- t(x_matrix[, -(1:10)])
        
        x_diff_1 <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
        x_diff_2 <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
        x_diff_3 <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)
        
        #Estimating beta
        x_diff_tmp <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]
        
        k <- 1
        for (b in different_b_){
          m_matrix[, 1] <- bump((1:t_len_)/t_len_) * b
          if ((type_of_m_ == "bump") & (b == 0)){
            m_matrix[, i] <- bump((1:t_len_)/t_len_) * 0.25
          }
          y_matrices[[k]][, i] <- alpha_vec[i] + m_matrix[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]

          #First differences
          y_diff     <- y_matrices[[k]][, i] - dplyr::lag(y_matrices[[k]][, i], n = 1, default = NA)
          y_diff_tmp <- as.matrix(y_diff)[-1, ]
          
          beta_hat_tmp       <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
          alpha_hat_tmp      <- mean(y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp))
          
          y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp) - alpha_hat_tmp
          
          AR.struc           <- estimate_lrv(data = y_augm_matrices[[k]][, i], q = q_,
                                             r_bar = r_, p = 1)
          sigma_hat_i        <- sqrt(AR.struc$lrv)
          sigmahat_list[[k]][i] <- sigma_hat_i 
          k <- k + 1
        }
      }
    } else {
      for (i in 1:n_ts_){
        error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                       innov = rnorm(t_len_, 0, sigma_),
                                       n = t_len_)
        
        k <- 1
        for (b in different_b_){
          m_matrix[, 1]        <- bump((1:t_len_)/t_len_) * b
          y_matrices[[k]][, i] <- alpha_vec[i] + m_matrix[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
        
          alpha_hat_tmp             <- mean(y_matrices[[k]][, i])
          y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - alpha_hat_tmp
          
          AR.struc              <- estimate_lrv(data = y_augm_matrices[[k]][, i], q = q_,
                                                r_bar = r_, p = 1)
          sigma_hat_i           <- sqrt(AR.struc$lrv)
          sigmahat_list[[k]][i] <- sigma_hat_i 
          k <- k + 1
        }
      }
    }
    results <- c()
    for (k in 1:length(different_b_)){
      psi     <- compute_statistics(data = y_augm_matrices[[k]],
                                    sigma_vec = sigmahat_list[[k]],
                                    n_ts = n_ts_, ijset = ijset_, grid = grid_)    
      results <- c(results, as.vector(psi$stat_pairwise))
    }
  }
  return(results)
}

#Function that simulates 3 covariates as VAR(3) process with the given
#coefficients (a_x_mat_ and sigma_x_mat_),
#the error terms as AR(1) also with the given coefficient a_ and sigma_,
#the fixed effect term alpha_ as a normally distributed random vector,
#N(0, Sigma_a_mat_), the time series as
#y = alpha_ + beta_ %*% covariates + m_matrix_ + errors,
#estimates the parameters, and then computes the test statistics
repl_SV <- function(rep_, n_ts_, t_len_, grid_, ijset_ = NULL, a_ = 0,
                    beta_ = NULL, a_x_vec_ = c(0, 0, 0), phi_ = 0,
                    rho_ = 0, different_b_ = c(0),
                    q_ = 25, r_ = 10, type_of_m_ = "", gaussian_sim = FALSE){
  
  library(MSinference)
  library(dplyr)
  
  sigma_tv_vec <- -0.15 * ((1:t_len_)/t_len_ - 0.5)^2 + 0.075
  
  if (gaussian_sim){
    z_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    z_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    sigma_vector  <- rep(1, n_ts_)
    
    for (i in 1:n_ts_){
      z_matrix[, i]      <- rnorm(t_len_, 0, 1)
      z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    }
    
    psi <- compute_statistics(data = z_augm_matrix,
                              sigma_vec = sigma_vector,
                              n_ts = n_ts_, ijset = ijset_, grid = grid_)
    results <- c(as.vector(psi$stat_pairwise))
    
    # z_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    # z_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    # 
    # for (i in 1:n_ts_){
    #   z_matrix[, i]      <- rnorm(t_len_, 0, sqrt(sigma_tv_vec))
    #   z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    # }
    # 
    # psi <- compute_statistics(data = z_augm_matrix,
    #                           sigma_vec = sqrt(sigma_tv_vec),
    #                           n_ts = n_ts_, ijset = ijset_, grid = grid_)
    # results <- c(as.vector(psi$stat_pairwise))
  } else {
    m_matrix        <- matrix(0, nrow = t_len_, ncol = n_ts_)    
    y_matrices      <- list()
    y_augm_matrices <- list()
    sigmahat_list   <- list()
    for (k in 1:length(different_b_)){
      y_matrices[[k]]      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      y_augm_matrices[[k]] <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      sigmahat_list[[k]]   <- rep(NA, n_ts_)  
    }
    
    error_matrix     <- matrix(0, nrow = t_len_, ncol = n_ts_)
    error_matrix_tmp <- matrix(0, nrow = t_len_ + 10, ncol = n_ts_)
        
    library(mvtnorm)    
    big_sigma_matrix       <- matrix(rho_, nrow = n_ts_, ncol = n_ts_)
    diag(big_sigma_matrix) <- 1
    alpha_vec              <- rmvnorm(1, mean = rep(0, n_ts_), sigma = big_sigma_matrix)
    
    if (!is.null(beta_)){
      phi_matrix       <- matrix(phi_, nrow = 3, ncol = 3)
      diag(phi_matrix) <- 1
      a_matrix         <- diag(a_x_vec_)      
      
      for (i in 1:n_ts_){
        eta      <- rnorm(t_len_ + 10, mean = 0, sd = 1)
        eta_SV   <- c(rep(1, 10), sqrt(sigma_tv_vec)) * eta
        nu       <- rmvnorm(t_len_ + 10, mean = c(0, 0, 0), sigma = phi_matrix)
        x_matrix <- matrix(0, 3, t_len_ + 10)
        
        for (t in 2:(t_len_ + 10)){
          x_matrix[, t] <- a_matrix %*% x_matrix[, t - 1] + nu[t, ]
          error_matrix_tmp[t, i] <- a_ * error_matrix_tmp[t - 1, i] + eta_SV[t]
        }
        
        error_matrix[, i] <- error_matrix_tmp[-(1:10), i]
        x_matrix <- t(x_matrix[, -(1:10)])
        
        x_diff_1 <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
        x_diff_2 <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
        x_diff_3 <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)
        
        #Estimating beta
        x_diff_tmp <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]
        
        k <- 1
        for (b in different_b_){
          m_matrix[, 1] <- bump((1:t_len_)/t_len_) * b
          if ((type_of_m_ == "bump") & (b == 0)){
            m_matrix[, i] <- bump((1:t_len_)/t_len_) * 0.25
          }
          y_matrices[[k]][, i] <- alpha_vec[i] + m_matrix[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
          
          #First differences
          y_diff     <- y_matrices[[k]][, i] - dplyr::lag(y_matrices[[k]][, i], n = 1, default = NA)
          y_diff_tmp <- as.matrix(y_diff)[-1, ]
          
          beta_hat_tmp       <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
          alpha_hat_tmp      <- mean(y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp))
          
          y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp) - alpha_hat_tmp
          
          AR.struc           <- estimate_lrv(data = y_augm_matrices[[k]][, i], q = q_,
                                             r_bar = r_, p = 1)
          sigma_hat_i        <- sqrt(AR.struc$lrv)
          sigmahat_list[[k]][i] <- sigma_hat_i 
          k <- k + 1
        }
      }
    } else {
      for (i in 1:n_ts_){
        eta      <- rnorm(t_len_ + 10, mean = 0, sd = 1)
        eta_SV   <- c(rep(1, 10), sqrt(sigma_tv_vec)) * eta

        for (t in 2:(t_len_ + 10)){
          error_matrix_tmp[t, i] <- a_ * error_matrix_tmp[t - 1, i] + eta_SV[t]
        }
        error_matrix[, i] <- error_matrix_tmp[-(1:10), i]
        
        k <- 1
        for (b in different_b_){
          m_matrix[, 1]        <- bump((1:t_len_)/t_len_) * b
          y_matrices[[k]][, i] <- alpha_vec[i] + m_matrix[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
          
          alpha_hat_tmp             <- mean(y_matrices[[k]][, i])
          y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - alpha_hat_tmp
          
          AR.struc              <- estimate_lrv(data = y_augm_matrices[[k]][, i], q = q_,
                                                r_bar = r_, p = 1)
          sigma_hat_i           <- sqrt(AR.struc$lrv)
          sigmahat_list[[k]][i] <- sigma_hat_i 
          k <- k + 1
        }
      }
    }
    results <- c()
    for (k in 1:length(different_b_)){
      psi     <- compute_statistics(data = y_augm_matrices[[k]],
                                    sigma_vec = sigmahat_list[[k]],
                                    n_ts = n_ts_, ijset = ijset_, grid = grid_)    
      results <- c(results, as.vector(psi$stat_pairwise))
    }
  }
  return(results)
}

repl_spurious <- function(rep_, n_ts_, t_len_, grid1_, grid2_,
                          ijset1_ = NULL, ijset2_ = NULL, a_ = 0, sigma_ = 1,
                          beta_ = NULL, a_x_vec_ = c(0, 0, 0), phi_ = 0,
                          rho_ = 0, different_b_ = c(0),
                          q_ = 25, r_ = 10, type_of_m_ = ""){
  
  library(MSinference)
  library(dplyr)
  
  m_matrix        <- matrix(0, nrow = t_len_, ncol = n_ts_)    
  y_matrices      <- list()
  y_augm_matrices <- list()
  sigmahat_list   <- list()
  for (k in 1:length(different_b_)){
    y_matrices[[k]]      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    y_augm_matrices[[k]] <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    sigmahat_list[[k]]   <- rep(NA, n_ts_)  
  }
  
  error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)
  
  library(mvtnorm)    
  big_sigma_matrix       <- matrix(rho_, nrow = n_ts_, ncol = n_ts_)
  diag(big_sigma_matrix) <- 1
  alpha_vec              <- rmvnorm(1, mean = rep(0, n_ts_), sigma = big_sigma_matrix)
  
  if (!is.null(beta_)){
    phi_matrix       <- matrix(phi_, nrow = 3, ncol = 3)
    diag(phi_matrix) <- 1
    a_matrix         <- diag(a_x_vec_)      
    
    for (i in 1:n_ts_){
      error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                     innov = rnorm(t_len_, 0, sigma_),
                                     n = t_len_)
      
      nu       <- rmvnorm(t_len_ + 10, mean = c(0, 0, 0), sigma = phi_matrix)
      x_matrix <- matrix(0, 3, t_len_ + 10)
      
      for (t in 2:(t_len_ + 10)){
        x_matrix[, t] <- a_matrix %*% x_matrix[, t - 1] + nu[t, ]
      }
      x_matrix <- t(x_matrix[, -(1:10)])
      
      x_diff_1 <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
      x_diff_2 <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
      x_diff_3 <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)
      
      #Estimating beta
      x_diff_tmp <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]
      
      k <- 1
      for (b in different_b_){
        m_matrix[, 1] <- bump((1:t_len_)/t_len_) * b
        if ((type_of_m_ == "bump") & (b == 0)){
          m_matrix[, i] <- bump((1:t_len_)/t_len_) * 0.25
        }
        y_matrices[[k]][, i] <- alpha_vec[i] + m_matrix[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
        
        #First differences
        y_diff     <- y_matrices[[k]][, i] - dplyr::lag(y_matrices[[k]][, i], n = 1, default = NA)
        y_diff_tmp <- as.matrix(y_diff)[-1, ]
        
        beta_hat_tmp       <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
        alpha_hat_tmp      <- mean(y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp))
        
        y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp) - alpha_hat_tmp
        
        AR.struc           <- estimate_lrv(data = y_augm_matrices[[k]][, i], q = q_,
                                           r_bar = r_, p = 1)
        sigma_hat_i        <- sqrt(AR.struc$lrv)
        sigmahat_list[[k]][i] <- sigma_hat_i 
        k <- k + 1
      }
    }
  } else {
    for (i in 1:n_ts_){
      error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                     innov = rnorm(t_len_, 0, sigma_),
                                     n = t_len_)
      
      k <- 1
      for (b in different_b_){
        m_matrix[, 1]        <- bump((1:t_len_)/t_len_) * b
        y_matrices[[k]][, i] <- alpha_vec[i] + m_matrix[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
        
        alpha_hat_tmp             <- mean(y_matrices[[k]][, i])
        y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - alpha_hat_tmp
        
        AR.struc              <- estimate_lrv(data = y_augm_matrices[[k]][, i], q = q_,
                                              r_bar = r_, p = 1)
        sigma_hat_i           <- sqrt(AR.struc$lrv)
        sigmahat_list[[k]][i] <- sigma_hat_i 
        k <- k + 1
      }
    }
  }
  results <- c()
  for (k in 1:length(different_b_)){
    psi1    <- compute_statistics(data = y_augm_matrices[[k]],
                                  sigma_vec = sigmahat_list[[k]],
                                  n_ts = n_ts_, ijset = ijset1_, grid = grid1_)
    psi2    <- compute_statistics(data = y_augm_matrices[[k]],
                                  sigma_vec = sigmahat_list[[k]],
                                  n_ts = n_ts_, ijset = ijset2_, grid = grid2_) 
    results <- c(results, as.vector(psi1$stat_pairwise), as.vector(psi2$stat_pairwise))
  }
  return(results)
}

repl_precision <- function(rep_, n_ts_, t_len_, 
                           a_ = 0, sigma_ = 1,
                           beta_ = NULL, a_x_vec_ = c(0, 0, 0), phi_ = 0,
                           rho_ = 0, different_b_ = c(0),
                           q_ = 25, r_ = 10, type_of_m_ = ""){
  
  library(MSinference)
  library(dplyr)
  
  m_matrix        <- matrix(0, nrow = t_len_, ncol = n_ts_)    
  y_matrices      <- list()
  y_augm_matrices <- list()
  sigmahat_list   <- list()
  betahat_list    <- list()
  for (k in 1:length(different_b_)){
    y_matrices[[k]]      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    y_augm_matrices[[k]] <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    betahat_list[[k]]    <- matrix(NA, nrow = length(beta_), ncol = n_ts_)
    sigmahat_list[[k]]   <- rep(NA, n_ts_)  
  }
  
  error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)
  
  library(mvtnorm)    
  big_sigma_matrix       <- matrix(rho_, nrow = n_ts_, ncol = n_ts_)
  diag(big_sigma_matrix) <- 1
  alpha_vec              <- rmvnorm(1, mean = rep(0, n_ts_), sigma = big_sigma_matrix)
  
  phi_matrix       <- matrix(phi_, nrow = 3, ncol = 3)
  diag(phi_matrix) <- 1
  a_matrix         <- diag(a_x_vec_)

  for (i in 1:n_ts_){
    error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                   innov = rnorm(t_len_, 0, sigma_),
                                   n = t_len_)
    
    nu       <- rmvnorm(t_len_ + 10, mean = c(0, 0, 0), sigma = phi_matrix)
    x_matrix <- matrix(0, 3, t_len_ + 10)
    
    for (t in 2:(t_len_ + 10)){
      x_matrix[, t] <- a_matrix %*% x_matrix[, t - 1] + nu[t, ]
    }
    x_matrix <- t(x_matrix[, -(1:10)])
    
    x_diff_1 <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
    x_diff_2 <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
    x_diff_3 <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)
    
    #Estimating beta
    x_diff_tmp <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]
    
    k <- 1
    for (b in different_b_){
      m_matrix[, 1] <- bump((1:t_len_)/t_len_) * b
      if ((type_of_m_ == "bump") & (b == 0)){
        m_matrix[, i] <- bump((1:t_len_)/t_len_) * 0.25
      }
      y_matrices[[k]][, i] <- alpha_vec[i] + m_matrix[, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
      
      #First differences
      y_diff     <- y_matrices[[k]][, i] - dplyr::lag(y_matrices[[k]][, i], n = 1, default = NA)
      y_diff_tmp <- as.matrix(y_diff)[-1, ]
      
      beta_hat_tmp           <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp
      alpha_hat_tmp          <- mean(y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp))
      betahat_list[[k]][, i] <- beta_hat_tmp
      
      y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp) - alpha_hat_tmp
      
      AR.struc           <- estimate_lrv(data = y_augm_matrices[[k]][, i], q = q_,
                                         r_bar = r_, p = 1)
      sigma_hat_i        <- sqrt(AR.struc$lrv)
      sigmahat_list[[k]][i] <- sigma_hat_i 
      k <- k + 1
    }
  }
  results <- c()
  for (k in 1:length(different_b_)){
    results <- c(results, as.vector(rowMeans(betahat_list[[k]])), as.vector(mean(sigmahat_list[[k]])))
  }
  return(results)
}


#Function that simulates 3 covariates as VAR(3) process with the given
#coefficients (a_x_mat_ and sigma_x_mat_),
#the error terms as AR(1) also with the given coefficient a_ and sigma_,
#the fixed effect term alpha_ as a normally distributed random vector,
#N(0, Sigma_a_mat_), the time series as
#y = alpha_ + beta_ %*% covariates + m_matrix_ + errors,
#estimates the parameters, and then computes the test statistics
repl_UCB <- function(rep_, n_ts_, t_len_, bw_ = 0.1, a_ = 0, sigma_ = 1,
                     beta_ = NULL, a_x_vec_ = c(0, 0, 0), phi_ = 0,
                     different_b_ = c(0), quant_UCB_ = 0,
                     gaussian_sim = FALSE){
  
  grid_points <- seq(from = 1 / t_len_, to = 1, by = 1 / t_len_)  
  
  if (gaussian_sim){
    z_vec        <- rnorm(t_len_, 0, 1)
    gaussian_UCB <- mapply(UCB_estimation, grid_points,
                           MoreArgs = list(data_p = z_vec,
                                           grid_p = grid_points,
                                           bw = bw_))
    results <- c(max(abs(gaussian_UCB)))
  } else {
    library(mvtnorm)
    library(dplyr)
    
    results <- c()
    
    k_n <- floor(t_len_^(1/3))
    m   <- floor(t_len_ / k_n)
    lower_boundary <- ceiling(0.05 * t_len_)
    upper_boundary <- floor(0.95 * t_len_)

    m_matrices          <- list()
    y_matrices          <- list()
    y_augm_matrices     <- list()
    estimated_trend_UCB <- list()
    upper_UCB           <- list()
    lower_UCB           <- list()
    
    for (k in 1:length(different_b_)){
      y_matrices[[k]]          <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      m_matrices[[k]]          <- matrix(0, nrow = t_len_, ncol = n_ts_)
      m_matrices[[k]][, 1]     <- bump((1:t_len_)/t_len_) * different_b_[k]
      y_augm_matrices[[k]]     <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      estimated_trend_UCB[[k]] <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      upper_UCB[[k]]           <- matrix(NA, nrow = t_len_, ncol = n_ts_)
      lower_UCB[[k]]           <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    }
    
    error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)

    phi_matrix       <- matrix(phi_, nrow = 3, ncol = 3)
    diag(phi_matrix) <- 1
    a_matrix         <- diag(a_x_vec_)      
      
    for (i in 1:n_ts_){
      error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                     innov = rnorm(t_len_, 0, sigma_),
                                     n = t_len_)
      
      nu       <- rmvnorm(t_len_ + 10, mean = c(0, 0, 0), sigma = phi_matrix)
      x_matrix <- matrix(0, 3, t_len_ + 10)
      
      for (t in 2:(t_len_ + 10)){
        x_matrix[, t] <- a_matrix %*% x_matrix[, t - 1] + nu[t, ]
      }
      x_matrix <- t(x_matrix[, -(1:10)])
        
      x_diff_1 <- x_matrix[, 1] - dplyr::lag(x_matrix[, 1], n = 1, default = NA)
      x_diff_2 <- x_matrix[, 2] - dplyr::lag(x_matrix[, 2], n = 1, default = NA)
      x_diff_3 <- x_matrix[, 3] - dplyr::lag(x_matrix[, 3], n = 1, default = NA)
        
      #Estimating beta
      x_diff_tmp <- as.matrix(cbind(x_diff_1, x_diff_2, x_diff_3))[-1, ]
        
      k <- 1
      for (b in different_b_){
        y_matrices[[k]][, i] <- m_matrices[[k]][, i] + beta_ %*% t(x_matrix) + error_matrix[, i]
          
        #First differences
        y_diff     <- y_matrices[[k]][, i] - dplyr::lag(y_matrices[[k]][, i], n = 1, default = NA)
        y_diff_tmp <- as.matrix(y_diff)[-1, ]
          
        beta_hat_tmp <- solve(t(x_diff_tmp) %*% x_diff_tmp) %*% t(x_diff_tmp) %*% y_diff_tmp

        y_augm_matrices[[k]][, i] <- y_matrices[[k]][, i] - x_matrix %*% as.vector(beta_hat_tmp)
          
        sigma_hat_UCB_i <- sqrt(sigma_estimation_UCB(y_ = y_matrices[[k]][, i],
                                                     x_matrix_ = x_matrix,
                                                     beta_est_ = beta_hat_tmp,
                                                     m_ = m, k_n_ = k_n)) 

        estimated_trend_UCB[[k]][, i] <- mapply(UCB_estimation, grid_points,
                                                MoreArgs = list(data_p = y_augm_matrices[[k]][, i],
                                                                grid_p = grid_points,
                                                                bw = bw_))
        upper_UCB[[k]][, i] <- estimated_trend_UCB[[k]][, i] + sigma_hat_UCB_i * quant_UCB_
        lower_UCB[[k]][, i] <- estimated_trend_UCB[[k]][, i] - sigma_hat_UCB_i * quant_UCB_
        
        k <- k + 1
      }
    }

    for (k in 1:length(different_b_)){
      pairwise_intersection_UCB <- matrix(NA, ncol = n_ts_, nrow = n_ts_)
      
      #Ignoring the boundary issues
      tmp_u <- upper_UCB[[k]][lower_boundary:upper_boundary, ]
      tmp_l <- lower_UCB[[k]][lower_boundary:upper_boundary, ]

      for (i in 1:n_ts_){
        for (j in 1:n_ts_){
          pairwise_intersection_UCB[i, j] <- sum((tmp_u[, i] < tmp_l[, j]) | (tmp_u[, i] < tmp_l[, j]))
        }
      }
      results <- c(results, as.vector(pairwise_intersection_UCB))
    }
  }
  return(results)
}

#Function used for simulated the three groups of time series and calculate 
#the corresponding test statistics, needed for parallel computations
repl_clustering <- function(rep, t_len_, n_ts_, grid_,
                            m_matrix_ = NULL, 
                            a_ = 0, sigma_ = 1,
                            beta_ = NULL, a_x_vec_ = c(0, 0, 0),
                            phi_ = 0,  rho_ = 0,
                            q_ = 25, r_ = 10, 
                            gaussian_sim = FALSE, comparison = FALSE){
  library(MSinference)

  if (gaussian_sim){
    z_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    z_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)

    for (i in 1:n_ts_){
      z_matrix[, i]      <- rnorm(t_len_, 0, sigma_)
      z_augm_matrix[, i] <- z_matrix[, i] - mean(z_matrix[, i])
    }
    
    psi <- compute_statistics(data = z_augm_matrix,
                              sigma_vec = rep(sigma_, n_ts_),
                              n_ts = n_ts_, grid = grid_)
    results <- c(as.vector(psi$stat_pairwise))
  } else {
    library(mvtnorm)
    y_matrix      <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    y_augm_matrix <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    error_matrix  <- matrix(NA, nrow = t_len_, ncol = n_ts_)
    
    sigmahat_vector <- c()
    
    if (!is.null(beta_)){    
      big_sigma_matrix       <- matrix(rho_, nrow = n_ts_, ncol = n_ts_)
      diag(big_sigma_matrix) <- 1
      alpha_vec              <- rmvnorm(1, mean = rep(0, n_ts_), sigma = big_sigma_matrix)
      
      phi_matrix       <- matrix(phi_, nrow = 3, ncol = 3)
      diag(phi_matrix) <- 1
      a_matrix         <- diag(a_x_vec_)
      
      for (i in 1:n_ts_){
        error_matrix[, i] <- arima.sim(model = list(ar = a_),
                                       innov = rnorm(t_len_, 0, sigma_),
                                       n = t_len_)
        
        nu       <- rmvnorm(t_len_ + 10, mean = c(0, 0, 0), sigma = phi_matrix)
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
        error_matrix[, i]  <- arima.sim(model = list(ar = a_),
                                        innov = rnorm(t_len_, 0, sigma_),
                                        n = t_len_)
        y_matrix[, i]      <- m_matrix_[, i] + error_matrix[, i]
        
        #Estimating the fixed effects
        alpha_hat_tmp      <- mean(y_matrix[, i])
        y_augm_matrix[, i] <- y_matrix[, i] - alpha_hat_tmp
        AR.struc           <- estimate_lrv(data = y_augm_matrix[, i], q = q_,
                                           r_bar = r_, p = 1)
        sigma_hat_i        <- sqrt(AR.struc$lrv)
        sigmahat_vector    <- c(sigmahat_vector, sigma_hat_i) 
      }
    }    
    psi <- compute_statistics(data = y_augm_matrix, sigma_vec = sigmahat_vector,
                              n_ts = n_ts_, grid = grid_)
    
    if (comparison){
      #psi_estimated <- compute_statistics(data = y_augm_matrix,
      #                                    sigma_vec = sigmahat_vector,
      #                                    n_ts = n_ts_, grid = grid_)
      #sigma_vector  <- rep(sqrt(sigma_^2/((1 - a_)^2)), n_ts_)      
      #psi_true      <- compute_statistics(data = y_augm_matrix,
      #                                    sigma_vec = sigma_vector,
      #                                    n_ts = n_ts_, grid = grid_)
      
      grid_points    <- seq(1/t_len_, 1, by = 1/t_len_)
      u_grid         <- unique(grid_$gset[, 1])
      smoothed_data  <- matrix(NA, nrow = length(u_grid), ncol = n_ts_)
      smoothed_data2 <- matrix(NA, nrow = length(u_grid), ncol = n_ts_)
      smoothed_data3 <- matrix(NA, nrow = length(u_grid), ncol = n_ts_)
      smoothed_data4 <- matrix(NA, nrow = length(u_grid), ncol = n_ts_)
      smoothed_data5 <- matrix(NA, nrow = length(u_grid), ncol = n_ts_)
      smoothed_data6 <- matrix(NA, nrow = length(u_grid), ncol = n_ts_)
      
      for (i in 1:n_ts_){
        h_min <- min(grid_$bws)
        h_max <- max(grid_$bws)
        smoothed_data[, i]   <- mapply(local_linear_smoothing, u_grid, MoreArgs = list(y_augm_matrix[, i], grid_points, h_min))
        smoothed_data2[, i]  <- mapply(local_linear_smoothing, u_grid, MoreArgs = list(y_augm_matrix[, i], grid_points, h_min + 0.2 * (h_max - h_min)))
        smoothed_data3[, i]  <- mapply(local_linear_smoothing, u_grid, MoreArgs = list(y_augm_matrix[, i], grid_points, h_min + 0.4 * (h_max - h_min)))
        smoothed_data4[, i]  <- mapply(local_linear_smoothing, u_grid, MoreArgs = list(y_augm_matrix[, i], grid_points, h_min + 0.6 * (h_max - h_min)))
        smoothed_data5[, i]  <- mapply(local_linear_smoothing, u_grid, MoreArgs = list(y_augm_matrix[, i], grid_points, h_min + 0.8 * (h_max - h_min)))
        smoothed_data6[, i]  <- mapply(local_linear_smoothing, u_grid, MoreArgs = list(y_augm_matrix[, i], grid_points, h_max))
      }
      
      #Calculating the benchmark model
      benchmark_results  <- matrix(0, nrow = n_ts_, ncol = n_ts_)
      benchmark_results2 <- matrix(0, nrow = n_ts_, ncol = n_ts_)
      benchmark_results3 <- matrix(0, nrow = n_ts_, ncol = n_ts_)
      benchmark_results4 <- matrix(0, nrow = n_ts_, ncol = n_ts_)
      benchmark_results5 <- matrix(0, nrow = n_ts_, ncol = n_ts_)
      benchmark_results6 <- matrix(0, nrow = n_ts_, ncol = n_ts_)
      for (i in 1:(n_ts_ - 1)){
        for (j in (i + 1):n_ts_){
          #benchmark_results[i, j]  <- sum((smoothed_data[, i] - smoothed_data[, j])^2) * (1/t_len_)
          #benchmark_results2[i, j] <- sum((smoothed_data2[, i] - smoothed_data2[, j])^2) * (1/t_len_)
          benchmark_results[i, j]  <- max(abs(smoothed_data[, i] - smoothed_data[, j]))
          benchmark_results2[i, j] <- max(abs(smoothed_data2[, i] - smoothed_data2[, j]))
          benchmark_results3[i, j] <- max(abs(smoothed_data3[, i] - smoothed_data3[, j]))
          benchmark_results4[i, j] <- max(abs(smoothed_data4[, i] - smoothed_data4[, j]))
          benchmark_results5[i, j] <- max(abs(smoothed_data5[, i] - smoothed_data5[, j]))
          benchmark_results6[i, j] <- max(abs(smoothed_data6[, i] - smoothed_data6[, j]))
        }
      }
      results <- c(as.vector(psi$stat_pairwise),
                   #as.vector(psi_estimated$stat_pairwise),
                   #as.vector(psi_true$stat_pairwise),
                   as.vector(benchmark_results), as.vector(benchmark_results2),
                   as.vector(benchmark_results3), as.vector(benchmark_results4),
                   as.vector(benchmark_results5), as.vector(benchmark_results6))
    } else {
      results <- c(as.vector(psi$stat_pairwise))
    }
  }
  return(results)
}

cluster_analysis <- function(t_len_, n_rep_, alpha_, results_matrix_,
                             correct_specification_){
  library(car)
  correct_number_of_groups   <- 0 #Starting the counter from zero
  correctly_specified_groups <- 0
  num_of_errors              <- c()
  
  for (i in 1:n_rep_){
    if (results_matrix_[1, i] == 3) {
      correct_number_of_groups = correct_number_of_groups + 1
    }
    if ((results_matrix_[1, i] == 2) | (results_matrix_[1, i] == 3)){
      groups123  <- results_matrix_[2:(n_ts + 1), i]
      groups132  <- car::recode(groups123, "2=3;3=2")
      groups213  <- car::recode(groups123, "1=2;2=1")
      groups231  <- car::recode(groups123, "1=2;2=3;3=1")
      groups312  <- car::recode(groups123, "1=3;2=1;3=2")
      groups321  <- car::recode(groups123, "1=3;3=1")
      difference <- min(sum(correct_specification_ != groups132),
                        sum(correct_specification_ != groups213),
                        sum(correct_specification_ != groups231),
                        sum(correct_specification_ != groups312),
                        sum(correct_specification_ != groups321),
                        sum(correct_specification_ != groups123))
    }
    if ((results_matrix_[1, i] == 1) | (results_matrix_[1, i] > 4)){
      difference <- 10
    }
    if (results_matrix_[1, i] == 4) {
      groups1234  <- results_matrix_[2:(n_ts + 1), i]
      groups1243  <- car::recode(groups1234, "4=3;3=4")
      groups1342  <- car::recode(groups1234, "2=3;3=4;4=2")
      groups1324  <- car::recode(groups1234, "2=3;3=2")
      groups1423  <- car::recode(groups1234, "2=4;3=2;4=3")
      groups1432  <- car::recode(groups1234, "2=4;4=2")
      
      groups2134  <- car::recode(groups1234, "1=2;2=1")
      groups2143  <- car::recode(groups1234, "1=2;2=1;3=4;4=3")
      groups2314  <- car::recode(groups1234, "1=2;2=3;3=1")
      groups2341  <- car::recode(groups1234, "1=2;2=3;3=4;4=1")
      groups2413  <- car::recode(groups1234, "1=2;2=4;3=1;4=3")
      groups2431  <- car::recode(groups1234, "1=2;2=4;4=1")
      
      groups3124  <- car::recode(groups1234, "1=3;2=1;3=2")
      groups3142  <- car::recode(groups1234, "1=3;2=1;3=4;4=2")
      groups3214  <- car::recode(groups1234, "1=3;3=1")
      groups3241  <- car::recode(groups1234, "1=3;3=4;4=1")
      groups3412  <- car::recode(groups1234, "1=3;2=4;3=1;4=2")
      groups3421  <- car::recode(groups1234, "1=3;2=4;3=2;4=1")
      
      groups4123  <- car::recode(groups1234, "1=4;2=1;3=2;4=3")
      groups4132  <- car::recode(groups1234, "1=4;2=1;4=2")
      groups4213  <- car::recode(groups1234, "1=4;3=1;4=3")
      groups4231  <- car::recode(groups1234, "1=4;4=1")
      groups4312  <- car::recode(groups1234, "1=4;2=3;3=1;4=2")
      groups4321  <- car::recode(groups1234, "1=4;2=3;3=2;4=1")
      
      difference <- min(sum(correct_specification_ != groups1234),
                        sum(correct_specification_ != groups1243),
                        sum(correct_specification_ != groups1342),
                        sum(correct_specification_ != groups1324),
                        sum(correct_specification_ != groups1423),
                        sum(correct_specification_ != groups1432),
                        
                        sum(correct_specification_ != groups2134),
                        sum(correct_specification_ != groups2143),
                        sum(correct_specification_ != groups2314),
                        sum(correct_specification_ != groups2341),
                        sum(correct_specification_ != groups2413),
                        sum(correct_specification_ != groups2431),
                        
                        sum(correct_specification_ != groups3124),
                        sum(correct_specification_ != groups3142),
                        sum(correct_specification_ != groups3214),
                        sum(correct_specification_ != groups3241),
                        sum(correct_specification_ != groups3412),
                        sum(correct_specification_ != groups3421),
                        
                        sum(correct_specification_ != groups4123),
                        sum(correct_specification_ != groups4132),
                        sum(correct_specification_ != groups4213),
                        sum(correct_specification_ != groups4231),
                        sum(correct_specification_ != groups4312),
                        sum(correct_specification_ != groups4321))
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


#Create a matrix (for size and power table for example) and write them in the tex file
output_matrix <- function(matrix_, filename_, numcols_){
  addtorow     <- list()
  addtorow$pos <- list(0, 0)
  if (numcols_ == 4) {
    addtorow$command <- c("& \\multicolumn{3}{c}{nominal size $\\alpha$} \\\\\n",
                          "$T$ & 0.01 & 0.05 & 0.1 \\\\\n") 
    print.xtable(xtable(matrix_, digits = c(3), align = "cccc"), type = "latex",
                 file = filename_, add.to.row = addtorow, include.colnames = FALSE,
                 sanitize.text.function=function(x){x})
  } else if (numcols_ == 7) {
    addtorow$command <- c("& \\multicolumn{6}{c}{nominal size $\\alpha$} \\\\\n",
                          "$T$ & 0.01 & 0.05 & 0.1 & 0.01 & 0.05 & 0.1 \\\\\n") 
    print.xtable(xtable(matrix_, digits = c(3), align = "ccccccc"), type = "latex",
                 file = filename_, add.to.row = addtorow, include.colnames = FALSE,
                 sanitize.text.function=function(x){x})
  } else if (numcols_ == 8){
    addtorow$pos     <- list(0)
    addtorow$command <- c("& $\\mathcal{T}_{\\text{MS}}$ & $\\mathcal{T}_{\\text{bmk, }1}$ & $\\mathcal{T}_{\\text{bmk, }2}$ & $\\mathcal{T}_{\\text{bmk, }3}$ & $\\mathcal{T}_{\\text{bmk, }4}$ & $\\mathcal{T}_{\\text{bmk, }5}$ & $\\mathcal{T}_{\\text{bmk, }6}$\\\\\n")
    
    print.xtable(xtable(matrix_, digits = c(3), align = "cccccccc"), type = "latex",
                 file = filename_, add.to.row = addtorow, include.colnames = FALSE,
                 sanitize.text.function=function(x){x})
  } else {cat("Number of columns not supported")}
}

produce_hist_plots <- function(file_extension_, different_T_, n_rep_,
                               group_count_, error_count_){ 
  labels_ <- c("(a)", "(b)", "(c)")
  
  pdf(paste0("output/revision/hist_groups", file_extension_, ".pdf"), width = 8, height = 2.9, paper="special")
  par(mfrow = c(1,3))
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
  par(mfrow = c(1,3))
  par(mar = c(3, 2, 0.5, 1)) #Margins for each plot
  par(oma = c(1.4, 1.5, 0.5, 0.2)) #Outer margins
  
  for (j in 1:length(different_T_)){
    t_len <- different_T_[j]
    bp2 <- barplot(error_count_[[j]], ylim = c(0, 1.1 *  n_rep_), xlab = "",
                   main = "", ylab = "", xaxt = 'n', space = 0)
    text(x = bp2, y = error_count_[[j]], label = error_count_[[j]], cex = 0.8, pos = 3)
    mtext(side = 1, text = paste0(labels_[j], " T = ", t_len), line = 3.4)
    axis(1, at = bp2, labels = 0:10, tick = FALSE, line = -0.5, cex.axis = 0.8)
    title(xlab = "number of errors", mgp = c(1.5,1,0), cex.lab = 1)
  }
  dev.off() 
}

