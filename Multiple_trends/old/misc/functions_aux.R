#' Calculates estimators of the standard deviations for a number of iid samples. 
#' It uses the standard first differences method.
#' @param data      Matrix of the samples. Each sample is located in one of the columns.
#' @param method    Method of estimation, either 'first', then first differences are used, or 'second',
#'                  then the two neighbours method is used, as described in Brown and Levine (2007)
#' @return lrv      Estimator of the long run variance of the error terms.
#' @return ahat     Vector of length p of estimated AR coefficients.
#' @return vareta   Estimator of the variance of the innovation term

estimate_iid_sd <- function(data, method = 'first'){
  N    = ncol(data)
  Tlen = nrow(data)
  sigmahat_vec <- c()
  if (method == 'first') {
    for (i in 1:N){
      variance_i   <- sum((data[2:Tlen, i] - data[1:(Tlen - 1), i])^2)/(2 * Tlen - 2)
      sigma_hat_i  <- sqrt(variance_i)
      sigmahat_vec <- c(sigmahat_vec, sigma_hat_i)
    }
  } else if (method == 'second'){
    for (i in 1:N){
      variance_i   <- 2 * sum((data[3:Tlen, i]/2 - data[2:(Tlen - 1), i] + data[1:(Tlen - 2), i]/2)^2)/(3 * Tlen - 6)
      sigma_hat_i  <- sqrt(variance_i)
      sigmahat_vec <- c(sigmahat_vec, sigma_hat_i)
    }
  } else {
    sigmahat_vec <- NULL
  }
  return(sigmahat_vec)
}

#Create a matrix (for size and power table for example) and write them in the tex file
creating_matrix_and_texing <- function(vect, vect_t, vect_alpha, filename){
  matrix_ <- matrix(vect, nrow = length(vect_t), ncol = length(vect_alpha), byrow = TRUE)
  rownames(matrix_) <- vect_t
  colnames(matrix_) <- vect_alpha
  
  addtorow     <- list()
  addtorow$pos <- list(0, 0)
  addtorow$command <- c("& \\multicolumn{3}{c}{nominal size $\\alpha$} \\\\\n",
                        "$T$ & 0.01 & 0.05 & 0.1 \\\\\n") 
  print.xtable(xtable(matrix_, digits = c(3), align = "cccc"), type = "latex",
               file = filename, add.to.row = addtorow, include.colnames = FALSE)
}

#Truncated functions for speeding up the process
statistics <- function(data, sigma_vec = 1, n_ts = 2, grid = NULL,
                       ijset = NULL, alpha = 0.05, sim_runs = 1000) {
  
  t_len <- nrow(data)
  
  #If grid is not supplied, we construct it by default
  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }
  
  #If ijset is not supplied, we compare all
  #possible pairs of time series.
  if (is.null(ijset)) {
    ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
    ijset <- ijset[ijset$i < ijset$j, ]
  }
  
  psi   <- compute_statistics(data = data, sigma = 1, sigma_vec = sigma_vec,
                              n_ts = n_ts, grid = grid, deriv_order = 0,
                              epidem = FALSE)
  stat  <- psi$stat
  gset_with_values <- psi$gset_with_values
  
  return(list(stat = stat, stat_pairwise = psi$stat_pairwise,
              ijset = ijset, gset_with_values = gset_with_values))
}

statistics_full <- function(data, sigma_vec = 1, n_ts = 2, grid = NULL,
                            ijset = NULL, alpha = c(0.05), sim_runs = 1000) {
  
  t_len <- nrow(data)
  
  #If grid is not supplied, we construct it by default
  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }
  
  #If ijset is not supplied, we compare all
  #possible pairs of time series.
  if (is.null(ijset)) {
    ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
    ijset <- ijset[ijset$i < ijset$j, ]
  }
  
  # Select (1-alpha) quantile of the multiscale statistic under the null
  quantiles <- compute_quantiles(t_len = t_len, grid = grid, n_ts = n_ts,
                                 ijset = ijset, sigma = 1,
                                 sim_runs = sim_runs,
                                 deriv_order = 0,
                                 correction = TRUE, epidem = FALSE)
  
  probs <- as.vector(quantiles$quant[1, ])
  quant <- as.vector(quantiles$quant[2, ])
  
  quant_vec <- c()
  for (alpha_ind in alpha){
    if (sum(probs == (1 - alpha_ind)) == 0)
      pos <- which.min(abs(probs - (1 - alpha_ind)))
    if (sum(probs == (1 - alpha_ind)) != 0)
      pos <- which.max(probs == (1 - alpha_ind))    
    quant_vec <- c(quant_vec, quant[pos])
  }
  
  psi   <- compute_statistics(data = data, sigma = 1, 
                              sigma_vec = sigma_vec, n_ts = n_ts,
                              grid = grid, deriv_order = 0,
                              epidem = FALSE)
  stat  <- psi$stat
  
  return(list(quant = quant_vec, stat = stat, stat_pairwise = psi$stat_pairwise,
              ijset = ijset, sim_runs = sim_runs, grid = grid))
}

#Local linear estimator of the trend function 
#using the rectangular kernel. 
nadaraya_watson_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    result = sum((abs((grid_p - u) / bw) <= 1) * data_p)
    norm = sum((abs((grid_p - u) / bw) <= 1))
    return(result/norm)
  }
}


produce_smoothed_plots <- function(matrix, pdfname, y_min, y_max, ticks_at,
                                   ticks_labels, yaxp_){
  t_len <- nrow(matrix)
  grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
  
  pdf(pdfname, width=10, height=10, paper="special")
  par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
  par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins
  
  for (h in c(0.05, 0.1, 0.15)){
    plot(NA, ylab = "", xlab = "", xlim = c(0,1), ylim = c(y_min, y_max),
         yaxp = yaxp_, xaxt = 'n', mgp = c(2,0.5,0), cex = 1.2, tck = -0.025)
    #i <- 1
    for (column in colnames(matrix)){
      smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                               MoreArgs = list(matrix[, column], grid_points, h))
      lines(grid_points, smoothed_curve)#, col = i)
      #  i <- i + 1
    }
    #axis(2, at = y_ticks_at)
    
    if (h == 0.15) {axis(1, at = grid_points[ticks_at], labels = ticks_labels)}
    #else {axis(1, at = grid_points[seq(5, 125, by = 20)], labels = NA)}
    legend("bottomright", inset = 0.01, legend = c(paste0("h = ", h)), lty = 1,
           cex = 0.95, ncol = 1)
    #legend("topright", inset = 0.01, legend = colnames(matrix), lty = 1,
    #       cex = 0.95, ncol = 1, col = 1:ncol(matrix))
  }
  
  dev.off()
}

# pdfname <- "output/plots/gdp/smoothed_gdp_data.pdf"
# produce_smoothed_plots(gdp_mat_growth, pdfname,
#                        y_min = min(gdp_mat_growth) + 0.035,
#                        y_max = max(gdp_mat_growth) - 0.02,
#                        ticks_at =  at, ticks_labels = dates[at],
#                        yaxp_ = c(-0.02, 0.02, 4))
# 
# pdfname_augm <- "output/plots/gdp/smoothed_gdp_data_augmented.pdf"
# produce_smoothed_plots(gdp_mat_augm, pdfname_augm,
#                        y_min = min(gdp_mat_augm) + 0.03,
#                        y_max = max(gdp_mat_augm) - 0.02,
#                        ticks_at =  at, ticks_labels = dates[at],
#                        yaxp_ = c(-0.03, 0.01, 4))


# ############################
# #Plots for the presentation#
# ############################
# 
# #Original time series
# 
# country1 <- which(colnames(gdp_mat_original) == "AUT")
# country2 <- which(colnames(gdp_mat_original) == "DEU")
# 
# pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], ".pdf"),
#     width=5.5, height=3, paper="special")
# 
# par(cex = 1, tck = -0.025)
# par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# 
# plot(x = dates, y = gdp_mat_original[, country1],
#      ylim=c(min(gdp_mat_original[, country1], gdp_mat_original[, country2]),
#             max(gdp_mat_original[, country1], gdp_mat_original[, country2] + 0.03)), type="l",
#      col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
# lines(x = dates, y = gdp_mat_original[, country2], col="#604c38")
# legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
#        col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
# dev.off()
# 
# 
# pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], "_1.pdf"),
#     width=5.5, height=3, paper="special")
# 
# par(cex = 1, tck = -0.025)
# par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# 
# plot(x = dates, y = gdp_mat_original[, country1],
#      ylim=c(min(gdp_mat_original[, country1], gdp_mat_original[, country2]),
#             max(gdp_mat_original[, country1], gdp_mat_original[, country2] + 0.03)), type="l",
#      col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
# rect(xleft=as.Date("1988-04-01", format = "%Y-%m-%d"),
#      xright = as.Date("1994-10-01", format = "%Y-%m-%d"),
#      ybottom=par("usr")[3], ytop=par("usr")[4],
#      density=40, col = "#d3d3d3", border = NA)
# lines(x = dates, y = gdp_mat_original[, country1], col="#EB811B")
# lines(x = dates, y = gdp_mat_original[, country2], col="#604c38")
# segments(as.Date("1988-04-01", format = "%Y-%m-%d"), par("usr")[3],
#          as.Date("1994-10-01", format = "%Y-%m-%d"), par("usr")[3],
#          col="#604c38", lwd = 3)
# legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
#        col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
# dev.off()
# 
# 
# pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], "_2.pdf"),
#     width=5.5, height=3, paper="special")
# 
# par(cex = 1, tck = -0.025)
# par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# 
# plot(x = dates, y = gdp_mat_original[, country1],
#      ylim=c(min(gdp_mat_original[, country1], gdp_mat_original[, country2]),
#             max(gdp_mat_original[, country1], gdp_mat_original[, country2] + 0.03)), type="l",
#      col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
# rect(xleft=as.Date("1988-04-01", format = "%Y-%m-%d"),
#      xright = as.Date("1994-10-01", format = "%Y-%m-%d"),
#      ybottom=par("usr")[3], ytop=par("usr")[4],
#      density=40, col = "#d3d3d3", border = NA)
# rect(xleft=as.Date("2006-01-01", format = "%Y-%m-%d"),
#      xright = as.Date("2010-04-01", format = "%Y-%m-%d"),
#      ybottom=par("usr")[3], ytop=par("usr")[4],
#      density=40, col = "#d3d3d3", border = NA)
# lines(x = dates, y = gdp_mat_original[, country1], col="#EB811B")
# lines(x = dates, y = gdp_mat_original[, country2], col="#604c38")
# segments(as.Date("1988-04-01", format = "%Y-%m-%d"), par("usr")[3],
#          as.Date("1994-10-01", format = "%Y-%m-%d"), par("usr")[3],
#          col="#604c38", lwd = 3)
# segments(as.Date("2006-01-01", format = "%Y-%m-%d"), par("usr")[3],
#          as.Date("2010-04-01", format = "%Y-%m-%d"), par("usr")[3], col="#604c38", lwd = 3)
# legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
#        col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
# dev.off()
# 
# #Adjusted time series
# 
# country1 <- which(colnames(gdp_mat_original) == "AUT")
# country2 <- which(colnames(gdp_mat_original) == "DEU")
# 
# pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], "_adj.pdf"),
#     width=5.5, height=3, paper="special")
# 
# par(cex = 1, tck = -0.025)
# par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# 
# plot(x = dates, y = gdp_mat[, country1],
#      ylim=c(min(gdp_mat[, country1], gdp_mat[, country2]),
#             max(gdp_mat[, country1], gdp_mat[, country2] + 0.04)), type="l",
#      col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
# lines(x = dates, y = gdp_mat[, country2], col="#604c38")
# legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
#        col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
# dev.off()
# 
# #Original time series for Canada/USA
# 
# country1 <- which(colnames(gdp_mat_original) == "CAN")
# country2 <- which(colnames(gdp_mat_original) == "USA")
# 
# pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], ".pdf"),
#     width=5.5, height=3, paper="special")
# 
# par(cex = 1, tck = -0.025)
# par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# 
# plot(x = dates, y = gdp_mat_original[, country1],
#      ylim=c(min(gdp_mat_original[, country1], gdp_mat_original[, country2]),
#             max(gdp_mat_original[, country1], gdp_mat_original[, country2] + 0.03)), type="l",
#      col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
# lines(x = dates, y = gdp_mat_original[, country2], col="#604c38")
# legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
#        col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
# dev.off()
# 
# #Adjusted time series for Canada/USA
# 
# country1 <- which(colnames(gdp_mat) == "CAN")
# country2 <- which(colnames(gdp_mat) == "USA")
# 
# pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], "_adj.pdf"),
#     width=5.5, height=3, paper="special")
# 
# par(cex = 1, tck = -0.025)
# par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
# par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
# 
# plot(x = dates, y = gdp_mat[, country1],
#      ylim=c(min(gdp_mat[, country1], gdp_mat[, country2]),
#             max(gdp_mat[, country1], gdp_mat[, country2] + 0.03)), type="l",
#      col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
# lines(x = dates, y = gdp_mat[, country2], col="#604c38")
# legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
#        col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
# dev.off()

##########################
#Plotting the time series#
##########################

# pdf("plots/real_housing_prices.pdf", width=10, height=10, paper="special")
# par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
# par(mar = c(0, 0.5, 1.5, 0)) #Margins for each plot
# par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins
# 
# plot(NA, ylab="", xlab = "", xlim = c(0, t_len),
#      ylim = c(0, max(hp_data$hpreal, na.rm = TRUE)), xaxt = 'n',
#      mgp=c(2,0.5,0), cex = 1.2, tck = -0.025, main = "Real house prices")
# for (country in countries){
#   tmp <- hp_data[hp_data$iso == country, ]
#   tmp <- tmp[order(tmp$year),]
#   lines(1:t_len, tmp$hpreal)
# }
# 
# plot(NA, ylab="", xlab = "", xlim = c(0,t_len),
#      ylim = c(1, max(hp_data$log_hp, na.rm = TRUE)), xaxt = 'n',
#      mgp=c(2,0.5,0), cex = 1.2, tck = -0.025,
#      main = "Logarithm of the real house prices")
# for (country in countries){
#   tmp <- hp_data[hp_data$iso == country, ]
#   tmp <- tmp[order(tmp$year),]
#   lines(1:t_len, tmp$log_hp)
# }
# 
# plot(NA, ylab="", xlab = "", xlim = c(0,t_len), ylim = c(-0.5, 0.75),
#      xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025,
#      main = "Growth rate of the logarithm of the real house prices")
# for (country in countries){
#   tmp <- hp_data[hp_data$iso == country, ]
#   tmp <- tmp[order(tmp$year),]
#   lines(1:t_len, tmp$delta_log_hp)
# }
# axis(1, at = ticks, labels = dates[ticks])
# 
# dev.off()
#
# #Producing the smoothed curves using local linear estimator
# 
# produce_smoothed_plots(matrix = hp_log,
#                        pdfname = "plots/smoothed_hp_data.pdf",
#                        y_min = min(hp_log) + 0.02, y_max = max(hp_log) - 0.01,
#                        ticks_at =  ticks, ticks_labels = dates[ticks],
#                        yaxp_ = c(0, 6, 6))
# 
# produce_smoothed_plots(matrix = hp_log_augm,
#                        pdfname = "plots/smoothed_hp_data_augmented.pdf",
#                        y_min = min(hp_log_augm) + 0.01,
#                        y_max = max(hp_log_augm) - 0.01,
#                        ticks_at =  ticks, ticks_labels = dates[ticks],
#                        yaxp_ = c(-3, 3, 6))



# ################################################################################
# #################  GROWTH RATE of the log of house prices  #####################
# ################################################################################
# 
# hp_growth_rate           <- matrix(NA, ncol = n_ts, nrow = t_len)
# colnames(hp_growth_rate) <- countries
# 
# hp_growth_rate_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
# colnames(hp_growth_rate_augm) <- countries
# 
# beta_growth_rate  <- matrix(data = NA, nrow = 4, ncol = n_ts)
# alpha_growth_rate <- c()
# 
# for (country in countries){
#   tmp <- hp_data[hp_data$iso == country, ]
#   tmp <- tmp[order(tmp$year),]
#   
#   #Calculating first difference of the growth rate
#   tmp <-
#     tmp %>%
#     mutate(delta_log_hp = coalesce(delta_log_hp, 0)) %>%
#     mutate(delta_log_gdp = coalesce(delta_log_gdp, 0)) %>%
#     mutate(delta_ltrate = coalesce(delta_ltrate, 0)) %>%
#     mutate(delta_log_pop = coalesce(delta_log_pop, 0)) %>%
#     mutate(delta_infl = coalesce(delta_infl, 0)) %>%
#     mutate(delta_delta_hp = delta_log_hp - dplyr::lag(delta_log_hp, n = 1, default = NA))%>%
#     mutate(delta_delta_gdp = delta_log_gdp - dplyr::lag(delta_log_gdp, n = 1, default = NA))%>%
#     mutate(delta_delta_pop = delta_log_pop - dplyr::lag(delta_log_pop, n = 1, default = NA))%>%
#     mutate(delta_delta_rate = delta_ltrate - dplyr::lag(delta_ltrate, n = 1, default = NA)) %>%
#     mutate(delta_delta_infl = delta_infl - dplyr::lag(delta_infl, n = 1, default = NA))
#   
#   #Estimating beta_i
#   y_vec_tmp <- as.matrix(tmp[-1, 'delta_log_hp'])
#   x_mat_tmp <- as.matrix(tmp[-1, c('delta_log_gdp', 'delta_log_pop',
#                                    'delta_ltrate', 'delta_infl')])
#   
#   y_vec_growth_rate_tmp <- as.matrix(tmp[-1, 'delta_delta_hp'])
#   x_mat_growth_rate_tmp <- as.matrix(tmp[-1 , c('delta_delta_gdp', 'delta_delta_pop',
#                                                 'delta_delta_rate', 'delta_delta_infl')])
#   
#   
#   beta_tmp  <- solve(t(x_mat_tmp) %*% x_mat_tmp) %*% t(x_mat_tmp) %*% y_vec_tmp
#   beta_log[, i] <- beta_tmp
#   
#   beta_growth_rate_tmp  <- solve(t(x_mat_growth_rate_tmp) %*% x_mat_growth_rate_tmp) %*% t(x_mat_growth_rate_tmp) %*% y_vec_growth_rate_tmp
#   beta_growth_rate[, i] <- beta_growth_rate_tmp
#   
#   #Estimating alpha_i
#   alpha_growth_rate_tmp <- mean(tmp$delta_log_hp - as.vector(as.matrix(tmp[, c('delta_log_gdp',
#                                                                                'delta_log_pop',
#                                                                                'delta_ltrate',
#                                                                                'delta_infl')]) %*% beta_growth_rate_tmp))
#   alpha_growth_rate[i]  <- alpha_growth_rate_tmp
#   
#   alpha_tmp     <- mean(tmp$log_hp - as.vector(as.matrix(tmp[, c('log_gdp',
#                                                                  'log_pop',
#                                                                  'ltrate',
#                                                                  'infl')]) %*% beta_tmp))
#   alpha_log[i]  <- alpha_tmp
#   
#   
#   #Working with adjusted time series and storing the original one
#   y_vec_growth_rate_adj    <- tmp$delta_log_hp - as.vector(as.matrix(tmp[, c('delta_log_gdp',
#                                                                              'delta_log_pop',
#                                                                              'delta_ltrate',
#                                                                              'delta_infl')]) %*% beta_growth_rate_tmp) - alpha_growth_rate_tmp
#   hp_growth_rate[, i]      <- tmp$delta_log_hp
#   hp_growth_rate_augm[, i] <- y_vec_growth_rate_adj
#   
#   y_vec_adj        <- tmp$log_hp - as.vector(as.matrix(tmp[, c('log_gdp',
#                                                                'log_pop',
#                                                                'ltrate',
#                                                                'infl')]) %*% beta_tmp) - alpha_tmp
#   hp_log[, i]      <- tmp$log_hp
#   hp_log_augm[, i] <- y_vec_adj
#   i = i + 1
# }


# #####################
# #Estimating variance#
# #####################
# 
# #Order selection
# q_vec <- 10:20
# r_vec <- 10:15
# order_results_2 <- c()
# 
# for (j in 1:n_ts){
#   criterion_matrix <- expand.grid(q = q_vec, r = r_vec)
#   
#   criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
#   
#   for (i in 1:nrow(criterion_matrix)){
#     FPE <- c()
#     AIC <- c()
#     AICC <- c()
#     SIC <- c()
#     HQ <- c()
#     
#     different_orders <- (1:9)
#     
#     for (order in different_orders){
#       AR.struc      <- estimate_lrv(data = hp_growth_rate_augm[, j],
#                                     q = criterion_matrix$q[[i]],
#                                     r_bar = criterion_matrix$r[[i]], p = order)
#       sigma_eta_hat <- sqrt(AR.struc$vareta)
#       FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
#       AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
#       AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order + 2)/t_len))
#       SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
#       HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
#     }
#     criterion_matrix$FPE[[i]]  <- which.min(FPE)
#     criterion_matrix$AIC[[i]]  <- which.min(AIC)
#     criterion_matrix$AICC[[i]] <- which.min(AICC)
#     criterion_matrix$SIC[[i]]  <- which.min(SIC)
#     criterion_matrix$HQ[[i]]   <- which.min(HQ)
#   }
#   maxim <- max(criterion_matrix[, 3:7])
#   order_results_2 <- c(order_results_2, maxim)
#   cat("For the country ", colnames(hp_growth_rate_augm)[j],
#       "for the growth rate of house prices the results are as follows: ",
#       max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ",
#       max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ",
#       max(criterion_matrix$HQ), " \n")
# }
# 
# #Calculating each sigma_i separately
# sigmahat_vector_2 <- c()
# for (i in 1:n_ts){
#   AR.struc        <- estimate_lrv(data = hp_log_augm[, i], q = q, r_bar = r,
#                                   #p = order_results_2[i])  
#                                   p = 1)
#   sigma_hat_i     <- sqrt(AR.struc$lrv)
#   sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_hat_i)
# }
# 
# 
# #########
# #Testing#
# #########
# 
# #Constructing the grid
# u_grid      <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
# h_grid      <- seq(from = 5 / t_len, to = 1 / 4, by = 1 / t_len)
# h_grid      <- h_grid[h_grid > log(t_len) / t_len]
# grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
# 
# result <- multiscale_test(data = hp_growth_rate_augm,
#                           sigma_vec = sigmahat_vector_2, alpha = alpha,
#                           n_ts = n_ts, grid = grid, sim_runs = sim_runs,
#                           epidem = FALSE)
# 
# 
# #####################
# #Plots for the paper#
# #####################
# 
# #Producing the smoothed curves using local linear estimator
# 
# 
# produce_smoothed_plots(matrix = hp_growth_rate,
#                        pdfname = "plots/growth_rate/smoothed_hp_data.pdf",
#                        y_min = min(hp_growth_rate), y_max = max(hp_growth_rate),
#                        ticks_at =  ticks, ticks_labels = dates[ticks])
# 
# produce_smoothed_plots(matrix = hp_growth_rate_augm,
#                        pdfname = "plots/growth_rate/smoothed_hp_data_augmented.pdf",
#                        y_min = min(hp_growth_rate_augm),
#                        y_max = max(hp_growth_rate_augm),
#                        ticks_at =  ticks, ticks_labels = dates[ticks])
# 
# 
# #Producing plots with the final results
# for (l in seq_len(nrow(result$ijset))){
#   i <- result$ijset[l, 1]
#   j <- result$ijset[l, 2]
#   if (result$stat_pairwise[i, j] > result$quant){
#     produce_plots(results = result, data_i = hp_growth_rate_augm[, i],
#                   data_j = hp_growth_rate_augm[, j], at_ = ticks,
#                   labels_ = dates[ticks], name_i = countries[i],
#                   name_j = countries[j])
#   }
# }
# 
