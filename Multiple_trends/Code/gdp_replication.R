rm(list=ls())

library(tidyr)
library(multiscale)
library(dplyr)

##############
#Data loading#
##############

source("gdp_replication_data.R")

X_mat_filled <- 
  X_mat %>%
  subset(date <= as.Date("30-09-2010", format = "%d-%m-%Y")) %>%
  subset(WBcode %in% c("AUT", "AUS", "CAN", "CHE", "DEU",
                       "FIN", "FRA", "GBR", "JPN", "NOR",
                       "USA")) %>%
  # subset(WBcode %in% c("AUT", "AUS", "CAN", "CHE", "DEU",
  #                      "FIN", "FRA", "GBR",
  #                      "USA")) %>%  
  group_by(WBcode)  %>%
  fill(h_it, .direction = 'down')  #Extrapolating educational attainment for the last two quarters 


countries <- unique(X_mat_filled$WBcode)
dates     <- unique(X_mat_filled$date)
dates     <- dates[dates >= as.Date("01-01-1976", format = "%d-%m-%Y")]
n_ts      <- length(unique(X_mat_filled$WBcode))
t_len     <- nrow(X_mat_filled) / n_ts

#############################
#Estimating the coefficients#
#############################

gdp_mat_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(gdp_mat_augm) <- countries

gdp_mat_growth           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(gdp_mat_growth) <- countries

beta      <- matrix(data = NA, nrow = 3, ncol = n_ts)
alpha_vec <- c()
i         <- 1

for (country in countries){
  tmp <- X_mat_filled[X_mat_filled$WBcode == country, ]
  tmp <- tmp[order(tmp$date),]
  
  #Calculating first differences
  tmp <- 
    tmp %>%
    mutate(delta_k_it = log(k_it) - log(dplyr::lag(k_it, n = 1, default = NA))) %>%
    mutate(delta_h_it = log(h_it) - log(dplyr::lag(h_it, n = 1, default = NA))) %>%
    mutate(delta_gdp_it = log(gdp_it) - log(dplyr::lag(gdp_it, n = 1, default = NA))) %>%
    mutate(delta_l_it = log(l_it) - log(dplyr::lag(l_it, n = 1, default = NA))) %>%
    mutate(delta_k_it = coalesce(delta_k_it, 0)) %>%
    mutate(delta_h_it = coalesce(delta_h_it, 0)) %>%
    mutate(delta_gdp_it = coalesce(delta_gdp_it, 0)) %>%
    mutate(delta_l_it = coalesce(delta_l_it, 0)) %>%
    mutate(delta_delta_gdp = delta_gdp_it- dplyr::lag(delta_gdp_it, n = 1, default = NA))%>%
    mutate(delta_delta_h = delta_h_it- dplyr::lag(delta_h_it, n = 1, default = NA))%>%
    mutate(delta_delta_l = delta_l_it- dplyr::lag(delta_l_it, n = 1, default = NA))%>%
    mutate(delta_delta_k = delta_k_it- dplyr::lag(delta_k_it, n = 1, default = NA))

  #Estimating beta_i
  y_vec_tmp <- as.matrix(tmp[-1, 'delta_delta_gdp'])
  x_mat_tmp <- as.matrix(tmp[-1 , c('delta_delta_h', 'delta_delta_l', 'delta_delta_k')])
  beta_tmp  <- solve(t(x_mat_tmp) %*% x_mat_tmp) %*% t(x_mat_tmp) %*% y_vec_tmp
  beta[, i] <- beta_tmp
  
  #Estimating alpha_i
  alpha_tmp    <- mean(tmp$delta_gdp_it - as.vector(as.matrix(tmp[, c('delta_h_it', 'delta_l_it', 'delta_k_it')]) %*% beta_tmp))
  alpha_vec[i] <- alpha_tmp

  #Working with adjusted time series and storing the original one
  y_vec_adj             <- tmp$delta_gdp_it - as.vector(as.matrix(tmp[, c('delta_h_it', 'delta_l_it', 'delta_k_it')]) %*% beta_tmp) - alpha_tmp
  gdp_mat_growth[, i]   <- tmp$delta_gdp_it
  gdp_mat_augm[, i]     <- y_vec_adj
  i = i + 1
}

rm(tmp, y_vec_tmp, x_mat_tmp, beta_tmp, alpha_tmp)


#####################
#Estimating variance#
#####################

#Order selection
q <- 10:20
r <- 10:15
order_results <- c()

for (j in 1:n_ts){
  criterion_matrix <- expand.grid(q = q, r = r)
  
  criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
  
  for (i in 1:nrow(criterion_matrix)){
    FPE <- c()
    AIC <- c()
    AICC <- c()
    SIC <- c()
    HQ <- c()
    
    different_orders <- (1:9)
    
    for (order in different_orders){
      AR.struc      <- estimate_lrv(data = gdp_mat_augm[, j], q = criterion_matrix$q[[i]],
                                    r_bar = criterion_matrix$r[[i]], p = order)
      sigma_eta_hat <- sqrt(AR.struc$vareta)
      FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
      AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
      AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order +2)/t_len))
      SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
      HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
    }
    criterion_matrix$FPE[[i]]  <- which.min(FPE)
    criterion_matrix$AIC[[i]]  <- which.min(AIC)
    criterion_matrix$AICC[[i]] <- which.min(AICC)
    criterion_matrix$SIC[[i]]  <- which.min(SIC)
    criterion_matrix$HQ[[i]]   <- which.min(HQ)
  }
  maxim <- max(criterion_matrix[, 3:7])
  order_results <- c(order_results, maxim)
  cat("For the country ", colnames(gdp_mat_augm)[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
}


#Setting tuning parameters for testing
q     <- 15
r_bar <- 10

#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 1:n_ts){
  AR.struc        <- estimate_lrv(data = gdp_mat_augm[, i], q = q, r_bar = r_bar, p = order_results[i])
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}

sigmahat_vector_2 <- rep(mean(sigmahat_vector), n_ts) #Averaging all lrvs

# #Calculating each sigma_i separately
# sigmahat_vector2 <- c()
# for (i in 1:n_ts){
#   AR.struc        <- estimate_lrv(data = gdp_mat[, i], q = q, r_bar = r_bar, p = 1)
#   sigma_hat_i     <- sqrt(AR.struc$lrv)
#   sigmahat_vector2 <- c(sigmahat_vector2, sigma_hat_i)
# }

#################################
#Parameters and necessary values#
################################$

alpha    <- 0.05
sim_runs <- 5000

#Constructing the grid
u_grid      <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
h_grid      <- seq(from = 5 / t_len, to = 1 / 4, by = 5 / t_len)
h_grid      <- h_grid[h_grid > log(t_len) / t_len]
grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

#For plotting
grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
ticks <- c(0, 40, 80, 120)

#####################
#Plots for the paper#
#####################

#Calculating smoothed curve for the data using local linear estimator#

source("functions.R")

pdf("plots/smoothed_gdp_data.pdf", width=10, height=10, paper="special")
par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins

for (h in c(0.05, 0.1, 0.15)){
  plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(-0.03, 0.03), yaxp  = c(-0.02, 0.02, 2), xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
  for (column in colnames(gdp_mat_growth)){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(gdp_mat_growth[, column], grid_points, h))
    lines(grid_points, smoothed_curve)
  }
  
  if (h == 0.15) {axis(1, at = grid_points[seq(5, 125, by = 20)], labels = dates[seq(5, 125, by = 20)])}
  #else {axis(1, at = grid_points[seq(5, 125, by = 20)], labels = NA)}
  legend("topright", inset = 0.01, legend=c(paste0("h = ", h)), lty = 1, cex = 0.95, ncol=1)
}

dev.off()

pdf("plots/smoothed_gdp_data_augmented.pdf", width=10, height=10, paper="special")
par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins

for (h in c(0.05, 0.1, 0.15)){
  plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(-0.03, 0.03), yaxp  = c(-0.02, 0.02, 2), xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
  for (column in colnames(gdp_mat_augm)){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(gdp_mat_augm[, column], grid_points, h))
    lines(grid_points, smoothed_curve)
  }
  
  if (h == 0.15) {axis(1, at = grid_points[seq(5, 125, by = 20)], labels = dates[seq(5, 125, by = 20)])}
  #else {axis(1, at = grid_points[seq(5, 125, by = 20)], labels = NA)}
  legend("topright", inset = 0.01, legend=c(paste0("h = ", h)), lty = 1, cex = 0.95, ncol=1)
}

dev.off()


#Calculating the statistic for real data
result <- multiscale_test(data = gdp_mat_augm,
                          sigma_vec = sigmahat_vector_2,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

#Rename the countries for the plots
#countries_names <- c("Australia", "Austria", "Canada", "Switzerland", "Germany",
#                     "Finland", "France", "UK", "Japan", "Norway", "USA")

for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  if (result$stat_pairwise[i, j] > result$quant){
    filename = paste0("plots/", countries[i], "_vs_", countries[j], "2.pdf")
    pdf(filename, width=5.5, height=10.5, paper="special")
    layout(matrix(c(1, 2, 3),ncol=1), widths=c(2.2, 2.2, 2.2),
           heights=c(1.5, 1.5, 1.8), TRUE)
    
    #Setting the layout of the graphs
    
    par(cex = 1, tck = -0.025)
    par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
    par(oma = c(0.2, 1.5, 2, 0.2)) #Outer margins
    
    plot(gdp_mat_augm[, i], ylim=c(min(gdp_mat_augm[, i], gdp_mat_augm[, j]),
                              max(gdp_mat_augm[, i], gdp_mat_augm[, j])),
         type="l", col="blue", ylab="", xlab="", xaxt = "n", mgp=c(1, 0.5, 0))
    lines(gdp_mat_augm[, j], col="red")
    axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
         cex.axis = 0.95, mgp=c(1, 0.5, 0))
    
    title(main = "(a) adjusted GDP", font.main = 1, line = 0.5)
    legend("topright", inset = 0.02, legend=c(countries[i], countries[j]),
           col = c("blue", "red"), lty = 1, cex = 0.95, ncol = 1)
    
    par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot
    
    #Plotting the smoothed version of the time series that we have
    smoothed_i  <- mapply(nadaraya_watson_smoothing, grid_points,
                          MoreArgs = list(gdp_mat_augm[, i], grid_points, bw = 5 / t_len))
    smoothed_j  <- mapply(nadaraya_watson_smoothing, grid_points,
                          MoreArgs = list(gdp_mat_augm[, j], grid_points, bw = 5 / t_len))
    
    plot(smoothed_i, ylim=c(min(gdp_mat_augm[, i], gdp_mat_augm[, j]),
                            max(gdp_mat_augm[, i], gdp_mat_augm[, j])), type="l",
         col="black", ylab="", xlab = "", xaxt = "n", mgp=c(1,0.5,0))
    axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
         cex.axis = 0.95, mgp=c(1, 0.5, 0))
    title(main = "(b) smoothed curves from (a)", font.main = 1, line = 0.5)
    lines(smoothed_j, col="red")
    
    par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot
    gset    <- result$gset_with_values[[l]]
    a_t_set <- subset(gset, test == TRUE, select = c(u, h))
    if (nrow(a_t_set) > 0){
      p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * t_len + 0.5,
                            'endpoint' = (a_t_set$u + a_t_set$h) * t_len - 0.5, 'values' = 0)
      p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
      
      #Produce minimal intervals
      p_t_set2  <- compute_minimal_intervals(p_t_set)
      
      plot(NA, xlim=c(0, t_len),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab="", xaxt = "n",
           mgp=c(2, 0.5, 0), yaxt = "n")
      axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
           cex.axis = 0.95, mgp=c(1, 0.5, 0))
      title(main = "(c) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
      #title(xlab = "quarter", line = 1.7, cex.lab = 0.9)
      segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint, p_t_set2$values, lwd = 2)
      segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint, p_t_set$values, col = "gray")
    } else {
      #If there are no intervals where the test rejects, we produce empty plots
      plot(NA, xlim=c(0, t_len),  ylim = c(0, 1), xlab="", ylab = "", xaxt = "n",
           mgp=c(2,0.5,0), yaxt = "n")
      axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
           cex.axis = 0.95, mgp=c(1, 0.5, 0))
      title(main = "(c) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
      #title(xlab = "quarter", line = 1.7, cex.lab = 0.9)
    }
    mtext(paste0("Comparison of ", countries[i], " and ", countries[j]), side = 3,
          line = 0, outer = TRUE, font = 1, cex = 1.2)
    dev.off()
  }
}

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
