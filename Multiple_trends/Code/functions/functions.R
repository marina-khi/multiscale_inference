add.quarters <- function(n, date_) {
  seq(date_, by = paste (3 * n, "months"), length = 2)[2]
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
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
s_t_1 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

#' Function needed for local linear smoothing
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
s_t_2 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

#' Function needed for local linear smoothing
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
s_t_0 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

#Local Linear estimator using the epanechnikov kernel. 
local_linear_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    s_t_2_value = s_t_2(u, bw, T_size)
    s_t_1_value = s_t_1(u, bw, T_size) 
    for (i in 1:T_size){
      k = (s_t_2_value - s_t_1_value * ((grid_p[i] - u) / bw)) * epanechnikov_kernel((grid_p[i] - u) / bw)
      result = result + k * data_p[i]
      norm = norm + k 
    }
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
    for (column in colnames(matrix)){
      smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                               MoreArgs = list(matrix[, column], grid_points, h))
      lines(grid_points, smoothed_curve)
    }
    #axis(2, at = y_ticks_at)
    
    if (h == 0.15) {axis(1, at = grid_points[ticks_at], labels = ticks_labels)}
    #else {axis(1, at = grid_points[seq(5, 125, by = 20)], labels = NA)}
    legend("bottomright", inset = 0.01, legend = c(paste0("h = ", h)), lty = 1,
           cex = 0.95, ncol=1)
  }
  
  dev.off()
}


produce_plots <- function(results, data_i, data_j, ticks_, labels_,
                          name_i, name_j){
  filename <- paste0("plots/", name_i, "_vs_", name_j, ".pdf")
  t_len <- length(data_i)
  grid_points <- seq(1/t_len, 1, by = 1/t_len)
  
  pdf(filename, width=5.5, height=10.5, paper="special")
  layout(matrix(c(1, 2, 3),ncol=1), widths=c(2.2, 2.2, 2.2),
         heights=c(1.5, 1.5, 1.8), TRUE)
    
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 2, 0.2)) #Outer margins
  
  if ((name_i %in% c("FRA", "USA")) & (name_j %in% c("FRA", "USA"))) {
    shift <- 0
  } else {
    shift <- 0.012
  }
    
  plot(data_i, ylim = c(min(data_i, data_j), max(data_i, data_j) + shift),
       type="l", col = "blue", ylab = "", xlab="", xaxt = "n",
       mgp = c(1, 0.5, 0))
  lines(data_j, col = "red")
  axis(side = 1, at = ticks_, cex.axis = 0.95, mgp = c(1, 0.5, 0), 
       labels = labels_)
    
  title(main = "(a) adjusted GDP growth rate", font.main = 1, line = 0.5)
  legend("topright", inset = 0.01, legend=c(name_i, name_j),
         col = c("blue", "red"), lty = 1, cex = 0.95, ncol = 2)
    
  par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot
    
  #Plotting the smoothed version of the time series that we have
  smoothed_i  <- mapply(nadaraya_watson_smoothing, grid_points,
                        MoreArgs = list(data_i, grid_points, bw = 0.1))
  smoothed_j  <- mapply(nadaraya_watson_smoothing, grid_points,
                        MoreArgs = list(data_j, grid_points, bw = 0.1))
    
  plot(smoothed_i, ylim = c(min(data_i, data_j), max(data_i, data_j)),
       type = "l", col = "black", ylab = "", xlab = "", xaxt = "n",
       mgp=c(1,0.5,0))
  axis(side = 1, at = ticks_, labels = labels_, cex.axis = 0.95,
       mgp=c(1, 0.5, 0))
  title(main = "(b) smoothed curves from (a)", font.main = 1, line = 0.5)
  lines(smoothed_j, col = "red")
    
  par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot
  gset    <- result$gset_with_values[[l]]
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
                 type="latex", file=paste0("plots/", name_i, "_vs_", name_j, ".tex"),
                 include.colnames = FALSE)
    print.xtable(xtable(p_t_set2_tex[order(p_t_set2_tex$from), ], digits = c(0),
                        align = paste(replicate(3, "c"), collapse = "")),
                 type="latex", file=paste0("plots/", name_i, "_vs_", name_j, "_min_intervals.tex"),
                 include.colnames = FALSE)
    
  } else {
    #If there are no intervals where the test rejects, we produce empty plots
    plot(NA, xlim=c(0, t_len),  ylim = c(0, 1), xlab="", ylab = "", xaxt = "n",
         mgp=c(2,0.5,0), yaxt = "n")
    axis(side = 1, at = ticks_, labels = labels_, cex.axis = 0.95,
         mgp = c(1, 0.5, 0))
    title(main = "(c) (minimal) intervals produced by our test", font.main = 1,
          line = 0.5)
    #title(xlab = "quarter", line = 1.7, cex.lab = 0.9)
  }
  mtext(paste0("Comparison of ", name_i, " and ", name_j), side = 3,
        line = 0, outer = TRUE, font = 1, cex = 1.2)
  dev.off()
}


produce_plots_talk <- function(results, l, data_i, data_j,
                               dates, name_i, name_j, filename){
  Tlen <- length(data_i)
  gset <- results$gset_with_values[[l]]
  
  pdf(filename, width = 5, height = 6.5, paper="special")
  layout(matrix(c(1, 2), ncol=1), widths=c(2.4, 2.4),
         heights=c(1.5, 1.8), TRUE)
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
  
  
  
  plot(x = dates, y = data_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
       col = "#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
  lines(x = dates, y = data_j, col="#604c38")
  title(main = "(a) observed daily exchange rates", font.main = 1, line = 0.5)
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
                             at_, labels_, name_i, name_j){
  filename    <- paste0("plots/hp_", name_i, "_vs_", name_j, ".pdf")
  t_len       <- length(data_i)
  grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
  
  pdf(filename, width=5.5, height=10.5, paper="special")
  layout(matrix(c(1, 2, 3),ncol=1), widths=c(2.2, 2.2, 2.2),
         heights=c(1.5, 1.5, 1.8), TRUE)
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 2, 0.2)) #Outer margins
  
  if ((name_i %in% c("AUS", "NLD")) & (name_j %in% c("AUS", "NLD"))) {
    shift <- 0.4
  } else {
    shift <- 0
  }
  
  plot(data_i, ylim = c(min(data_i, data_j), max(data_i, data_j) + shift),
       type = "l", col = "blue", ylab = "", xlab = "", xaxt = "n",
       mgp = c(1, 0.5, 0))
  lines(data_j, col="red")
  axis(side = 1, at = at_, labels = labels_,
       cex.axis = 0.95, mgp=c(1, 0.5, 0))
  
  title(main = "(a) adjusted log of housing prices", font.main = 1, line = 0.5)
  legend("topright", inset = 0.02, legend = c(name_i, name_j),
         col = c("blue", "red"), lty = 1, cex = 0.95, ncol = 2)
  
  par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot
  
  #Plotting the smoothed version of the time series that we have
  smoothed_i  <- mapply(nadaraya_watson_smoothing, grid_points,
                        MoreArgs = list(data_i, grid_points, bw = 10/t_len))
  smoothed_j  <- mapply(nadaraya_watson_smoothing, grid_points,
                        MoreArgs = list(data_j, grid_points, bw = 10/t_len))
  
  plot(smoothed_i, ylim = c(min(data_i, data_j), max(data_i, data_j)), type = "l",
       col="blue", ylab = "", xlab = "", xaxt = "n", mgp = c(1,0.5,0))
  axis(side = 1, at = at_, labels = labels_, cex.axis = 0.95,
       mgp = c(1, 0.5, 0))
  title(main = "(b) smoothed curves from (a)", font.main = 1, line = 0.5)
  lines(smoothed_j, col="red")
  
  par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot
  gset    <- result$gset_with_values[[l]]
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
    #title(xlab = "quarter", line = 1.7, cex.lab = 0.9)
    segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint,
             p_t_set2$values, lwd = 2)
    segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint,
             p_t_set$values, col = "gray")
  } else {
    #If there are no intervals where the test rejects, we produce empty plots
    plot(NA, xlim = c(0, t_len),  ylim = c(0, 1), xlab = "", ylab = "",
         xaxt = "n", mgp = c(2,0.5,0), yaxt = "n")
    axis(side = 1, at = at_, labels = labels_, cex.axis = 0.95,
         mgp = c(1, 0.5, 0))
    title(main = "(c) (minimal) intervals produced by our test", font.main = 1,
          line = 0.5)
  }
  mtext(paste0("Comparison of ", name_i, " and ", name_j), side = 3, line = 0,
        outer = TRUE, font = 1, cex = 1.2)
  dev.off()
}


produce_plots_talk <- function(results, l, data_i, data_j,
                               dates, name_i, name_j, filename){
  Tlen <- length(data_i)
  gset <- results$gset_with_values[[l]]
  
  pdf(filename, width = 5, height = 6.5, paper="special")
  layout(matrix(c(1, 2), ncol=1), widths=c(2.4, 2.4),
         heights=c(1.5, 1.8), TRUE)
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
  
  plot(x = dates, y = data_i, ylim=c(min(data_i, data_j), max(data_i, data_j)), type="l",
       col = "#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
  lines(x = dates, y = data_j, col="#604c38")
  title(main = "(a) observed daily exchange rates", font.main = 1, line = 0.5)
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
