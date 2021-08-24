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