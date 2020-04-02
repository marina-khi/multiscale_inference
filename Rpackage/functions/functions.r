#Function that plots the graph with MSE of the lrv/coefficient a_1 estimators for different model specifications.
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

#Function that plots the graph with histograms of the the lrv/coefficient a_1 estimators for different model specifications.
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


epan <- function(x){
  return(0.75*(1-x^2)*((sign(1-x^2)+1)/2))
}

SiZermap <- function(u.grid, h.grid, test.results, plot.title = NA){
  # computes SiZer map from the test results 
  
  col.vec <- c("red", "purple", "blue", "gray") 
  #col.vec <- c("#F7F7F7", "#969696", "#525252", "#636363") 
  temp    <- sort(unique(as.vector(test.results))) + 2
  temp    <- seq(min(temp),max(temp),by=1)
  col.vec <- col.vec[temp]
  
  image(x=u.grid, y=log(h.grid,10), z=t(test.results), col=col.vec, xlab = '', ylab=expression(log[10](h)), main = plot.title, xaxt = 'n', mgp=c(1,0.5,0))
}

# Epanechnikov kernel function, which is defined as f(x) = 3/4(1-x^2)
# for |x|<=1 and 0 elsewhere
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

#Functions needed for local linear smoothing
s_t_1 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

s_t_2 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

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

