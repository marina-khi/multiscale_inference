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