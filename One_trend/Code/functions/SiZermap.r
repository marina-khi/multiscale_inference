

SiZermap <- function(u.grid, h.grid, test.results, plot.title)

{  # computes SiZer map from the test results 

   col.vec <- c("red", "purple", "blue", "gray") 
   temp    <- sort(unique(as.vector(test.results))) + 2
   temp    <- seq(min(temp),max(temp),by=1)
   col.vec <- col.vec[temp]

   image(x=u.grid, y=log(h.grid,10), z=t(test.results), col=col.vec, xlab="u", ylab=expression(log[10](h)), main = plot.title)
}