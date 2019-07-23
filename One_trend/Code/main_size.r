#This is the main file for producing the simulation results for size, which are reported in Section 5.1.1.
rm(list=ls())

library(Rcpp)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

#The following file contains one main function that computes statistical size of different tests based on different specifications.
#All the necessary arguments to use this function are described in the file.
source("functions/CalculateSiZe.r")


##########################################
#Calculating size for small parameters a1#
##########################################
Nsim       <- 1000        # number of simulation runs for size/power calculations
sigma_eta  <- 1           # standard deviation of the innovation term in the AR model
SimRuns    <- 5000        # number of simulation runs to produce critical values
sigma.type <- 'estimated' # Estimating the long-run variance \sigma^2 or plugging the true theoretical value

different_T     <- c(250, 500, 1000)         #Sample sizes for calculating size 
different_a1    <- c(-0.5, -0.25, 0.25, 0.5) #Values of AR(1) parameters considered 
different_alpha <- c(0.01, 0.05, 0.1)        #Different significant levels

number_of_cols            <- length(different_a1) * (length(different_alpha) + 1)
size_matrix_ms            <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix_ms)  <- different_T

k <- 1
set.seed(111)
for (a1 in different_a1){
  i <- 1
  for (T in different_T){
    size <- CalculateSize(T, a1, sigma_eta, different_alpha, Nsim = Nsim, SimRuns =SimRuns, sigma.type = sigma.type, q_ = 25, remove.small.ess = FALSE)$size.ms
    size_matrix_ms[i, (k * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(k * (length(different_alpha) + 1))] <- size
    i <- i + 1
  }
  k <- k + 1
}

print.xtable(xtable(size_matrix_ms, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
  type="latex", file=paste0("plots/size_", sigma.type, "_sigma.tex"), include.colnames = FALSE)


######################################
#Calculating size for for a_1 = +-0.9#
######################################
Nsim       <- 1000        # number of simulation runs for size/power calculations
sigma_eta  <- 1           # standard deviation of the innovation term in the AR model
SimRuns    <- 5000        # number of simulation runs to produce critical values
sigma.type <- 'estimated' # Estimating the long-run variance \sigma^2 or plugging the true theoretical value

different_T     <- c(250, 500, 1000, 2000, 3000)
different_a1    <- c(-0.9, 0.9)
different_alpha <- c(0.01, 0.05, 0.1)        #Different significant levels

number_of_cols    <- length(different_a1) * (length(different_T) + 1)
size_ms           <- matrix(NA, nrow = length(different_alpha), ncol = number_of_cols)
rownames(size_ms) <- different_alpha

k <- 1
for (a1 in different_a1){
  set.seed(0)
  i <- 1
  for (T in different_T){
    result <- CalculateSize(T, a1, sigma_eta, different_alpha, Nsim = Nsim, SimRuns = SimRuns, sigma.type = sigma.type, q_ = 50, remove.small.ess = FALSE)$size.ms
    size_ms[, (k - 1) * (length(different_T) + 1) + i + 1] <- t(result)
    i <- i + 1
  }
  k <- k + 1
}

print.xtable(xtable(size_ms, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/size_persistent_", sigma.type, "_sigma.tex"), include.colnames = FALSE)


########################################
#Calculating global size for comparison#
########################################
Nsim       <- 1000        # number of simulation runs for size/power calculations
sigma_eta  <- 1           # standard deviation of the innovation term in the AR model
SimRuns    <- 5000        # number of simulation runs to produce critical values
sigma.type <- 'true'      # Estimating the long-run variance \sigma^2 or plugging the true theoretical value

different_T     <- c(250, 500, 1000)
different_a1    <- c(-0.5, 0.5)
different_alpha <- c(0.05)

number_of_cols         <- length(different_a1) * 5
size_matrix            <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix)  <- different_T

k <- 1
set.seed(0)
for (a1 in different_a1){
  #Here we are plugging the true long-run variance to make the comparison between the method fair.
  #Therefore, all the differences in size come from the methods themselves
  i <- 1
  for (T in different_T){
    size <- CalculateSize(T, a1, sigma_eta, different_alpha, Nsim = Nsim, SimRuns = SimRuns, sigma.type = sigma.type, remove.small.ess = TRUE)
    size_matrix[i, (k-1)*5 + 2] <- size$size.ms
    size_matrix[i, (k-1)*5 + 3] <- size$size.uncor
    size_matrix[i, (k-1)*5 + 4] <- size$size.rows
    size_matrix[i, (k-1)*5 + 5] <- size$size.SiZer
    i <- i + 1
  }
  k <- k + 1
}

print.xtable(xtable(size_matrix, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
             type="latex", file=paste0("plots/comparing_size.tex"), include.colnames = FALSE)


#################################################################
#Calculating rowwise size and creating parallel coordinate plots#
#################################################################
Nsim       <- 1000        # number of simulation runs for size/power calculations
sigma_eta  <- 1           # standard deviation of the innovation term in the AR model
SimRuns    <- 5000        # number of simulation runs to produce critical values
sigma.type <- 'true'      # Estimating the long-run variance \sigma^2 or plugging the true theoretical value

T               <- 1000
different_a1    <- c(-0.5)
different_alpha <- c(0.05)

set.seed(0)
for (a1 in different_a1){
  result <- CalculateSize(T, a1, sigma_eta, different_alpha,  Nsim = Nsim, SimRuns = SimRuns, sigma.type = sigma.type, remove.small.ess = TRUE)
  h.grid <- result$h.grid
  
  for (j in 1:length(different_alpha)){
    alpha <- different_alpha[j]  
    pdffilename <- paste0("plots/new/pcp_size_T_", T, "_a1_", a1*100, ".pdf")
    pdf(pdffilename, width=5.5, height=4.16, paper="special")
  
    par(mfrow = c(1, 1), cex = 1,  tck = -0.025) #Setting the layout of the graphs

    par(mar = c(3.5, 3.5, 0, 0)) #Margins for each plot
    par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins 
  
    plot(x = h.grid, y = result$size.rw.ms[[j]]*100, ylim=c(0, 18), yaxp = c(0, 15, 3), type="l", lty=1, xaxt='n',
        ylab = "Percentage (%)", xlab='bandwidth h',mgp=c(2,0.5,0))
    points(x = h.grid, y = result$size.rw.ms[[j]]*100, pch=19, cex = 0.8)
  
    lines(x = h.grid, y = result$size.rw.uncor[[j]]*100, lwd=1.5, lty = 'dashed')
    lines(x = h.grid, y = result$size.rw.rows[[j]]*100, lwd=1.5, lty = 'dotted') 
    lines(x = h.grid, y = result$size.rw.SiZer[[j]]*100, lwd=1.5)
  
    abline(h = alpha*100, lty = 'dashed', col = 'grey')
  
    axis(1, at=seq(0.05, 0.25, by = 0.05), mgp=c(1.8,0.5,0))
    legend('topleft', cex = 0.8, bty = "n", legend = c(expression(italic(T)[MS]), expression(italic(T)[UC]), expression(italic(T)[RW]), expression(italic(T)[SiZer])),
         pch = c(19, NA, NA, NA), lty = c('solid', 'dashed', 'dotted', 'solid'), y.intersp=1.25)
    dev.off()
  }
}