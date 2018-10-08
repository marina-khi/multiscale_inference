library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/simulations_variance.R")
source("Shape/C_code/estimating_sigma.R")
source("Shape/estimating_sigma_new.R")
dyn.load("Shape/C_code/estimating_sigma.dll")


###############################
#Defining necessary parameters#
###############################
p                <- 1 #Order of AR(p)
N_rep            <- 1000 #Number of replications for comparison of the estimates
sigma_eta        <- 1 # Sqrt root of the variance of the innovation \eta_t
different_a      <- c(-0.95,-0.75,-0.5,-0.25,0.25,0.5,0.75,0.95) # AR parameters
different_slope_facs <- c(1, 10) # slope of linear trend = slope_fac * sqrt(1/(1-a_star^2))

T_size <- 500

L1 <- 25 # tuning parameters for first-step estimator
L2 <- 25
K1 <- p + 4 # tuning parameters for second-step estimator
K2 <- 10
M1 <- 20  # tuning parameters for HvK estimator
M2 <- 30


#################################################
#Calculating histograms for different values a_1#
#################################################

set.seed(1)

for (slope_fac in different_slope_facs){
  
  name_spec <- paste0("T=",T_size,"_slope=",slope_fac*10,"_(L1,L2,K1,K2,M1,M2)=(",L1,",",L2,",",K1,",",K2,",",M1,",",M2,")")

  a_mat2 <- a_mat_HvK <- a_mat_oracle <- numeric()
  lrv_mat2 <- lrv_mat_HvK <- lrv_mat_oracle <- numeric()
  
  for (a_1 in different_a){
    pdfname_a_hat = paste0("Paper/Plots/a_hat_histograms_a1=", a_1*100,"_", name_spec, ".pdf")
    pdfname_lrv = paste0("Paper/Plots/lrv_histograms_a1=", a_1*100,"_", name_spec, ".pdf")
    slope <- slope_fac * sqrt(1/(1 - a_1^2))
    mse_results <- histograms_for_variance_estimators(a_1, sigma_eta, T_size, p, slope, N_rep, pdfname_a_hat, pdfname_lrv,
                                       L1, L2, K1, K2, M1, M2)
    
    a_mat2 <- c(a_mat2, mse_results[[1]])
    a_mat_HvK <- c(a_mat_HvK, mse_results[[2]])
    a_mat_oracle <- c(a_mat_oracle, mse_results[[3]])
    lrv_mat2 <- c(lrv_mat2, mse_results[[4]])
    lrv_mat_HvK <- c(lrv_mat_HvK, mse_results[[5]])
    lrv_mat_oracle <- c(lrv_mat_oracle, mse_results[[6]])
  }
  
  a_mat_min <- min(c(a_mat2,a_mat_HvK,a_mat_oracle))
  a_mat_max <- max(c(a_mat2,a_mat_HvK,a_mat_oracle))
  
  pdfname = paste("Paper/Plots/MSE_a1_",name_spec,".pdf",sep="")
  pdf(pdfname, width=5.5, height=4.16, paper="special")
  par(mar = c(4, 4, 1.5, 0)) #Margins for each plot
  par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins 
  plot(a_mat2,ylim=c(a_mat_min,a_mat_max),type="l",lty=1,xaxt='n',ylab="MSE",xlab=expression(a[1]))
  points(a_mat2,pch=19)
  lines(a_mat_HvK,lty="dashed",lwd=1.5)
  points(a_mat_HvK,pch=19)
  lines(a_mat_oracle,lty="dotted",lwd=1.5)
  points(a_mat_oracle,pch=19)
  axis(1, at=1:length(different_a), labels=different_a)
  legend( "topright", cex = 1, bty = "n", legend = c(expression(widehat(a)), expression(widehat(a)[HvK]), expression(widehat(a)[oracle])), lty=c("solid","dashed","dotted"),lwd=1.5,y.intersp=1.25)
  dev.off()
  
  lrv_mat2 <- log(lrv_mat2)
  lrv_mat_HvK <- log(lrv_mat_HvK)
  lrv_mat_oracle <- log(lrv_mat_oracle)
  
  lrv_mat_min <- min(c(lrv_mat2,lrv_mat_HvK,lrv_mat_oracle))
  lrv_mat_max <- max(c(lrv_mat2,lrv_mat_HvK,lrv_mat_oracle))
  
  pdfname = paste("Paper/Plots/MSE_lrv_",name_spec,".pdf",sep="")
  pdf(pdfname, width=5.5, height=4.16, paper="special")
  par(mar = c(4, 4, 1.5, 0)) #Margins for each plot
  par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins 
  plot(lrv_mat2,ylim=c(lrv_mat_min,lrv_mat_max),type="l",lty=1,xaxt='n',ylab="log(MSE)",xlab=expression(a[1]))
  points(lrv_mat2,pch=19)
  lines(lrv_mat_HvK,lty="dashed",lwd=1.5)
  points(lrv_mat_HvK,pch=19)
  lines(lrv_mat_oracle,lty="dotted",lwd=1.5)
  points(lrv_mat_oracle,pch=19)
  axis(1, at=1:length(different_a), labels=different_a)
  legend( "topleft", cex = 1, bty = "n", legend = c(expression(widehat(sigma)^2), expression(widehat(sigma)[HvK]^2), expression(widehat(sigma)[oracle]^2)), lty=c("solid","dashed","dotted"),lwd=1.5,y.intersp=1.25)
  dev.off()
}
