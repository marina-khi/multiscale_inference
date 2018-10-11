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
p         <- 1 #Order of AR(p)
N_rep     <- 1000 #Number of replications for comparison of the estimates
sigma_eta <- 1 # Sqrt root of the variance of the innovation \eta_t
T_size    <- 500

different_a          <- c(-0.95,-0.75,-0.5,-0.25,0.25,0.5,0.75,0.95) # AR parameters
different_slope_facs <- c(1, 10) # slope of linear trend = slope_fac * sqrt(sigma_eta^2/(1-a_star^2))
different_q          <- c(20, 25, 30, 35) # tuning parameters for first-step estimator

K1 <- p + 1 # tuning parameters for second-step estimator
K2 <- 10
M1 <- 20  # tuning parameters for HvK estimator
M2 <- 30


#######################################################################
#Comparing three methods for different slope factors and different a_1#
#######################################################################

set.seed(1)

for (q in different_q){
  for (K2 in c(10, 15)){
    for (slope_fac in different_slope_facs){
      name_spec <- paste0("T=",T_size,"_slope=",slope_fac,"_(q,K1,K2,M1,M2)=(",q,",",K1,",",K2,",",q - 5,",",q + 5,")")
    
      a_mat2   <- a_mat_HvK <- a_mat_oracle <- numeric()
      lrv_mat2 <- lrv_mat_HvK <- lrv_mat_oracle <- numeric()
      
      for (a_1 in different_a){
        pdfname_a_hat = paste0("Paper/Plots/Robustness/a_hat_histograms_a1=", a_1*100,"_", name_spec, ".pdf")
        pdfname_lrv = paste0("Paper/Plots/Robustness/lrv_histograms_a1=", a_1*100,"_", name_spec, ".pdf")
        slope <- slope_fac * sqrt(sigma_eta^2/(1 - a_1^2))
        mse_results <- histograms_for_variance_estimators(a_1, sigma_eta, T_size, p, slope, N_rep, pdfname_a_hat, pdfname_lrv,
                                           q, K1, K2, q - 5, q + 5, produce_plots = "no")
        
        a_mat2 <- c(a_mat2, mse_results[[1]])
        a_mat_HvK <- c(a_mat_HvK, mse_results[[2]])
        a_mat_oracle <- c(a_mat_oracle, mse_results[[3]])
        lrv_mat2 <- c(lrv_mat2, mse_results[[4]])
        lrv_mat_HvK <- c(lrv_mat_HvK, mse_results[[5]])
        lrv_mat_oracle <- c(lrv_mat_oracle, mse_results[[6]])
      }
      
      pdfname = paste("Paper/Plots/Robustness/MSE_a1_",name_spec,".pdf",sep="")
      plotting_MSE_graphs(a_mat2, a_mat_HvK, a_mat_oracle, pdfname, margin_ = 3, legend_position = "topright", 
                          ylab_ = "MSE", legend_ = c(expression(widehat(a)), expression(widehat(a)[HvK]), expression(widehat(a)[oracle])),
                          title_ = bquote("q = " *.(q)* ", r = " *.(K1)* ", (" *m[1]* ", " *m[2]* ") = (" *.(q - 5)* "," *.(q + 5)* ")"),
                          different_a)
      
      lrv_mat2 <- log(lrv_mat2)
      lrv_mat_HvK <- log(lrv_mat_HvK)
      lrv_mat_oracle <- log(lrv_mat_oracle)
    
      pdfname = paste("Paper/Plots/Robustness/MSE_lrv_",name_spec,".pdf",sep="")
      plotting_MSE_graphs(lrv_mat2, lrv_mat_HvK, lrv_mat_oracle, pdfname, margin_ = 3, legend_position = "topleft", ylab_ = "log(MSE)",
                          legend_ = c(expression(widehat(sigma)^2), expression(widehat(sigma)[HvK]^2), expression(widehat(sigma)[oracle]^2)),
                          title_ = bquote("q = " *.(q)* ", r = " *.(K1)* ", (" *m[1]* ", " *m[2]* ") = (" *.(q - 5)* "," *.(q + 5)* ")"),
                          different_a)
    }
  }
}