library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("functions/functions.R")
source("functions/long_run_variance.r")
source("functions/simulations_variance.R")
source("functions/sim.r")
source("functions/estimating_sigma.R")
dyn.load("functions/estimating_sigma.dll")

###############################
#Defining necessary parameters#
###############################
p          <- 1 #Order of AR(p)
N_rep      <- 1000 #Number of replications for comparison of the estimates
sigma_eta  <- 1 # Sqrt root of the variance of the innovation \eta_t
T_size     <- 500 #Sample size considered
sim.design <- "linear"  # trend specification: "constant", "blocks", "spike", ...

different_a          <- c(-0.95,-0.75,-0.5,-0.25,0.25,0.5,0.75,0.95) # AR parameters a_1
different_slope_facs <- c(1, 10) # slope of linear trend = slope_fac * sqrt(sigma_eta^2/(1-a_star^2))
different_q          <- c(25) # tuning parameters for first-step estimator
different_r          <- c(10) #tuning parameters for second-step estimator


#######################################################################
#Comparing three methods for different slope factors and different a_1#
#######################################################################
for (q in different_q){
  for (r in different_r){
    for (slope_fac in different_slope_facs){
      name_spec <- paste0("T=",T_size,"_slope=",slope_fac,"_(q,r,M1,M2)=(",q,",",r,",",q - 5,",",q + 5,")")
      
      set.seed(1) #This is for comparing different scenarios on the same data
      
      a_mat2   <- a_mat_HvK   <- a_mat_oracle   <- numeric()
      lrv_mat2 <- lrv_mat_HvK <- lrv_mat_oracle <- numeric()
      
      for (a_1 in different_a){
        pdfname_a_hat = paste0("plots/a_hat_histograms_a1=", a_1*100,"_", name_spec, ".pdf")
        pdfname_lrv   = paste0("plots/lrv_histograms_a1=", a_1*100,"_", name_spec, ".pdf")
        #slope <- slope_fac * sqrt(sigma_eta^2/(1 - a_1^2))

        mse_results <- estimating_variance(a_1, sigma_eta, T_size, p, slope_fac, sim.design = sim.design, N_rep, pdfname_a_hat, pdfname_lrv,
                                           q, r, q - 5, q + 5, produce_plots = "selected")
        
        a_mat2       <- c(a_mat2, mse_results[[1]])
        a_mat_HvK    <- c(a_mat_HvK, mse_results[[2]])
        a_mat_oracle <- c(a_mat_oracle, mse_results[[3]])
        lrv_mat2       <- c(lrv_mat2, mse_results[[4]])
        lrv_mat_HvK    <- c(lrv_mat_HvK, mse_results[[5]])
        lrv_mat_oracle <- c(lrv_mat_oracle, mse_results[[6]])
      }
      
      lrv_mat2       <- log(lrv_mat2) #Plotting log(MSE) for long run variance estimators
      lrv_mat_HvK    <- log(lrv_mat_HvK)
      lrv_mat_oracle <- log(lrv_mat_oracle)
      
      #We produce "zoomed in" graphs only for a_hat estimators and large slope
      if (slope_fac == 10){ 
        pdfname = paste("plots/Robustness/MSE_a1_zoomed_",name_spec,".pdf",sep="")
        plotting_MSE_graphs(a_mat2, a_mat_HvK, a_mat_oracle, pdfname, margin_ = 3, legend_position = "topright", 
                            ylab_ = "MSE", legend_ = c(expression(widehat(a)), expression(widehat(a)[HvK]), expression(widehat(a)[oracle])),
                            title_ = bquote("q = " *.(q)* ", " *bar(r)* " = " *.(r)* ", (" *m[1]* ", " *m[2]* ") = (" *.(q - 5)* "," *.(q + 5)* ")"),
                            different_a, zoomed = 'yes')
      }
      
      #Graphs for the paper, not the supplement
      if ((q == 25)&&(r == 10)){ 
         pdfname = paste("plots/MSE_a1_",name_spec,".pdf",sep="")
         plotting_MSE_graphs(a_mat2, a_mat_HvK, a_mat_oracle, pdfname, margin_ = 1.5, legend_position = "topright", 
                             ylab_ = "MSE", legend_ = c(expression(widehat(a)), expression(widehat(a)[HvK]), expression(widehat(a)[oracle])),
                             title_ = "",  different_a, zoomed = 'no')
         
         pdfname = paste("plots/MSE_lrv_",name_spec,".pdf",sep="")
         plotting_MSE_graphs(lrv_mat2, lrv_mat_HvK, lrv_mat_oracle, pdfname, margin_ = 1.5, legend_position = "topleft", ylab_ = "log(MSE)",
                             legend_ = c(expression(widehat(sigma)^2), expression(widehat(sigma)[HvK]^2), expression(widehat(sigma)[oracle]^2)),
                             title_ = "", different_a, zoomed = 'no')
      }
      
      #Graphs for robustness check (supplement)
      pdfname = paste("plots/Robustness/MSE_a1_",name_spec,".pdf",sep="")
      plotting_MSE_graphs(a_mat2, a_mat_HvK, a_mat_oracle, pdfname, margin_ = 3, legend_position = "topright", 
                          ylab_ = "MSE", legend_ = c(expression(widehat(a)), expression(widehat(a)[HvK]), expression(widehat(a)[oracle])),
                          title_ = bquote("q = " *.(q)* ", " *bar(r)* " = " *.(r)* ", (" *m[1]* ", " *m[2]* ") = (" *.(q - 5)* "," *.(q + 5)* ")"),
                          different_a, zoomed = 'no')
      
      pdfname = paste("plots/Robustness/MSE_lrv_",name_spec,".pdf",sep="")
      plotting_MSE_graphs(lrv_mat2, lrv_mat_HvK, lrv_mat_oracle, pdfname, margin_ = 3, legend_position = "topleft", ylab_ = "log(MSE)",
                          legend_ = c(expression(widehat(sigma)^2), expression(widehat(sigma)[HvK]^2), expression(widehat(sigma)[oracle]^2)),
                          title_ = bquote("q = " *.(q)* ", " *bar(r)* " = " *.(r)* ", (" *m[1]* ", " *m[2]* ") = (" *.(q - 5)* "," *.(q + 5)* ")"),
                          different_a, zoomed = 'no')
    }
  }
}