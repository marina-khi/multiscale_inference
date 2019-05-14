library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/estimating_sigma_new.R")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/data_analysis.R")


#############################
#Point 3 in Referee Report 2#
#############################

# #Defining necessary parameters
# alpha <- 0.05 #alpha for calculating quantiles
# h     <- c(0.05, 0.1, 0.15, 0.2) #Different bandwidth for plotting. Number must be <=6 in order for the plot to be readable
# 
# test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 
# pdffilename_global = paste0("Paper/Plots/temperature_data_global.pdf") #Filename for the graph
# 
# #Loading the real data for global yearly temperature
# temperature_global  <- read.table("Shape/data/global_temp.txt", header = TRUE, skip = 16)
# yearly_tempr_global <- temperature_global[temperature_global$ANNUAL > -99, 'ANNUAL']
# T_tempr_global      <- length(yearly_tempr_global)
# 
# #Order selection for global
# q <- 10:20
# r <- 5:15
# criterion_matrix_global <- expand.grid(q = q, r = r)
# 
# criterion_matrix_global$FPE <- numeric(length = nrow(criterion_matrix_global))
# criterion_matrix_global$AIC <- numeric(length = nrow(criterion_matrix_global))
# criterion_matrix_global$AICC <- numeric(length = nrow(criterion_matrix_global))
# criterion_matrix_global$SIC <- numeric(length = nrow(criterion_matrix_global))
# criterion_matrix_global$HQ  <- numeric(length = nrow(criterion_matrix_global))
# 
# for (i in 1:nrow(criterion_matrix_global)){
#   FPE <- c()
#   AIC <- c()
#   AICC <- c()
#   SIC <- c()
#   HQ <- c()
#   different_orders <- (1:9)
#   for (order in different_orders){
#     sigma_eta_hat_method2_global <- estimating_variance_new(yearly_tempr_global, criterion_matrix_global$q[[i]], order, criterion_matrix_global$r[[i]])[[3]]
#     FPE <- c(FPE, (sigma_eta_hat_method2_global^2 * (T_tempr_global + order)) / (T_tempr_global - order))
#     AIC <- c(AIC, T_tempr_global * log(sigma_eta_hat_method2_global^2) + 2 * order)
#     AICC <- c(AICC, T_tempr_global * log(sigma_eta_hat_method2_global^2) + T_tempr_global* (1 + order / T_tempr_global)/(1 - (order +2)/T_tempr_global))
#     SIC <- c(SIC, log(sigma_eta_hat_method2_global^2) + order * log(T_tempr_global) / T_tempr_global)
#     HQ <- c(HQ, log(sigma_eta_hat_method2_global^2) + 2 * order * log(log(T_tempr_global)) / T_tempr_global)
#   }
#   criterion_matrix_global$FPE[[i]] <- which.min(FPE)
#   criterion_matrix_global$AIC[[i]] <- which.min(AIC)
#   criterion_matrix_global$AICC[[i]] <- which.min(AICC)
#   criterion_matrix_global$SIC[[i]] <- which.min(SIC)
#   criterion_matrix_global$HQ[[i]]  <- which.min(HQ)
# }
# 
# #Setting tuning parameters for testing global temperature
# p_global <- 4
# q_global <- 18
# r_global <- 10
# 
# #Data analysis
# 
# sigma_hat_global <- estimating_variance_new(yearly_tempr_global, q_global, order = p_global, r_global)[[1]]
# data_analysis_global(alpha, yearly_tempr_global, test_problem, sigma_hat_global, pdffilename_global)
# 
# #plotting fitted data
# pdf("Paper/Plots/fitting_different_AR.pdf")
# par(mfrow = c(9,2), cex = 0.8, tck = -0.025) #Setting the layout of the graphs
# par(mar = c(0, 0.5, 1.0, 0)) #Margins for each plot
# par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
# 
# for (order_AR in 1:9){
#   coefficient <- estimating_variance_new(yearly_tempr_global, q_global, order = order_AR, r_global)
#   ts.sim <- arima.sim(list(c(order_AR, 0, 0), ar = coefficient[[2]]), n=T_tempr_global, rand.gen=function(n){rnorm(n, sd=coefficient[[3]])} )
#   
#   plot(ts.sim, ylab="", xlab = "", ylim = c(-1, 1), type = 'l',axes=FALSE, frame.plot=TRUE, cex = 1.2, tck = -0.025)
#   Axis(side=2, at  = c(-0.75,-0.25, 0.25, 0.75))
#   Axis(side=1, labels=FALSE)
#   plot(yearly_tempr_global - ts.sim, ylab="", xlab = "", ylim = c(-1, 1),yaxp  = c(-0.75, 0.75, 3), type = 'l', axes=FALSE, frame.plot=TRUE, cex = 1.2, tck = -0.025)
#   Axis(side=1, labels=FALSE)
#   Axis(side=2, labels = FALSE)
# }
# 
# dev.off()

#############################
#Point 5 in Referee Report 1#
#############################
set.seed(1) #For reproducibility

#Defining necessary parameters
sigma_eta           <- 1    #Standard deviation of the innovation term
alpha               <- 0.05  #Level of significance
different_T         <- c(250) #Different lengths of time series for which we compare SiZer and our method
different_a1        <- c(-0.25, 0.25) #Different a_1 in AR(1) model
slopes_for_negative <- c(1.0, 1.25, 1.5) #Slopes for power calculations for negative a_1
slopes_for_positive <- c(2.0, 2.25, 2.5) #Slopes for power calculations for positive a_1

colorlist  <- c('red', 'purple', 'blue', 'grey')
kernel_ind <- 2

pdffilename <- "JRSSB_submission/Plots/SiZer_comparison_without_adjustment.pdf" #Path 

pdf(pdffilename)

par(mfrow = c(length(different_a1)*length(different_T)*length(slopes_for_negative),2), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
par(mar = c(1, 1, 0, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins


for (a_1 in different_a1){
  sigmahat   <- sqrt(sigma_eta^2/((1 - a_1)^2))
  if (a_1 > 0){
    slopes <- slopes_for_positive
  } else {
    slopes <- slopes_for_negative
  }
  for (slope in slopes){
    for (T_size in different_T){
      different_i <- seq(from = 5/T_size, to = 1, by = 5/T_size)
      different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
      gamma = c()
      for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
        gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
      }
        
      #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
      true_var <- gamma[1] / T_size
      for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}
        
      T_star   <- gamma[1]/true_var
      
      SiZer_matrix      <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma, a_1, sigma_eta)  
          
      line_trend  <- numeric(T_size)
      for (i in 1:T_size) {line_trend[i] = (i - 0.5*T_size) * slope/T_size}
      y_data_ar_1_with_trend <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta)) + line_trend
          
      #g_t_set <- psihat_statistic_ll(y_data_ar_1_with_trend, SiZer_matrix, kernel_ind, sigmahat)[[1]]
      


      g_t_set_temp <- NULL
      for (bandwidth in different_h){
        SiZer_matrix_temp <- subset(SiZer_matrix, h == bandwidth, select = c(u, h, values, lambda))
        gaussian_quantile <- calculating_gaussian_quantile_ll(T_size, SiZer_matrix_temp, "comparison", kernel_ind, alpha)
        g_t_set_temp_temp <- psihat_statistic_ll(y_data_ar_1_with_trend, SiZer_matrix_temp, kernel_ind, sigmahat)[[1]]
        g_t_set_temp_temp$gaussian_quantile <- gaussian_quantile
        g_t_set_temp      <- rbind(g_t_set_temp, g_t_set_temp_temp)
      }
      
      SiZer_matrix$values <- NULL
      g_t_set <- merge(SiZer_matrix, g_t_set_temp, by = c('h', 'u', 'lambda'))      
                  
      for (row in 1:nrow(g_t_set)){
        i              = g_t_set[row, 'u']
        h              = g_t_set[row, 'h']
        q_h            = g_t_set[row, 'q_h']
        sd_m_hat_prime = g_t_set[row, 'sd']
        
        XtWX_inverse_XtW   = g_t_set$XtWX_inv_XtW[[row]]
            
        if (!is.null(XtWX_inverse_XtW)) {
          m_hat_prime <- (XtWX_inverse_XtW %*% y_data_ar_1_with_trend)[2]
          
          if (m_hat_prime - q_h * sd_m_hat_prime > 0){
            g_t_set$results_their[[row]] = 1
          } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
            g_t_set$results_their[[row]] = -1
          } else {
            g_t_set$results_their[[row]] = 0
          }
        } else {
          g_t_set$results_their[[row]] = 2
        }
        
        if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
          g_t_set$results_our[[row]] = 1
        } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + gaussian_quantile){
          g_t_set$results_our[[row]] = -1
        } else {
          g_t_set$results_our[[row]] = 0
        }
      }
    }

    
        
    plot.SiZer(g_t_set$results_their, different_i, different_h, ylab=expression(log[10](h)), 
                      colorlist=colorlist, title = paste0("SiZer results, T=", T_size, ", a1 = ", a_1, ", slope = ", slope))    

    plot.SiZer(g_t_set$results_our, different_i, different_h, ylab=expression(log[10](h)), 
                      colorlist=colorlist, title = paste0("Our results, T=", T_size, ", a1 = ", a_1, ", slope = ", slope))    
  }
}
dev.off()
