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
#different_T         <- c(250) #Different lengths of time series for which we compare SiZer and our method
#different_a1        <- c(-0.25, 0, 0.25) #Different a_1 in AR(1) model
#slopes_for_negative <- c(1.0) #Slopes for power calculations for negative a_1
#slopes_for_positive <- c(2.0) #Slopes for power calculations for positive a_1

kernel_ind <- 2

#method = 'rowwise'
method = 'global'


#THIS IS FOR PLOTTING SIZER MAPS
#h_for_blocks <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
#t_for_blocks <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.4, 0.65, 0.76, 0.78, 0.81)
#
#colorlist  <- c('red', 'purple', 'blue', 'grey')
#
# pdffilename <- paste0("JRSSB_submission/Plots/SiZer_comparison_blocks_", method, ".pdf")
# 
# pdf(pdffilename)
# 
# par(mfrow = c(length(different_a1)*length(different_T)*length(slopes_for_negative),2), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
# par(mar = c(1, 1, 2, 0)) #Margins for each plot
# par(oma = c(1.5, 1.5, 3, 0.2)) #Outer margins
# 
# for (a_1 in different_a1){
#   sigmahat   <- sqrt(sigma_eta^2/((1 - a_1)^2))
#   if (a_1 > 0){
#     slopes <- slopes_for_positive
#   } else {
#     slopes <- slopes_for_negative
#   }
#   for (slope in slopes){
#     for (T_size in different_T){
#       different_i <- seq(from = 5/T_size, to = 1, by = 5/T_size)
#       different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
# 
#       trend_function  <- numeric(T_size)
#       for (i in 1:T_size) {trend_function[i] = 4 * sum(h_for_blocks*(1 + sign(i/T_size - t_for_blocks)))/2}
#       #for (i in 1:T_size) {trend_function[i] = 0}
#       #for (i in 1:T_size) {trend_function[i] = sinpi(6*i/T_size)}
# 
#       results <- comparing_us_and_Sizer(different_i, different_h, alpha, T_size, a_1, sigma_eta, sigmahat, trend_function, method)
#       plot.SiZer(results[[1]]$results_their, different_i, different_h, ylab=expression(log[10](h)), 
#                       colorlist=colorlist, title = paste0("SiZer results, T=", T_size, ", a1 = ", a_1))    
# 
#       plot.SiZer(results[[2]]$results_our, different_i, different_h, ylab=expression(log[10](h)), 
#                       colorlist=colorlist, title = paste0("Our results, T=", T_size, ", a1 = ", a_1))    
#     }
#   }
# }
# mtext(paste0('Blocks plus AR(1), ', method, ' method'), outer = TRUE, cex = 1.0)
# dev.off()

num_of_reps <- 1000
T_size      <- 250
different_i <- seq(from = 5/T_size, to = 1, by = 5/T_size)
different_h <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)

#CALCULATING SIZE OF THE TESTS
trend_function  <- numeric(T_size)
for (i in 1:T_size) {trend_function[i] = 0}

a_1      <- 0.25
sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))

testing_SiZer_null_a1_pos <- expand.grid(u = different_i, h = different_h)
testing_ours_null_a1_pos  <- expand.grid(u = different_i, h = different_h)

for (i in 1:num_of_reps){
  results <- comparing_us_and_Sizer(different_i, different_h, alpha, T_size, a_1, sigma_eta, sigmahat, trend_function, method)
  
  testing_SiZer_null_a1_pos <- cbind(testing_SiZer_null_a1_pos, results[[1]]$results_their)
  colnames(testing_SiZer_null_a1_pos)[i+2] <- paste0('test', i)

  testing_ours_null_a1_pos  <- cbind(testing_ours_null_a1_pos, results[[2]]$results_our)
  colnames(testing_ours_null_a1_pos)[i+2] <- paste0('test', i)
}

a_1      <- -0.25
sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))

testing_SiZer_null_a1_neg <- expand.grid(u = different_i, h = different_h)
testing_ours_null_a1_neg  <- expand.grid(u = different_i, h = different_h)

for (i in 1:num_of_reps){
  results <- comparing_us_and_Sizer(different_i, different_h, alpha, T_size, a_1, sigma_eta, sigmahat, trend_function, method)
  
  testing_SiZer_null_a1_neg <- cbind(testing_SiZer_null_a1_neg, results[[1]]$results_their)
  colnames(testing_SiZer_null_a1_neg)[i+2] <- paste0('test', i)
  
  testing_ours_null_a1_neg  <- cbind(testing_ours_null_a1_neg, results[[2]]$results_our)
  colnames(testing_ours_null_a1_neg)[i+2] <- paste0('test', i)
}

#CALCULATING POWER OF THE TESTS

a_1      <- 0.25
sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))
trend_function  <- numeric(T_size)
for (i in 1:T_size) {trend_function[i] = (i - 0.5*T_size) * 2.25/T_size}

testing_SiZer_altern_a1_pos <- expand.grid(u = different_i, h = different_h)
testing_ours_altern_a1_pos  <- expand.grid(u = different_i, h = different_h)

for (i in 1:num_of_reps){
  results <- comparing_us_and_Sizer(different_i, different_h, alpha, T_size, a_1, sigma_eta, sigmahat, trend_function, method)
  
  testing_SiZer_altern_a1_pos <- cbind(testing_SiZer_altern_a1_pos, results[[1]]$results_their)
  colnames(testing_SiZer_altern_a1_pos)[i+2] <- paste0('test', i)
  
  testing_ours_altern_a1_pos  <- cbind(testing_ours_altern_a1_pos, results[[2]]$results_our)
  colnames(testing_ours_altern_a1_pos)[i+2] <- paste0('test', i)
}

a_1      <- -0.25
sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))
trend_function  <- numeric(T_size)
for (i in 1:T_size) {trend_function[i] = (i - 0.5*T_size) * 1.25/T_size}

testing_SiZer_altern_a1_neg <- expand.grid(u = different_i, h = different_h)
testing_ours_altern_a1_neg  <- expand.grid(u = different_i, h = different_h)

for (i in 1:num_of_reps){
  results <- comparing_us_and_Sizer(different_i, different_h, alpha, T_size, a_1, sigma_eta, sigmahat, trend_function, method)
  
  testing_SiZer_altern_a1_neg <- cbind(testing_SiZer_altern_a1_neg, results[[1]]$results_their)
  colnames(testing_SiZer_altern_a1_neg)[i+2] <- paste0('test', i)
  
  testing_ours_altern_a1_neg  <- cbind(testing_ours_altern_a1_neg, results[[2]]$results_our)
  colnames(testing_ours_altern_a1_neg)[i+2] <- paste0('test', i)
}