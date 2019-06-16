source("Shape/functions_for_referees.R")
source("Shape/estimating_sigma_new.R")
dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_without_lambda.dll")
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

#PLOTTING SIZER MAPS FOR COMPARISON

#Defining necessary parameters
sigma_eta           <- 0.1     #Standard deviation of the innovation term, for blocks = 0.1, for other signals = 1
alpha               <- 0.05    #Level of significance
different_T         <- c(1000) #Different lengths of time series for which we compare SiZer and our method
different_a1        <- c(-0.25, 0, 0.25) #Different a_1 in AR(1) model

kernel_ind <- 2

#method = 'rowwise'
method = 'global'


#THIS IS FOR PLOTTING SIZER MAPS
h_for_blocks <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)
t_for_blocks <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.4, 0.65, 0.76, 0.78, 0.81)

colorlist  <- c('red', 'purple', 'blue', 'grey')


for (T_size in different_T){
  different_i <- seq(from = 5/T_size, to = 1, by = 5/T_size)
  different_h <- seq(from = 3/T_size, to = 1/5+3/T_size, by = 5/T_size)

  trend_function  <- numeric(T_size)
  
  #for (i in 1:T_size) {trend_function[i] = (0.6 / 9.2) * (sum(h_for_blocks*(1 + sign(i/T_size - t_for_blocks))/2) + 2) + 0.2}
  #for (i in 1:T_size) {trend_function[i] = 0}
  for (i in 1:T_size) {trend_function[i] = sinpi(6*i/T_size)}
  pdffilename <- paste0("JRSSB_submission/Plots/SiZer_comparison_sin_", method, "_T_", T_size, ".pdf")
  
  pdf(pdffilename)
  
  par(mfrow = c(length(different_a1)*length(different_T),3), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(1, 2, 3, 0.5)) #Margins for each plot
  par(oma = c(1.5, 1.5, 3, 0.2)) #Outer margins
  
  for (a_1 in different_a1){
    sigmahat   <- sqrt(sigma_eta^2/((1 - a_1)^2))

    results <- comparing_us_and_Sizer(different_i, different_h, alpha, T_size, a_1, sigma_eta, sigmahat, trend_function, method)
    plot.SiZer(results[[1]]$results_their, different_i, different_h, ylab=expression(log[10](h)),
                    colorlist=colorlist, title = paste0("SiZer results, T=", T_size, ", a1 = ", a_1))

    plot.SiZer(results[[2]]$results_our, different_i, different_h, ylab=expression(log[10](h)),
                    colorlist=colorlist, title = paste0("Our results, T=", T_size, ", a1 = ", a_1))
  }
}

mtext(paste0('Sine curve plus AR(1), ', method, ' method'), outer = TRUE, cex = 1.0)
dev.off()

# 
# #CALCULATING "SIZE" AND "POWER" FOR COMPARISON
# 
# #Defining necessary parameters
# sigma_eta           <- 1    #Standard deviation of the innovation term
# alpha               <- 0.05  #Level of significance
# 
# kernel_ind <- 2
# 
# num_of_reps  <- 1000
# different_T  <- c(250)
# different_a1 <- c(-0.5, 0, 0.25) #Different a_1 in AR(1) model
# 
# slopes_for_negative <- c(0, 1.25) #Slopes for power calculations for negative a_1
# slopes_for_positive <- c(0, 2.25) #Slopes for power calculations for positive a_1
# 
# colorlist   <- c('red', 'purple', 'blue', 'grey')
# percentiles <- c(500, 750, 850, 950)
# 
# path <- "JRSSB_submission/Plots/"
# 
# 
# #1000 simulations for SiZer, our global method and our rowwise method simultaneously
# for (T_size in different_T){
#   different_i  <- seq(from = 5/T_size, to = 1, by = 5/T_size)
#   different_h  <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
# 
#   for (a_1 in different_a1){
#     sigmahat <- sqrt(sigma_eta^2/((1 - a_1)^2))
# 
#     if (a_1 > 0){
#       slopes <- slopes_for_positive
#     } else {
#       slopes <- slopes_for_negative
#     }
#     for (slope in slopes){
#       set.seed(1)
# 
#       trend_function  <- numeric(T_size)
#       for (i in 1:T_size) {trend_function[i] = (i - 0.5*T_size) * slope/T_size}
# 
#       testing_SiZer        <- expand.grid(u = different_i, h = different_h)
#       testing_ours_global  <- expand.grid(u = different_i, h = different_h)
#       testing_ours_rowwise <- expand.grid(u = different_i, h = different_h)
#       testing_ours_without_lambda <- expand.grid(u = different_i, h = different_h)
#       
#       for (i in 1:num_of_reps){
#         results <- comparing_us_and_Sizer_global_and_rowwise(different_i, different_h, alpha, T_size, a_1, sigma_eta, sigmahat, trend_function)
# 
#         testing_SiZer <-  merge(testing_SiZer, results[[1]], by = c('h', 'u'))
#         colnames(testing_SiZer)[i+2] <- paste0('test', i)
# 
#         testing_ours_global  <- merge(testing_ours_global, results[[2]], by = c('h', 'u'))
#         colnames(testing_ours_global)[i+2] <- paste0('test', i)
# 
#         testing_ours_rowwise  <- merge(testing_ours_rowwise, results[[3]], by = c('h', 'u'))
#         colnames(testing_ours_rowwise)[i+2] <- paste0('test', i)
#         
#         testing_ours_without_lambda <- merge(testing_ours_without_lambda, results[[4]], by = c('h', 'u'))
#         colnames(testing_ours_without_lambda)[i+2] <- paste0('test', i)
#       }
# 
#       save(testing_SiZer, file = paste0(path, "testing_SiZer_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
#       save(testing_ours_global, file = paste0(path, "testing_ours_global_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
#       save(testing_ours_rowwise, file = paste0(path, "testing_ours_rowwise_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
#       save(testing_ours_without_lambda, file = paste0(path, "testing_ours_without_lambda_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
# 
#       rm(testing_SiZer, testing_ours_global, testing_ours_rowwise, testing_ours_without_lambda)
#     }
#   }
# }
# 
# #Based on those simulations, we are interpreting the results.
# #Parallel coordinate plots and represnetative SiZer maps under the null
# for (T_size in different_T){
#   different_i  <- seq(from = 5/T_size, to = 1, by = 5/T_size)
#   different_h  <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)
# 
#   grid_for_plotting <- expand.grid(u = different_i, h = different_h)
# 
#   for (a_1 in different_a1){
#     if (a_1 > 0){
#       slopes <- slopes_for_positive
#     } else {
#       slopes <- slopes_for_negative
#     }
#     for (slope in slopes){
# 
#       load(file = paste0(path, "testing_SiZer_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
#       load(file = paste0(path, "testing_ours_global_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
#       load(file = paste0(path, "testing_ours_rowwise_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
#       load(file = paste0(path, "testing_ours_without_lambda_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
#       
#       results_SiZer                 <- aggregate(. ~ h, testing_SiZer, FUN = function(x) sum(abs(x)) > 0)
#       results_ours_global           <- aggregate(. ~ h, testing_ours_global, FUN = function(x) sum(abs(x)) > 0)
#       results_ours_rowwise          <- aggregate(. ~ h, testing_ours_rowwise, FUN = function(x) sum(abs(x)) > 0)
#       results_ours_without_lambda   <- aggregate(. ~ h, testing_ours_without_lambda, FUN = function(x) sum(abs(x)) > 0)
#       
#       results_SiZer$u               <- NULL
#       results_ours_global$u         <- NULL
#       results_ours_rowwise$u        <- NULL
#       results_ours_without_lambda$u <- NULL
#       
# 
#       pdffilename <- paste0(path, "rowwise_sig_comparison_T_", T_size, "_a1_", a_1*100, "_slope_", slope*100, ".pdf")
#       pdf(pdffilename)
# 
#       if (slope == 0) {ymax <- 15} else {ymax <- 100}
#       
#       plot(results_ours_global$h, rowMeans(results_ours_global[, -1])*100, ylim = c(0, ymax), xlab = 'bandwidth',
#         ylab = "Rowwise % sig.", main = paste0("Percentage of blue/red pixels for T = ", T_size, ", a_1 = ", a_1, ", slope = ", slope),
#         type = 'l')
#       points(results_ours_rowwise$h, rowMeans(results_ours_rowwise[, -1])*100, type = 'l', lty = 2)
#       points(results_ours_without_lambda$h, rowMeans(results_ours_without_lambda[, -1])*100, type = 'l', lty = 3)
#       points(results_SiZer$h, rowMeans(results_SiZer[, -1])*100, type = 'l', lty = 4)
#       legend(0.03, ymax, legend=c("Our global method", "Our rowwise method", "Our method without lambda", "SiZer"), lty=1:4, cex=0.8)
#       dev.off()
#       
#       if (slope == 0){
# 
#         ordered_SiZer               <- order(colSums(abs(testing_SiZer[-c(1, 2)])))
#         ordered_ours_global         <- order(colSums(abs(testing_ours_global[-c(1, 2)])))
#         ordered_ours_rowwise        <- order(colSums(abs(testing_ours_rowwise[-c(1, 2)])))
#         ordered_ours_without_lambda <- order(colSums(abs(testing_ours_without_lambda[-c(1, 2)])))
#         
#         pdffilename <- paste0(path, "representatives_T_", T_size, "_a1_", a_1*100, "_slope_0.pdf")
#         pdf(pdffilename)
#         
#         par(mfrow = c(length(percentiles),4), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
#         par(mar = c(1, 2, 3, 0.5)) #Margins for each plot
#         par(oma = c(1.5, 1.5, 3, 0.2)) #Outer margins
#         
#         for (percentile in percentiles) {
#           plot.SiZer.representatives(different_i, different_h, testing_SiZer, ordered_SiZer[percentile]+2, colorlist, title ="SiZer, ")
#           plot.SiZer.representatives(different_i, different_h, testing_ours_global, ordered_ours_global[percentile]+2, colorlist, title ="Global, ")
#           plot.SiZer.representatives(different_i, different_h, testing_ours_rowwise, ordered_ours_rowwise[percentile]+2, colorlist, title ="Rowwise, ")
#           plot.SiZer.representatives(different_i, different_h, testing_ours_without_lambda, ordered_ours_without_lambda[percentile]+2, colorlist, title ="Without lambda, ")
#         }
#         mtext(paste0('Under the null plus AR(1), T = ', T_size, ', a1 = ', a_1), outer = TRUE, cex = 1.0)
#         dev.off()
#       }
#       
#       rm(testing_SiZer, testing_ours_global, testing_ours_rowwise)
#     }
#   }
# }

#############################
#Point 7 in Referee Report 1#
#############################

#Defining necessary parameters
alpha <- 0.05 #alpha for calculating quantiles

test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported.
pdffilename = paste0("JRSSB_submission/Plots/temperature_data.pdf") #Filename for the graph
colorlist   <- c('red', 'purple', 'blue', 'grey')

#Recoding testing problem and type of kernel estimator 
if (test_problem == "zero"){
  kernel_ind = 1
} else if (test_problem == "constant"){
  kernel_ind = 2
} else {
  print('Given testing problem is currently not supported')
}

#Loading the real data for yearly temperature in England
temperature  <- read.table("Shape/data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
T_tempr      <- length(yearly_tempr)

different_i <- seq(from = 5/T_tempr, to = 1, by = 5/T_tempr)
different_h <- seq(from = 3/T_tempr, to = 1/4+3/T_tempr, by = 5/T_tempr)

#Setting tuning parameters for testing
p <- 2
q <- 25
r <- 10


#Data analysis
parameters <- estimating_variance_new(yearly_tempr, q, order = p, r) 

sigma_hat <- parameters[[1]]
a_hat     <- parameters[[2]]
sigma_eta <- parameters[[3]]

gamma = c()
for (k in 0:(T_tempr-1)){
  gamma_temp <- gamma
  gamma = c(gamma, autocovariance_function_AR2(k, a_hat[[1]], a_hat[[2]], sigma_eta, gamma_temp))  #Note that gamma[i] := \gamma(i-1)
}

#cat("Autocovariance function:", gamma, "\n")

#Calculating \Var(\bar{Y}) based on the true values of gamma(k)
true_var <- gamma[1] / T_tempr
for (k in 1:(T_tempr-1)){true_var = true_var + (2/T_tempr) * (1 - k/T_tempr) * gamma[k+1]}

T_star   <- gamma[1]/true_var

SiZer_matrix      <- calculating_SiZer_matrix(different_i, different_h, T_tempr, T_star, alpha, gamma, a_hat, sigma_eta)  

g_t_set <- psihat_statistic(yearly_tempr, SiZer_matrix, kernel_ind, sigma_hat)[[1]]
gaussian_quantile <- calculating_gaussian_quantile(T_tempr, SiZer_matrix, "data", kernel_ind, alpha)
g_t_set$gaussian_quantile <- gaussian_quantile

g_t_set_temp <- NULL

for (bandwidth in different_h){
  SiZer_matrix_temp                           <- subset(SiZer_matrix, h == bandwidth, select = c(u, h, lambda))
  if (nrow(SiZer_matrix_temp)>0){
    gaussian_quantile_rowwise                   <- calculating_gaussian_quantile(T_tempr, SiZer_matrix_temp, paste0("data_h_", bandwidth*1000, "_a1_", a_hat*100), kernel_ind, alpha)
    g_t_set_temp_temp                           <- psihat_statistic(yearly_tempr, SiZer_matrix_temp, kernel_ind, sigma_hat)[[1]]
    g_t_set_temp_temp$gaussian_quantile_rowwise <- gaussian_quantile_rowwise
    g_t_set_temp                                <- rbind(g_t_set_temp, g_t_set_temp_temp)
  }
}

g_t_set_rowwise <- merge(SiZer_matrix, g_t_set_temp, by = c('h', 'u', 'lambda'))


for (row in 1:nrow(g_t_set)){
  i              = g_t_set[row, 'u']
  h              = g_t_set[row, 'h']
  q_h            = g_t_set[row, 'q_h']
  sd_m_hat_prime = g_t_set[row, 'sd']
  
  XtWX_inverse_XtW = g_t_set$XtWX_inv_XtW[[row]]
  
  if (!is.null(XtWX_inverse_XtW)) {
    m_hat_prime <- (XtWX_inverse_XtW %*% yearly_tempr)[2]
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
  
  if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + g_t_set[row, 'gaussian_quantile']){
    g_t_set$results_our[[row]] = 1
  } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + g_t_set[row, 'gaussian_quantile']){
    g_t_set$results_our[[row]] = -1
  } else {
    g_t_set$results_our[[row]] = 0
  }
  
  if (g_t_set_rowwise[row, 'values_with_sign'] > g_t_set_rowwise[row, 'lambda'] + g_t_set_rowwise[row, 'gaussian_quantile_rowwise']){
    g_t_set_rowwise$results_our_rowwise[[row]] = 1
  } else if (-g_t_set_rowwise[row, 'values_with_sign'] > g_t_set_rowwise[row, 'lambda'] + g_t_set_rowwise[row, 'gaussian_quantile_rowwise']){
    g_t_set_rowwise$results_our_rowwise[[row]] = -1
  } else {
    g_t_set_rowwise$results_our_rowwise[[row]] = 0
  }
  
}
result_SiZer <- subset(g_t_set, select = c(u, h, results_their))
result_our   <- subset(g_t_set, select = c(u, h, results_our))
result_our_rowwise        <- subset(g_t_set_rowwise, select = c(u, h, results_our_rowwise))

plot.SiZer(result_SiZer$results_their, different_i, different_h, ylab=expression(log[10](h)),
           colorlist=colorlist, title = paste0("SiZer results, T=", T_tempr, ", a1 = ", a_hat))

plot.SiZer(result_our$results_our, different_i, different_h, ylab=expression(log[10](h)),
           colorlist=colorlist, title = paste0("Our global results, T=", T_tempr, ", a1 = ", a_hat))

plot.SiZer(result_our_rowwise$results_our_rowwise, different_i, different_h, ylab=expression(log[10](h)),
           colorlist=colorlist, title = paste0("Our rowwise results, T=", T_tempr, ", a1 = ", a_hat))

