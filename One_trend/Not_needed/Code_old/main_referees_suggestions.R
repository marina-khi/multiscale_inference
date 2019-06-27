library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/grid_construction.r")
dyn.load("functions/C_code/kernel_weights.dll")
source("functions/C_code/kernel_weights.r")
source("functions/multiscale_statistics.r")
source("functions/critical_value.r")
source("functions/multiscale_test.r")

source("functions/functions_new.R")
source("functions/simulating_data.r")

#############################
#Point 5 in Referee Report 1#
#############################

#PLOTTING SIZER MAPS FOR COMPARISON

#Defining necessary parameters
different_T      <- c(250) #Different lengths of time series for which we calculate size and power
alpha            <- 0.05 #Different alpha for which we calculate size and power
different_a1     <- c(-0.25) #Different AR(1) parameters
different_slopes <- c(1.5) #Slopes for power simulations
sigma_eta        <- 0.3 #Sqrt(variance) for the innovation \eta_t
kappa            <- 0.1                   # parameter to determine order statistic for the version
colorlist        <- c('red', 'purple', 'blue', 'grey')
sim.design       <- "bump"
SimRuns          <- 1000

for (T_size in different_T){
  grid  <- grid_construction(T_size)  # grid of location-bandwidth points 
  wghts <- kernel_weights(T=T_size,grid=grid)
  
  gset <- grid$gset
  N    <- dim(gset)[1]
  
  pdffilename <- paste0("Plots/SiZer_comparison_", sim.design, "_T_", T_size, ".pdf")
  
  pdf(pdffilename, width = 14, height = 7)
  
  par(mfrow = c(length(different_a1)*length(different_T),6), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
  par(mar = c(1, 2, 3, 0.5)) #Margins for each plot
  par(oma = c(1.5, 1.5, 3, 0.2)) #Outer margins
  
  for (a_1 in different_a1){
    set.seed(1)
    
    g_set <- grid$gset
    
    gamma = c()
    for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
      gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
    }
    
    #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
    true_var <- gamma[1] / T_size
    for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}
    
    T_star   <- gamma[1]/true_var
    
    SiZer_matrix   <- calculating_SiZer_matrix(g_set, T_size, T_star, alpha, gamma)  
    
    simulated_data <- simulating_data(T_size, a_1, sigma_eta, sim.design)
    trend          <- simulated_data[[1]]
    data           <- simulated_data[[2]]
    sigmahat       <- simulated_data[[3]]

    plot(seq(from = 1/T_size, to = 1, length.out = T_size), data, ylim = c(0.0, 1.0), ylab = 'data')
    lines(seq(from = 1/T_size, to = 1, length.out = T_size), trend, type = 'l')
    
    results       <- multiscale_test(alpha=alpha, data=data, weights=wghts, sigmahat=sigmahat, grid=grid, kappa=kappa, SimRuns=SimRuns)
    results_SiZer <- Sizer_test(gset, SiZer_matrix, data)

    plot.SiZer(results_SiZer, unique(gset$u), unique(gset$h), ylab=expression(log[10](h)),
               colorlist=colorlist, title = paste0("SiZer results, T=", T_size, ", a1 = ", a_1))
    
    plot.SiZer(results[[1]], unique(gset$u), unique(gset$h), ylab=expression(log[10](h)),
               colorlist=colorlist, title = paste0("MS results, T=", T_size, ", a1 = ", a_1))
    plot.SiZer(results[[2]], unique(gset$u), unique(gset$h), ylab=expression(log[10](h)),
               colorlist=colorlist, title = paste0("Uncorrected version, T=", T_size, ", a1 = ", a_1))
    plot.SiZer(results[[3]], unique(gset$u), unique(gset$h), ylab=expression(log[10](h)),
               colorlist=colorlist, title = paste0("Rowwise version, T=", T_size, ", a1 = ", a_1))
    plot.SiZer(results[[4]], unique(gset$u), unique(gset$h), ylab=expression(log[10](h)),
               colorlist=colorlist, title = paste0("Order statistic version, T=", T_size, ", a1 = ", a_1))
  }
}

mtext(paste0(sim.design, ' curve plus AR(1)'), outer = TRUE, cex = 1.0)
dev.off()


#CALCULATING "SIZE" AND "POWER" FOR COMPARISON

#Defining necessary parameters
num_of_reps  <- 1000

slopes_for_negative <- c(0, 1.25) #Slopes for power calculations for negative a_1
slopes_for_positive <- c(0, 2.25) #Slopes for power calculations for positive a_1

percentiles <- c(500, 750, 850, 950)

path <- "Plots/"

for (T_size in different_T){
  grid  <- grid_construction(T_size)  # grid of location-bandwidth points 
  wghts <- kernel_weights(T=T_size,grid=grid)
  
  gset <- grid$gset
  N    <- dim(gset)[1]
  
  for (a_1 in different_a1){
    g_set <- grid$gset
    
    gamma = c()
    for (k in 0:(T_size-1)){                                            #\gamma(k) = \sigma_\eta^2 * a_1^|k| / (1 - a_1^2)
      gamma = c(gamma, autocovariance_function_AR1(k, a_1, sigma_eta))  #Note that gamma[i] := \gamma(i-1)
    }
    
    #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
    true_var <- gamma[1] / T_size
    for (k in 1:(T_size-1)){true_var = true_var + (2/T_size) * (1 - k/T_size) * gamma[k+1]}
    
    T_star   <- gamma[1]/true_var
    
    SiZer_matrix      <- calculating_SiZer_matrix(g_set, T_size, T_star, alpha, gamma)  
    
    if (a_1 > 0){
      slopes <- slopes_for_positive
    } else {
      slopes <- slopes_for_negative
    }
    for (slope in slopes){
      testing_SiZer <- grid$gset
      testing_ms    <- grid$gset
      testing_uncor <- grid$gset
      testing_order <- grid$gset
      testing_rows  <- grid$gset 
      
      set.seed(1)
      for (i in 1:num_of_reps){
        simulated_data <- simulating_data(T_size, a_1, sigma_eta, sim.design, slope)
        trend          <- simulated_data[[1]]
        data           <- simulated_data[[2]]
        sigmahat       <- simulated_data[[3]]
      
        results       <- multiscale_test(alpha=alpha, data=data, weights=wghts, sigmahat=sigmahat, grid=grid, kappa=kappa, SimRuns=SimRuns)
        results_SiZer <- Sizer_test(gset, SiZer_matrix, data)
        
        testing_SiZer <-  cbind(testing_SiZer, results_SiZer)
        colnames(testing_SiZer)[i+2] <- paste0('test', i)
        
        testing_ms  <- cbind(testing_ms, results[[1]])
        colnames(testing_ms)[i+2] <- paste0('test', i)
        
        testing_uncor  <- cbind(testing_uncor, results[[2]])
        colnames(testing_uncor)[i+2] <- paste0('test', i)

        testing_order  <- cbind(testing_order, results[[3]])
        colnames(testing_order)[i+2] <- paste0('test', i)
        
        testing_rows  <- cbind(testing_rows, results[[4]])
        colnames(testing_rows)[i+2] <- paste0('test', i)
      }
      
      all_results <- list(testing_SiZer, testing_ms, testing_uncor, testing_order, testing_rows)
      rm(testing_SiZer, testing_ms, testing_uncor, testing_order, testing_rows)
      saveRDS(all_results, file = paste0(path, "testing_results_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
    }
  }
}

#Based on those simulations, we are interpreting the results.
#Parallel coordinate plots and represnetative SiZer maps under the null
for (T_size in different_T){
  different_i  <- seq(from = 5/T_size, to = 1, by = 5/T_size)
  different_h  <- seq(from = 3/T_size, to = 1/4+3/T_size, by = 5/T_size)

  grid_for_plotting <- expand.grid(u = different_i, h = different_h)

  for (a_1 in different_a1){
    if (a_1 > 0){
      slopes <- slopes_for_positive
    } else {
      slopes <- slopes_for_negative
    }
    for (slope in slopes){
      
      testing_results <- readRDS(file = paste0(path, "testing_results_T_", T_size, "_slope_", slope*100, "_a1_", a_1*100, ".RData"))
      
      testing_results_SiZer<-testing_results[[1]]
      testing_results_SiZer[testing_results_SiZer == 2] <- NA
      
      results_SiZer  <- aggregate(. ~ h, testing_results_SiZer, FUN = function(x) sum(abs(x)) > 0)
      results_ms     <- aggregate(. ~ h, testing_results[[2]], FUN = function(x) sum(abs(x)) > 0)
      results_uncor  <- aggregate(. ~ h, testing_results[[3]], FUN = function(x) sum(abs(x)) > 0)
      results_order  <- aggregate(. ~ h, testing_results[[4]], FUN = function(x) sum(abs(x)) > 0)
      results_rows   <- aggregate(. ~ h, testing_results[[5]], FUN = function(x) sum(abs(x)) > 0)
      
      results_SiZer$u <- NULL
      results_ms$u    <- NULL
      results_uncor$u <- NULL
      results_order$u <- NULL
      results_rows$u  <- NULL
      
      pdffilename <- paste0(path, "sig_comparison_T_", T_size, "_a1_", a_1*100, "_slope_", slope*100, ".pdf")
      pdf(pdffilename)

      if (slope == 0) {ymax <- 20} else {ymax <- 100}

      plot(results_ms$h, rowMeans(results_ms[, -1])*100, ylim = c(0, ymax), xlab = 'bandwidth',
        ylab = "Rowwise % sig.", main = paste0("Percentage of blue/red pixels for T = ", T_size, ", a_1 = ", a_1, ", slope = ", slope),
        type = 'l')
      points(results_uncor$h, rowMeans(results_uncor[, -1])*100, type = 'l', lty = 2)
      points(results_order$h, rowMeans(results_order[, -1])*100, type = 'l', lty = 3)
      points(results_rows$h, rowMeans(results_rows[, -1])*100, type = 'l', lty = 4)
      points(results_SiZer$h, rowMeans(results_SiZer[, -1])*100, type = 'l', lty = 5)
      
      legend(0.03, ymax, legend=c("MS method", "Uncorrected version", "Order version", "Rowwise version", "SiZer method"), lty=1:5, cex=0.8)
      dev.off()

      if (slope == 0){
        ordered_SiZer <- order(colSums(abs(testing_results[[1]][-c(1, 2)])))
        ordered_ms    <- order(colSums(abs(testing_results[[2]][-c(1, 2)])))
        ordered_uncor <- order(colSums(abs(testing_results[[3]][-c(1, 2)])))
        ordered_order <- order(colSums(abs(testing_results[[4]][-c(1, 2)])))
        ordered_rows  <- order(colSums(abs(testing_results[[5]][-c(1, 2)])))
        
        pdffilename <- paste0(path, "representatives_T_", T_size, "_a1_", a_1*100, "_slope_0.pdf")
        pdf(pdffilename)

        par(mfrow = c(length(percentiles),5), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
        par(mar = c(1, 2, 3, 0.5)) #Margins for each plot
        par(oma = c(1.5, 1.5, 3, 0.2)) #Outer margins

        for (percentile in percentiles) {
          plot.SiZer.representatives(different_i, different_h, testing_results[[1]], ordered_SiZer[percentile]+2, colorlist, title ="SiZer, ")
          plot.SiZer.representatives(different_i, different_h, testing_results[[2]], ordered_ms[percentile]+2, colorlist, title ="MS, ")
          plot.SiZer.representatives(different_i, different_h, testing_results[[3]], ordered_uncor[percentile]+2, colorlist, title ="Uncorrected, ")
          plot.SiZer.representatives(different_i, different_h, testing_results[[4]], ordered_order[percentile]+2, colorlist, title ="Ordered version, ")
          plot.SiZer.representatives(different_i, different_h, testing_results[[5]], ordered_rows[percentile]+2, colorlist, title ="Rowwise, ")
        }
        mtext(paste0('Under the null plus AR(1), T = ', T_size, ', a1 = ', a_1), outer = TRUE, cex = 1.0)
        dev.off()
      }

      rm(testing_results)
    }
  }
}



#
#THIS IS NOT YET DEBUGGED!
#
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
# 
# 
# 
# #############################
# #Point 7 in Referee Report 1#
# #############################
# 
# 
# source("Shape/functions_for_referees.R")
# 
# #Defining necessary parameters
# alpha <- 0.05 #alpha for calculating quantiles
# 
# test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported.
# pdffilename = paste0("JRSSB_submission/Plots/temperature_data_SiZer_maps.pdf") #Filename for the graph
# colorlist   <- c('red', 'purple', 'blue', 'grey')
# 
# #Recoding testing problem and type of kernel estimator 
# if (test_problem == "zero"){
#   kernel_ind = 1
# } else if (test_problem == "constant"){
#   kernel_ind = 2
# } else {
#   print('Given testing problem is currently not supported')
# }
# 
# #Loading the real data for yearly temperature in England
# temperature  <- read.table("Shape/data/cetml1659on.dat", header = TRUE, skip = 6)
# yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
# T_tempr      <- length(yearly_tempr)
# 
# different_i <- seq(from = 5/T_tempr, to = 1, by = 5/T_tempr)
# different_h <- seq(from = 3/T_tempr, to = 1/4+3/T_tempr, by = 5/T_tempr)
# 
# #Setting tuning parameters for testing
# p <- 2
# q <- 25
# r <- 10
# 
# 
# #Data analysis
# parameters <- estimating_variance_new(yearly_tempr, q, order = p, r) 
# 
# sigma_hat <- parameters[[1]]
# a_hat     <- parameters[[2]]
# sigma_eta <- parameters[[3]]
# 
# gamma = c()
# for (k in 0:(T_tempr-1)){
#   gamma_temp <- gamma
#   gamma = c(gamma, autocovariance_function_AR2(k, a_hat[[1]], a_hat[[2]], sigma_eta, gamma_temp))  #Note that gamma[i] := \gamma(i-1)
# }
# rm(gamma_temp)
# 
# 
# #Calculating \Var(\bar{Y}) based on the true values of gamma(k)
# true_var <- gamma[1] / T_tempr
# for (k in 1:(T_tempr-1)){true_var = true_var + (2/T_tempr) * (1 - k/T_tempr) * gamma[k+1]}
# 
# T_star   <- gamma[1]/true_var
# 
# SiZer_matrix      <- calculating_SiZer_matrix(different_i, different_h, T_tempr, T_star, alpha, gamma)  
# 
# g_t_set <- psihat_statistic(yearly_tempr, SiZer_matrix, kernel_ind, sigma_hat)[[1]]
# gaussian_quantile <- calculating_gaussian_quantile(T_tempr, SiZer_matrix, "data", kernel_ind, alpha)
# g_t_set$gaussian_quantile <- gaussian_quantile
# 
# g_t_set_temp <- NULL
# 
# for (bandwidth in different_h){
#   SiZer_matrix_temp                           <- subset(SiZer_matrix, h == bandwidth, select = c(u, h, lambda))
#   if (nrow(SiZer_matrix_temp)>0){
#     gaussian_quantile_rowwise                   <- calculating_gaussian_quantile(T_tempr, SiZer_matrix_temp, paste0("data_h_", bandwidth*1000), kernel_ind, alpha)
#     g_t_set_temp_temp                           <- psihat_statistic(yearly_tempr, SiZer_matrix_temp, kernel_ind, sigma_hat)[[1]]
#     g_t_set_temp_temp$gaussian_quantile_rowwise <- gaussian_quantile_rowwise
#     g_t_set_temp                                <- rbind(g_t_set_temp, g_t_set_temp_temp)
#   }
# }
# 
# g_t_set_rowwise <- merge(SiZer_matrix, g_t_set_temp, by = c('h', 'u', 'lambda'))
# 
# 
# for (row in 1:nrow(g_t_set)){
#   i              = g_t_set[row, 'u']
#   h              = g_t_set[row, 'h']
#   q_h            = g_t_set[row, 'q_h']
#   sd_m_hat_prime = g_t_set[row, 'sd']
#   
#   XtWX_inverse_XtW = g_t_set$XtWX_inv_XtW[[row]]
#   
#   if (!is.null(XtWX_inverse_XtW)) {
#     m_hat_prime <- (XtWX_inverse_XtW %*% yearly_tempr)[2]
#     if (m_hat_prime - q_h * sd_m_hat_prime > 0){
#       g_t_set$results_their[[row]] = 1
#     } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
#       g_t_set$results_their[[row]] = -1
#     } else {
#       g_t_set$results_their[[row]] = 0
#     }
#   } else {
#     g_t_set$results_their[[row]] = 2
#   }
#   
#   if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + g_t_set[row, 'gaussian_quantile']){
#     g_t_set$results_our[[row]] = 1
#   } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + g_t_set[row, 'gaussian_quantile']){
#     g_t_set$results_our[[row]] = -1
#   } else {
#     g_t_set$results_our[[row]] = 0
#   }
#   
#   if (g_t_set_rowwise[row, 'values_with_sign'] > g_t_set_rowwise[row, 'lambda'] + g_t_set_rowwise[row, 'gaussian_quantile_rowwise']){
#     g_t_set_rowwise$results_our_rowwise[[row]] = 1
#   } else if (-g_t_set_rowwise[row, 'values_with_sign'] > g_t_set_rowwise[row, 'lambda'] + g_t_set_rowwise[row, 'gaussian_quantile_rowwise']){
#     g_t_set_rowwise$results_our_rowwise[[row]] = -1
#   } else {
#     g_t_set_rowwise$results_our_rowwise[[row]] = 0
#   }
#   
# }
# result_SiZer <- subset(g_t_set, select = c(u, h, results_their))
# result_our   <- subset(g_t_set, select = c(u, h, results_our))
# result_our_rowwise        <- subset(g_t_set_rowwise, select = c(u, h, results_our_rowwise))
# 
# pdf(pdffilename, width = 7, height = 3)
# 
# par(mfrow = c(1, 3), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
# par(mar = c(1, 2, 3, 0.5)) #Margins for each plot
# par(oma = c(1.5, 1.5, 3, 0.2)) #Outer margins
# 
# plot.SiZer(result_SiZer$results_their, different_i, different_h, ylab=expression(log[10](h)),
#            colorlist=colorlist, title = "SiZer results for the data")
# 
# plot.SiZer(result_our$results_our, different_i, different_h, ylab=expression(log[10](h)),
#            colorlist=colorlist, title = "Our global results for the data")
# 
# plot.SiZer(result_our_rowwise$results_our_rowwise, different_i, different_h, ylab=expression(log[10](h)),
#            colorlist=colorlist, title = "Our rowwise results for the data")
# dev.off()

