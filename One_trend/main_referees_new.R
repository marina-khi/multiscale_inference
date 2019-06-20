library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/estimating_sigma_new.R")

source("Shape/functions_merging/grid_construction.r")
dyn.load("Shape/functions_merging/kernel_weights.dll")
source("Shape/functions_merging/kernel_weights.r")
source("Shape/functions_merging/multiscale_statistics.r")
source("Shape/functions_merging/critical_value.r")
source("Shape/functions_merging/multiscale_test.r")

source("Shape/functions_new.R")
source("Shape/simulating_data.r")

#############################
#Point 5 in Referee Report 1#
#############################

#PLOTTING SIZER MAPS FOR COMPARISON

#Defining necessary parameters
different_T      <- c(250) #Different lengths of time series for which we calculate size and power
alpha            <- 0.05 #Different alpha for which we calculate size and power
different_a1     <- c(-0.25, 0.25) #Different AR(1) parameters
different_slopes <- c(1.5) #Slopes for power simulations
sigma_eta        <- 1 #Sqrt(variance) for the innovation \eta_t
kappa            <- 0.1                   # parameter to determine order statistic for the version
colorlist        <- c('red', 'purple', 'blue', 'grey')
sim.design       <- "spike"


for (T_size in different_T){
  grid  <- grid_construction(T_size)  # grid of location-bandwidth points 
  wghts <- kernel_weights(T=T_size,grid=grid)
  
  gset <- grid$gset
  N    <- dim(gset)[1]
  
  pdffilename <- paste0("JRSSB_submission/Plots/SiZer_comparison_", sim.design, "_T_", T_size, "_new.pdf")
  
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
    
    SiZer_matrix      <- calculating_SiZer_matrix(g_set, T_size, T_star, alpha, gamma)  
    
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

path <- "JRSSB_submission/Plots/"

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

      results_SiZer  <- aggregate(. ~ h, testing_results[[1]], FUN = function(x) sum(abs(x)) > 0)
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
