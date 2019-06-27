rm(list=ls())

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/grid_construction.r")
dyn.load("functions/kernel_weights.dll")
source("functions/kernel_weights.r")
source("functions/multiscale_statistics.r")
source("functions/multiscale_quantiles.r")
source("functions/multiscale_testing.r")

source("functions/long_run_variance.r")
source("functions/sim.r")


# Parameters
different_T     <- c(250, 500, 1000)
different_a1    <- c(-0.25, -0.5, -0.9, 0.25, 0.5, 0.9)
different_alpha <- c(0.01, 0.05, 0.1)

Nsim       <- 1000                  # number of simulation runs for size/power calculations
sigma_eta  <- 1                     # standard deviation of the innovation term in the AR model
SimRuns    <- 1000                  # number of simulation runs to produce critical values
kappa      <- 0.1                   # parameter to determine order statistic for the version
                                    # of the multiscale statistic from Section ?? 


#Constructing necessary auxiliary parameters
number_of_cols            <- length(different_a1) * (length(different_alpha) + 1)
size_matrix_ms            <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix_ms)  <- different_T

k <- 1
for (a1 in different_a1){
  size_ms <- numeric(length(different_alpha))
  for (T in different_T){
    #set.seed(1)
    # Construct grid with effective sample size ESS.star(u,h) >= 5 for any (u,h)
    grid      <- grid_construction(T)
    gset      <- grid$gset
    u.grid    <- sort(unique(gset[,1]))
    h.grid    <- sort(unique(gset[,2]))

    # Compute kernel weights and critical value for multiscale test  
    wghts  <- kernel_weights(T=T, grid=grid)
    
    filename = paste0("quantiles/distr_T_", T,"_ms.RData")
    if(!file.exists(filename)) {
      quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=kappa, SimRuns=SimRuns)
      save(quants, file = filename)
    } else {
      load(filename)
    }
    
    size_matrix_temp <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))

    for (i in 1:Nsim){
      #Simulating the time series
      data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
      data           <- data.simulated$data
      trend          <- data.simulated$trend
      sigma_true     <- data.simulated$sigma
          
      #Estimating the coefficients for the ts and the long-run variance
      AR.struc  <- AR_lrv(data=data, q=25, r.bar=10, p=1)    
      a.hat     <- AR.struc$ahat 
      vareta    <- AR.struc$vareta   
      sigma_hat <- sqrt(AR.struc$lrv)
      autocov   <- AR_acf(coefs=a.hat, var.eta=vareta, len=T)
      
      stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid) 
      vals     <- stats$values
      
      #Based on the same values of the test statistic, perform the test at dife
      for (j in 1:length(different_alpha)){
        alpha <- different_alpha[j]
        test.res <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
        if (sum(abs(test.res$test_ms)) >0) {size_matrix_temp[i, j] <- 1} else {size_matrix_temp[i, j] <- 0}
      }
    }
    size_ms <- cbind(size_ms, colSums(size_matrix_temp)/Nsim)
    rm(size_matrix_temp)
  }
  size_matrix_ms[, (k * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(k * (length(different_alpha) + 1))] <- t(size_ms[, -1])
  k <- k + 1
}

for (i in 1:2){
  print.xtable(xtable(size_matrix_ms[, (number_of_cols/2 * (i-1) + 1):(number_of_cols/2 * i)], digits = c(3), align = paste(replicate(number_of_cols/2 + 1, "c"), collapse = "")),
               type="latex", file=paste0("plots/size_part", i, ".tex"), include.colnames = FALSE)
  
}


# Parameters
different_T     <- seq(1000, 10000, by = 1000)
different_a1    <- c(-0.9, 0.9)


#Constructing necessary auxiliary parameters
number_of_cols                       <- length(different_a1) * (length(different_alpha) + 1)
size_matrix_ms_big_samples           <- matrix(NA, nrow = length(different_T), ncol = number_of_cols)
rownames(size_matrix_ms_big_samples) <- different_T

k <- 1
for (a1 in different_a1){
  size_ms <- numeric(length(different_alpha))
  for (T in different_T){
    #set.seed(1)
    # Construct grid with effective sample size ESS.star(u,h) >= 5 for any (u,h)
    grid      <- grid_construction(T)
    gset      <- grid$gset
    u.grid    <- sort(unique(gset[,1]))
    h.grid    <- sort(unique(gset[,2]))
    
    # Compute kernel weights and critical value for multiscale test  
    wghts  <- kernel_weights(T=T, grid=grid)
    
    filename = paste0("quantiles/distr_T_", T,"_ms.RData")
    if(!file.exists(filename)) {
      quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=kappa, SimRuns=SimRuns)
      save(quants, file = filename)
    } else {
      load(filename)
    }
    
    size_matrix_temp <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
    
    for (i in 1:Nsim){
      #Simulating the time series
      data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
      data           <- data.simulated$data
      trend          <- data.simulated$trend
      sigma_true     <- data.simulated$sigma
      
      #Estimating the coefficients for the ts and the long-run variance
      AR.struc  <- AR_lrv(data=data, q=25, r.bar=10, p=1)    
      a.hat     <- AR.struc$ahat 
      vareta    <- AR.struc$vareta   
      sigma_hat <- sqrt(AR.struc$lrv)
      autocov   <- AR_acf(coefs=a.hat, var.eta=vareta, len=T)
      
      stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid) 
      vals     <- stats$values
      
      #Based on the same values of the test statistic, perform the test at dife
      for (j in 1:length(different_alpha)){
        alpha <- different_alpha[j]
        test.res <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
        if (sum(abs(test.res$test_ms)) >0) {size_matrix_temp[i, j] <- 1} else {size_matrix_temp[i, j] <- 0}
      }
    }
    size_ms <- cbind(size_ms, colSums(size_matrix_temp)/Nsim)
    rm(size_matrix_temp)
  }
  size_matrix_ms_big_samples[, (k * (length(different_alpha) + 1) - (length(different_alpha) - 1)):(k * (length(different_alpha) + 1))] <- t(size_ms[, -1])
  k <- k + 1
}

print.xtable(xtable(size_matrix_ms_big_samples, digits = c(3), align = paste(replicate(number_of_cols + 1, "c"), collapse = "")),
               type="latex", file="plots/size_big_samples.tex", include.colnames = FALSE)

