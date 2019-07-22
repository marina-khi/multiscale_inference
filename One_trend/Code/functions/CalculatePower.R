CalculatePower <- function(T, a1, sigma_eta, alpha = 0.05, Nsim = 1000, SimRuns =1000,
                              sigma.type = 'true', q_ = 25,
                              remove.small.ess = FALSE, sim.design = 'bump', bump.height = 1,
                              region = 'increase',
                              power.type = ''){
  # Computes the statistical power of the multiscale test (T_ms), the uncorrected version of the multiscale test (T_uc),
  #   the rowwise version of the multiscale test (T_rw) compared to the power of SiZer test (T_SiZer). The 
  #   The power is calculated only for the case where the errors follow AR(1) process.
  #
  # Args:
  #   T:  Sample size of the time series simulated.
  #   a1: AR(1) coefficient for simulating the error distribution.
  #   sigma_eta: standard deviation of the innovation term in the AR(1) process.
  #   alpha : Confidence level. Default is 5%.
  #   Nsim: number of simulations for power calculations. Default is 1000.
  #   SimRuns: number of simulations for calculating gaussian quantile for T_ms. Default is 1000.
  #   sigma.type: If 'estimated', then the long-run variance used in the test statistic is first estimated by
  #     AR_lrv() function. Otherwise, if 'true' the true theoretical value of the long-run variance is used.
  #   q_: tuning parameter for estimating the long-run variance from Section 4. Default is 25. 
  #   remove.small.ess: If TRUE, then we restrict attention only to those points (u, h) for which 
  #     the effective sample size for correlated data ESS^*(u, h) is not smaller than 5.
  #   sim.design: Type of the trend function under alternative. Can be "constant", "linear", "brokenline", "bump", "blocks" or "sine".
  #     Default is "bump".
  #   bump.height: Height of the bump signal or the slope of "linear" sim.design. Default is 1.
  #   region: If "increase", than the power is calculated only for increases of the trend function. The other option is "decrease".
  #     Default is "increase".
  #   power.type: If 'spurious', then the results are for the sputious power. Otherwise, the results are for "normal" power.
  #     See section 5.1.2 on this.
  #
  # Returns:
  #   power: One vector of length 4 with each entry corresponding to the power of each testing procedure. The order is as follows.
  #     T_ms, T_uc, T_rw, T_SiZer.
  #   power_ms:    A vector of length equal to the number of bandwidths analysed, each entry corresponding to the rowwise power 
  #     of the multiscale testing procedure (T_ms) for this bandwidth.
  #   power_uncor: A vector of length equal to the number of bandwidths analysed, each entry corresponding to the rowwise power 
  #     of the uncorrected version of the multiscale testing procedure (T_uc) for this bandwidth.
  #   power_rows:  A vector of length equal to the number of bandwidths analysed, each entry corresponding to the rowwise power 
  #     of the the rowwise version of the multiscale testing procedure (T_rw) for this bandwidth.
  #   power_SiZer: A vector of length equal to the number of bandwidths analysed, each entry corresponding to the rowwise power 
  #     of the SiZer testing procedure (T_SiZer) for this bandwidth.
  #   h.grid.new:  A vector of the bandwidths analysed.

  #Load necessary functions  
  source("functions/ConstructGrid.r")
  source("functions/multiscale_statistics.r")
  source("functions/multiscale_quantiles.r")
  source("functions/multiscale_testing.r")
  source("functions/long_run_variance.r")
  source("functions/sim.r")
  source("functions/SiZer_functions.r")
  sourceCpp("functions/kernel_weights.cpp")
  sourceCpp("functions/SiZer_functions.cpp")
  
  #Construct grid
  grid      <- grid_construction(T)
  gset      <- grid$gset
  u.grid    <- sort(unique(gset[,1]))
  h.grid    <- sort(unique(gset[,2]))
  autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))  
  
  if (remove.small.ess){
    ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
    deletions <- ess$del
    grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
  }
    
  gset        <- grid$gset
  N           <- as.integer(dim(gset)[1])
  gset.full   <- grid$gset_full
  h.grid.full <- sort(unique(gset.full[,2]))
  h.len       <- length(h.grid.full)
  u.len       <- length(unique(gset.full[,1]))
  pos.full    <- grid$pos_full
  h.grid.new  <- sort(unique(grid$gset[,2]))
  
  u.left  <- (gset[,1] - gset[,2])
  u.right <- (gset[,1] + gset[,2])
    
  #Simulate data one time in order to get the int.plus and int.minus regions
  data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = sim.design)
  
  if (region == 'decrease'){
    int <- data.simulated$int.minus
    coded.value <- -1
  } else if (region == 'increase'){
    int <- data.simulated$int.plus
    coded.value <- 1
  } else {
    cat("Region is not supported \n")
    }
  
  if(is.null(int)){
    pos <- matrix(1, ncol=u.len, nrow=h.len)
  } else { 
    temp                                       <- rep(1,N)
    temp[u.left >= int[2] | u.right <= int[1]] <- 100 
    pos                                        <- rep(NA,length(pos.full)) 
    pos[!is.na(pos.full)]                      <- temp
    pos                                        <- matrix(pos, ncol=u.len, byrow=TRUE)
  }

  if (power.type == 'spurious'){
    coded.value <- 100 * coded.value
  }

  # Compute kernel weights and critical value for multiscale test
  T                      <- as.integer(T) 
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp) 
  storage.mode(gset_cpp) <- "double"
  
  wghts <- matrix(kernel_weights(T, gset_cpp, N), ncol = T, byrow = TRUE)

  filename = paste0("quantiles/distr_T_", T,".RData")
  if((!file.exists(filename)) || (remove.small.ess)) {
    quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
    if (!remove.small.ess){
      save(quants, file = filename)
    }
  } else {
    load(filename)
  }
  
  sizer.wghts  <- SiZer_weights(T=T, grid=grid)
  sizer.std    <- SiZer_std(weights=sizer.wghts, autocov=autocov, T)
  sizer.quants <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
  
  power_matrix_temp <- matrix(NA, nrow = Nsim, ncol = 4)
  
  power_ms_temp     <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
  power_uncor_temp  <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
  power_rows_temp   <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
  power_SiZer_temp  <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
  

  for (i in 1:Nsim){
    #Simulating the time series
    data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = sim.design, slope.fac = bump.height)
    data           <- data.simulated$data
    trend          <- data.simulated$trend
    sigma_true     <- data.simulated$sigma
      
    if (sigma.type == 'estimated'){
      #Estimating the coefficients for the ts and the long-run variance
      AR.struc  <- AR_lrv(data=data, q=q_, r.bar=10, p=1)
      sigma_hat <- sqrt(AR.struc$lrv)
      stats     <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid)
    } else if (sigma.type == 'true'){
      stats     <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_true, grid=grid)
    }
    
    SiZer.values  <- sizer.wghts %*% data
    sizer.vals    <- as.vector(SiZer.values)
    SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)
      
    vals     <- stats$values
    test.res <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
    
    test_ms    <- test.res$test_ms * pos
    test_uncor <- test.res$test_uncor * pos
    test_rows  <- test.res$test_rows * pos
    test_SiZer <- SiZer_results$test * pos
    
    if (sum(test_ms == coded.value, na.rm = TRUE) > 0) {power_matrix_temp[i, 1] <- 1} else {power_matrix_temp[i, 1] <- 0}
    if (sum(test_uncor == coded.value, na.rm = TRUE) > 0) {power_matrix_temp[i, 2] <- 1} else {power_matrix_temp[i, 2] <- 0}
    if (sum(test_rows == coded.value, na.rm = TRUE) > 0) {power_matrix_temp[i, 3] <- 1} else {power_matrix_temp[i, 3] <- 0}
    if (sum(test_SiZer == coded.value, na.rm = TRUE) > 0) {power_matrix_temp[i, 4] <- 1} else {power_matrix_temp[i, 4] <- 0}
    
    for (j in 1:length(h.grid.new)){
      bw  <- h.grid.new[j]
      h_index <- match(bw, h.grid)
      
      if (sum(test_ms[h_index, ] == coded.value, na.rm = TRUE) >0) {power_ms_temp[i, j] <- 1} else {power_ms_temp[i, j] <- 0}
      if (sum(test_uncor[h_index, ] == coded.value, na.rm = TRUE) >0) {power_uncor_temp[i, j] <- 1} else {power_uncor_temp[i, j] <- 0}
      if (sum(test_rows[h_index, ] == coded.value, na.rm = TRUE) >0) {power_rows_temp[i, j] <- 1} else {power_rows_temp[i, j] <- 0}
      if (sum(test_SiZer[h_index, ] == coded.value, na.rm = TRUE) >0) {power_SiZer_temp[i, j] <- 1} else {power_SiZer_temp[i, j] <- 0}
    }
  }
  power <- colSums(power_matrix_temp)/Nsim
  
  power_ms    <- colSums(power_ms_temp)/Nsim
  power_uncor <- colSums(power_uncor_temp)/Nsim
  power_rows  <- colSums(power_rows_temp)/Nsim
  power_SiZer <- colSums(power_SiZer_temp)/Nsim
  
  return(list(power = power, power_ms = power_ms, power_uncor = power_uncor, power_rows = power_rows, power_SiZer = power_SiZer, h.grid = h.grid.new))
}

