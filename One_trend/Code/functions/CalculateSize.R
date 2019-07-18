CalculateSize <- function(T, a1, sigma_eta, different_alpha, Nsim = 1000,
                          SimRuns =1000, sigma.type = 'estimated', q_ = 25, 
                          remove.small.ess =FALSE){
  # Computes the statistical size of the multiscale test (T_ms), the uncorrected version of the multiscale test (T_uc),
  #   the rowwise version of the multiscale test (T_rw) compared to the size of SiZer test (T_SiZer). The 
  #   The size is calculated only for the case where the errors follow AR(1) process.
  #   IMPORTANT: For the computational purposes, the size of SiZer is calculated only for sample sies <= 1000.
  #
  # Args:
  #   T:  Sample size of the time series simulated.
  #   a1: AR(1) coefficient for simulating the error distribution.
  #   sigma_eta: standard deviation of the innovation term in the AR(1) process.
  #   different_alpha: Different confidence levels. Default is 5%.
  #   Nsim: number of simulations for power calculations. Default is 1000.
  #   SimRuns: number of simulations for calculating gaussian quantile for T_ms. Default is 5000.
  #   sigma.type: If 'estimated', then the long-run variance used in the test statistic is first estimated by
  #     AR_lrv() function. Otherwise, if 'true' the true theoretical value of the long-run variance is used.
  #   q_: tuning parameter for estimating the long-run variance from Section 4. Default is 25. 
  #   remove.small.ess: If TRUE, then we restrict attention only to those points (u, h) for which 
  #     the effective sample size for correlated data ESS^*(u, h) is not smaller than 5. Default is FALSE.
  #
  # Returns:
  #   size.ms:    A vector of length equal to the number of different significant levels alpha, each entry corresponding to the size 
  #     of the multiscale testing procedure (T_ms) for this alpha.
  #   size.uncor: A vector of length equal to the number of different significant levels alpha, each entry corresponding to the size 
  #     of the uncorrected version of the multiscale testing procedure (T_uc) for this alpha.
  #   size.rows:  A vector of length equal to the number of different significant levels alpha, each entry corresponding to the size 
  #     of the the rowwise version of the multiscale testing procedure (T_rw) for this alpha.
  #   size.SiZer: A vector of length equal to the number of different significant levels alpha, each entry corresponding to the size 
  #     of the SiZer testing procedure (T_SiZer) for this alpha.
  #   size.rw.ms:    A list of length equal to the number of different significant levels alpha, each entry being a vector of length
  #     equal to the number of bandwidths analysed. The vector contains rowwise size of the multiscale testing procedure (T_ms).
  #   size.rw.uncor: A list of length equal to the number of different significant levels alpha, each entry being a vector of length
  #     equal to the number of bandwidths analysed. The vector contains rowwise size of the uncorrected version
  #     of the multiscale testing procedure (T_uc).
  #   size.rw.rows:  A list of length equal to the number of different significant levels alpha, each entry being a vector of length
  #     equal to the number of bandwidths analysed. The vector contains rowwise size of the multiscale testing procedure (T_rw).
  #   size.rw.SiZer: A list of length equal to the number of different significant levels alpha, each entry being a vector of length
  #     equal to the number of bandwidths analysed. The vector contains rowwise size of the SiZer testing procedure (T_SiZer).
  #   h.grid.new:  A vector of the bandwidths analysed.
  
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
  h.grid.new  <- sort(unique(grid$gset[,2]))
  
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
  
  if (T <= 1000){
    sizer.wghts  <- SiZer_weights(T=T, grid=grid)
    sizer.std    <- SiZer_std(weights=sizer.wghts, autocov=autocov, T)
  
    sizer.quants <- vector("list", length(different_alpha))
    for (k in 1:length(different_alpha)){
      sizer.quants[[k]] <- SiZer_quantiles(alpha=different_alpha[k], T=T, grid=grid, autocov=autocov)
    }
  }
  
  size_matrix_temp_ms    <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
  size_matrix_temp_uncor <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
  size_matrix_temp_rows  <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
  size_matrix_temp_SiZer <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))

  size_rowwise_temp_ms    <- vector("list", length(different_alpha))
  size_rowwise_temp_uncor <- vector("list", length(different_alpha))
  size_rowwise_temp_rows  <- vector("list", length(different_alpha))
  size_rowwise_temp_SiZer <- vector("list", length(different_alpha))
  
  for (j in 1:length(different_alpha)){
    size_rowwise_temp_ms[[j]]    <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
    size_rowwise_temp_uncor[[j]] <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
    size_rowwise_temp_rows[[j]]  <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
    size_rowwise_temp_SiZer[[j]] <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
  }    
  
  for (i in 1:Nsim){
    #Simulating the time series
    data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
    data           <- data.simulated$data
    trend          <- data.simulated$trend
    sigma_true     <- data.simulated$sigma
      
    if (sigma.type == 'estimated'){
      #Estimating the coefficients for the ts and the long-run variance
      AR.struc  <- AR_lrv(data=data, q=q_, r.bar=10, p=1)
      a.hat     <- AR.struc$ahat
      vareta    <- AR.struc$vareta
      sigma_hat <- sqrt(AR.struc$lrv)
      stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid)
    } else if (sigma.type == 'true'){
      stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_true, grid=grid)
    }
    #Values of the multiscale statistic  
    vals     <- stats$values

    if (T <= 1000){
      #Values of SiZer
      values     <- sizer.wghts %*% data
      sizer.vals <- as.vector(values)
    }
    
    #Based on the same values of the test statistic, perform the test at different significance levels
    for (j in 1:length(different_alpha)){
      alpha         <- different_alpha[j]
      test.res      <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
      
      if (sum(abs(test.res$test_ms) == 1, na.rm = TRUE) > 0) {size_matrix_temp_ms[i, j] <- 1} else {size_matrix_temp_ms[i, j] <- 0}
      if (sum(abs(test.res$test_uncor) == 1, na.rm = TRUE) > 0) {size_matrix_temp_uncor[i, j] <- 1} else {size_matrix_temp_uncor[i, j] <- 0}
      if (sum(abs(test.res$test_rows) == 1, na.rm = TRUE) > 0) {size_matrix_temp_rows[i, j] <- 1} else {size_matrix_temp_rows[i, j] <- 0}
      if (T <= 1000){
        SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants[[j]], grid=grid)
        if (sum(abs(SiZer_results$test) == 1, na.rm = TRUE) >0) {size_matrix_temp_SiZer[i, j] <- 1} else {size_matrix_temp_SiZer[i, j] <- 0}
      }
      
      for (k in 1:length(h.grid.new)){
        bw      <- h.grid.new[k]
        h_index <- match(bw, h.grid)
        
        if (sum(abs(test.res$test_ms[h_index, ]) == 1, na.rm = TRUE) > 0) {size_rowwise_temp_ms[[j]][i, k] <- 1} else {size_rowwise_temp_ms[[j]][i, k] <- 0}
        if (sum(abs(test.res$test_uncor[h_index, ]) == 1, na.rm = TRUE) > 0) {size_rowwise_temp_uncor[[j]][i, k] <- 1} else {size_rowwise_temp_uncor[[j]][i, k] <- 0}
        if (sum(abs(test.res$test_rows[h_index, ]) == 1, na.rm = TRUE) > 0) {size_rowwise_temp_rows[[j]][i, k] <- 1} else {size_rowwise_temp_rows[[j]][i, k] <- 0}
        if (T <= 1000){
          if (sum(abs(SiZer_results$test[h_index, ])== 1, na.rm = TRUE) > 0) {size_rowwise_temp_SiZer[[j]][i, k] <- 1} else {size_rowwise_temp_SiZer[[j]][i, k] <- 0}
        }
      }
    }
  }
  
  size_rowwise_ms    <- lapply(size_rowwise_temp_ms, colMeans)
  size_rowwise_uncor <- lapply(size_rowwise_temp_uncor, colMeans)
  size_rowwise_rows  <- lapply(size_rowwise_temp_rows, colMeans)
  size_rowwise_SiZer <- lapply(size_rowwise_temp_SiZer, colMeans)

  size_ms    <- colSums(size_matrix_temp_ms)/Nsim
  size_uncor <- colSums(size_matrix_temp_uncor)/Nsim
  size_rows  <- colSums(size_matrix_temp_rows)/Nsim
  size_SiZer <- colSums(size_matrix_temp_SiZer)/Nsim
  
  rm(size_matrix_temp_ms, size_matrix_temp_uncor, size_matrix_temp_rows, size_matrix_temp_SiZer)

  return(list(size.ms = size_ms, size.uncor = size_uncor, size.rows = size_rows, size.SiZer = size_SiZer,
              size.rw.ms = size_rowwise_ms, size.rw.uncor = size_rowwise_uncor, size.rw.rows = size_rowwise_rows, size.rw.SiZer = size_rowwise_SiZer,
              h.grid = h.grid.new))
}