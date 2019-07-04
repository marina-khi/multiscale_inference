calculating_size <- function(a1, different_T, different_alpha, sigma_eta, Nsim = 1000, kappa =0.1,
                             SimRuns =1000, type_of_sigma = 'estimated', q_ = 25, remove.small.ess = 'false'){
  size_ms    <- numeric(length(different_alpha))
  size_uncor <- numeric(length(different_alpha))
  size_rows  <- numeric(length(different_alpha))
  
  for (T in different_T){
    # Construct grid
    grid      <- grid_construction(T)
    gset      <- grid$gset
    u.grid    <- sort(unique(gset[,1]))
    h.grid    <- sort(unique(gset[,2]))
    
    if (remove.small.ess == 'true'){
      autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))
      ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
      deletions <- ess$del
      grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
    }
    
    # Compute kernel weights and critical value for multiscale test
    wghts  <- kernel_weights(T=T, grid=grid)
    
    filename = paste0("quantiles/distr_T_", T,".RData")
    if((!file.exists(filename)) || (remove.small.ess == 'true')) {
      quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=kappa, SimRuns=SimRuns)
      if (remove.small.ess == 'false'){
        save(quants, file = filename)
      }
    } else {
      load(filename)
    }
    
    size_matrix_temp       <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
    size_matrix_temp_uncor <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
    size_matrix_temp_rows  <- matrix(NA, nrow = Nsim, ncol = length(different_alpha))
    
    set.seed(321)     
    for (i in 1:Nsim){
      #Simulating the time series
      data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
      data           <- data.simulated$data
      trend          <- data.simulated$trend
      sigma_true     <- data.simulated$sigma
      
      if (type_of_sigma == 'estimated'){
        #Estimating the coefficients for the ts and the long-run variance
        AR.struc  <- AR_lrv(data=data, q=q_, r.bar=10, p=1)
        a.hat     <- AR.struc$ahat
        vareta    <- AR.struc$vareta
        sigma_hat <- sqrt(AR.struc$lrv)
        stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid)
      } else if (type_of_sigma == 'true'){
        stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_true, grid=grid)
      }
      
      vals     <- stats$values
      
      #Based on the same values of the test statistic, perform the test at dife
      for (j in 1:length(different_alpha)){
        alpha    <- different_alpha[j]
        test.res <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
        
        test.res$test_ms[test.res$test_ms == 2] <- NA
        test.res$test_uncor[test.res$test_uncor == 2] <- NA
        test.res$test_rows[test.res$test_rows == 2] <- NA        
        
        if (sum(abs(test.res$test_ms), na.rm = TRUE) >0) {size_matrix_temp[i, j] <- 1} else {size_matrix_temp[i, j] <- 0}
        if (sum(abs(test.res$test_uncor), na.rm = TRUE) >0) {size_matrix_temp_uncor[i, j] <- 1} else {size_matrix_temp_uncor[i, j] <- 0}
        if (sum(abs(test.res$test_rows), na.rm = TRUE) >0) {size_matrix_temp_rows[i, j] <- 1} else {size_matrix_temp_rows[i, j] <- 0}
      }
    }
    size_ms <- cbind(size_ms, colSums(size_matrix_temp)/Nsim)
    size_uncor <- cbind(size_uncor, colSums(size_matrix_temp_uncor)/Nsim)
    size_rows <- cbind(size_rows, colSums(size_matrix_temp_rows)/Nsim)
    
    #rm(size_matrix_temp, size_matrix_temp_uncor, size_matrix_temp_rows)
  }
  
  return(list(size.results.ms = t(size_ms[, -1]),size.results.uncor = t(size_uncor[, -1]),size.results.rows = t(size_rows[, -1]) ))
}


calculating_size_for_SiZer <- function(a1, different_T, alpha, sigma_eta, Nsim = 1000, SimRuns =1000){
  size    <- c()
  
  for (T in different_T){
    # Construct grid
    grid      <- grid_construction(T)
    gset      <- grid$gset
    u.grid    <- sort(unique(gset[,1]))
    h.grid    <- sort(unique(gset[,2]))

    autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))
    ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
    deletions <- ess$del
    grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)

    sizer.wghts  <- SiZer_weights(T=T, grid=grid)
    sizer.std    <- SiZer_std(weights=sizer.wghts, autocov=autocov, T)
    sizer.quants <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
    
    size_temp <- numeric(Nsim)
    
    set.seed(321)
    for (i in 1:Nsim){
      #Simulating the time series
      data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
      data           <- data.simulated$data
      trend          <- data.simulated$trend
      sigma_true     <- data.simulated$sigma
    
      values     <- sizer.wghts %*% data
      sizer.vals <- as.vector(values)
      
      SiZer_results    <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)
      
      SiZer_results$test[SiZer_results$test == 2] <- NA
      if (sum(abs(SiZer_results$test), na.rm = TRUE) >0) {size_temp[i] <- 1} else {size_temp[i] <- 0}
    }
    size <- c(size, sum(size_temp)/Nsim)
    cat("Size of SiZer is", sum(size_temp)/Nsim, " for T = ", T, "\n")
    rm(size_temp)
  }
  return(size.results.sizer = size)
}

calculating_size_rowwise <- function(a1, T, alpha, sigma_eta, Nsim = 1000, SimRuns =1000){
  # Construct grid

  grid      <- grid_construction(T)
  gset      <- grid$gset
  u.grid    <- sort(unique(gset[,1]))
  h.grid    <- sort(unique(gset[,2]))
    
  autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))
  ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
  deletions <- ess$del
  grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
  
  h.grid.new <- sort(unique(grid$gset[,2]))
  
  size_ms_temp     <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
  size_uncor_temp  <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
  size_rows_temp   <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
  size_SiZer_temp  <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))

  # Compute kernel weights and critical value for multiscale test
  wghts  <- kernel_weights(T=T, grid=grid)
  quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)

  sizer.wghts  <- SiZer_weights(T=T, grid=grid)
  sizer.std    <- SiZer_std(weights=sizer.wghts, autocov=autocov, T)
  sizer.quants <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
  
  set.seed(321)
  for (i in 1:Nsim){
    #Simulating the time series
    data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
    data           <- data.simulated$data
    trend          <- data.simulated$trend
    sigma_true     <- data.simulated$sigma
      
    stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_true, grid=grid)
    vals     <- stats$values
      
    test.res <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
        
    test.res$test_ms[test.res$test_ms == 2] <- NA
    test.res$test_uncor[test.res$test_uncor == 2] <- NA
    test.res$test_rows[test.res$test_rows == 2] <- NA        

    
    sizer.values  <- sizer.wghts %*% data
    sizer.vals    <- as.vector(sizer.values)
    
    SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)
    
    SiZer_results$test[SiZer_results$test == 2] <- NA
    
    
    for (j in 1:length(h.grid.new)){
      bw  <- h.grid.new[j]
      h_index <- match(bw, h.grid)
        
      if (sum(abs(test.res$test_ms[h_index, ]), na.rm = TRUE) >0) {size_ms_temp[i, j] <- 1} else {size_ms_temp[i, j] <- 0}
      if (sum(abs(test.res$test_uncor[h_index, ]), na.rm = TRUE) >0) {size_uncor_temp[i, j] <- 1} else {size_uncor_temp[i, j] <- 0}
      if (sum(abs(test.res$test_rows[h_index, ]), na.rm = TRUE) >0) {size_rows_temp[i, j] <- 1} else {size_rows_temp[i, j] <- 0}
      if (sum(abs(SiZer_results$test[h_index, ]), na.rm = TRUE) >0) {size_SiZer_temp[i, j] <- 1} else {size_SiZer_temp[i, j] <- 0}
    }
  }
  size_ms    <- colSums(size_ms_temp)/Nsim
  size_uncor <- colSums(size_uncor_temp)/Nsim
  size_rows  <- colSums(size_rows_temp)/Nsim
  size_SiZer <- colSums(size_SiZer_temp)/Nsim

  return(list(size_ms = size_ms, size_uncor = size_uncor, size_rows = size_rows, size_SiZer = size_SiZer, h.grid = h.grid.new))
}