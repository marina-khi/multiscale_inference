calculating_power <- function(a1, T, alpha, sigma_eta, Nsim = 1000, SimRuns =1000,
                              type_of_sigma = 'true', q_ = 25,
                              remove.small.ess = 'true', sim.design = 'bump', bump.height = 1,
                              region = 'increase',
                              type_of_power = 'global'){
  
  #Construct grid
  grid      <- grid_construction(T)
  gset      <- grid$gset
  u.grid    <- sort(unique(gset[,1]))
  h.grid    <- sort(unique(gset[,2]))
  autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))  
  
  if (remove.small.ess == 'true'){
    ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
    deletions <- ess$del
    grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
  }
    
  gset        <- grid$gset
  N           <- dim(gset)[1]
  gset.full   <- grid$gset_full
  h.grid.full <- sort(unique(gset.full[,2]))
  h.len       <- length(h.grid.full)
  u.len       <- length(unique(gset.full[,1]))
  pos.full    <- grid$pos_full
    
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
  } else if (region == 'all'){
    int.plus  <- data.simulated$int.plus
    int.minus <- data.simulated$int.minus
  } else {cat("Region is not supported \n")}
  
  if(is.null(int)){
    pos <- matrix(1, ncol=u.len, nrow=h.len)
  } else { 
    temp                                       <- rep(1,N)
    temp[u.left >= int[2] | u.right <= int[1]] <- 100 
    pos                                        <- rep(NA,length(pos.full)) 
    pos[!is.na(pos.full)]                      <- temp
    pos                                        <- matrix(pos, ncol=u.len, byrow=TRUE)
  }

  if (type_of_power == 'spurious'){
    coded.value <- 100 * coded.value
  }

  # Compute kernel weights and critical value for multiscale test
  wghts  <- kernel_weights(T=T, grid=grid)
    
  filename = paste0("quantiles/distr_T_", T,".RData")
  if((!file.exists(filename)) || (remove.small.ess == 'true')) {
    quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
    if (remove.small.ess == 'false'){
      save(quants, file = filename)
    }
  } else {
    load(filename)
  }
  
  sizer.wghts  <- SiZer_weights(T=T, grid=grid)
  sizer.std    <- SiZer_std(weights=sizer.wghts, autocov=autocov, T)
  sizer.quants <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
  
  power_matrix_temp          <- matrix(NA, nrow = Nsim, ncol = 4)

  for (i in 1:Nsim){
    #Simulating the time series
    data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = sim.design, slope.fac = bump.height)
    data           <- data.simulated$data
    trend          <- data.simulated$trend
    sigma_true     <- data.simulated$sigma
      
    if (type_of_sigma == 'estimated'){
      #Estimating the coefficients for the ts and the long-run variance
      AR.struc  <- AR_lrv(data=data, q=q_, r.bar=10, p=1)
      sigma_hat <- sqrt(AR.struc$lrv)
      stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid)
    } else if (type_of_sigma == 'true'){
      stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_true, grid=grid)
    }
    
    SiZer.values  <- sizer.wghts %*% data
    sizer.vals    <- as.vector(SiZer.values)
    SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)
      
    vals     <- stats$values
    test.res <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
        
    if (sum(test.res$test_ms * pos == coded.value, na.rm = TRUE) >0) {power_matrix_temp[i, 1] <- 1} else {power_matrix_temp[i, 1] <- 0}
    if (sum(test.res$test_uncor * pos == coded.value, na.rm = TRUE) >0) {power_matrix_temp[i, 2] <- 1} else {power_matrix_temp[i, 2] <- 0}
    if (sum(test.res$test_rows * pos == coded.value, na.rm = TRUE) >0) {power_matrix_temp[i, 3] <- 1} else {power_matrix_temp[i, 3] <- 0}
    if (sum(SiZer_results$test * pos == coded.value, na.rm = TRUE) >0) {power_matrix_temp[i, 4] <- 1} else {power_matrix_temp[i, 4] <- 0}
  }
  power <- colSums(power_matrix_temp)/Nsim
  return(power)
}

calculating_power_rowwise <- function(a1, T, alpha, sigma_eta, Nsim = 1000, SimRuns =1000,
                              type_of_sigma = 'true', q_ = 25,
                              remove.small.ess = 'true', sim.design = 'bump', bump.height = 1,
                              region = 'increase',
                              type_of_power = 'global'){
  
  #Construct grid
  grid      <- grid_construction(T)
  gset      <- grid$gset
  u.grid    <- sort(unique(gset[,1]))
  h.grid    <- sort(unique(gset[,2]))
  autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))  
  
  if (remove.small.ess == 'true'){
    ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
    deletions <- ess$del
    grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
  }
  
  gset        <- grid$gset
  N           <- dim(gset)[1]
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
  } else if (region == 'all'){
    int.plus  <- data.simulated$int.plus
    int.minus <- data.simulated$int.minus
  } else {cat("Region is not supported \n")}
  
  if(is.null(int)){
    pos <- matrix(1, ncol=u.len, nrow=h.len)
  } else { 
    temp                                       <- rep(1,N)
    temp[u.left >= int[2] | u.right <= int[1]] <- 100 
    pos                                        <- rep(NA,length(pos.full)) 
    pos[!is.na(pos.full)]                      <- temp
    pos                                        <- matrix(pos, ncol=u.len, byrow=TRUE)
  }
  
  if (type_of_power == 'spurious'){
    coded.value <- 100 * coded.value
  }
  
  # Compute kernel weights and critical value for multiscale test
  wghts  <- kernel_weights(T=T, grid=grid)
  
  filename = paste0("quantiles/distr_T_", T,".RData")
  if((!file.exists(filename)) || (remove.small.ess == 'true')) {
    quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
    if (remove.small.ess == 'false'){
      save(quants, file = filename)
    }
  } else {
    load(filename)
  }
  
  sizer.wghts  <- SiZer_weights(T=T, grid=grid)
  sizer.std    <- SiZer_std(weights=sizer.wghts, autocov=autocov, T)
  sizer.quants <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
  
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
    
    if (type_of_sigma == 'estimated'){
      #Estimating the coefficients for the ts and the long-run variance
      AR.struc  <- AR_lrv(data=data, q=q_, r.bar=10, p=1)
      sigma_hat <- sqrt(AR.struc$lrv)
      stats     <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_hat, grid=grid)
    } else if (type_of_sigma == 'true'){
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
    
    for (j in 1:length(h.grid.new)){
      bw  <- h.grid.new[j]
      h_index <- match(bw, h.grid)
      
      if (sum(test_ms[h_index, ] == coded.value, na.rm = TRUE) >0) {power_ms_temp[i, j] <- 1} else {power_ms_temp[i, j] <- 0}
      if (sum(test_uncor[h_index, ] == coded.value, na.rm = TRUE) >0) {power_uncor_temp[i, j] <- 1} else {power_uncor_temp[i, j] <- 0}
      if (sum(test_rows[h_index, ] == coded.value, na.rm = TRUE) >0) {power_rows_temp[i, j] <- 1} else {power_rows_temp[i, j] <- 0}
      if (sum(test_SiZer[h_index, ] == coded.value, na.rm = TRUE) >0) {power_SiZer_temp[i, j] <- 1} else {power_SiZer_temp[i, j] <- 0}
    }
  }
  power_ms    <- colSums(power_ms_temp)/Nsim
  power_uncor <- colSums(power_uncor_temp)/Nsim
  power_rows  <- colSums(power_rows_temp)/Nsim
  power_SiZer <- colSums(power_SiZer_temp)/Nsim
  
  return(list(power_ms = power_ms, power_uncor = power_uncor, power_rows = power_rows, power_SiZer = power_SiZer, h.grid = h.grid.new))
}