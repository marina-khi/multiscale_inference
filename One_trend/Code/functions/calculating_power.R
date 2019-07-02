calculating_power <- function(a1, T, alpha, sigma_eta, Nsim = 1000, SimRuns =1000, type_of_sigma = 'true',
                              q_ = 25, remove.small.ess = 'true', sim.design = 'bump'){
  
  set.seed(1)
  
  spurious_power <- numeric(4)

  # Construct grid
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
  int.plus       <- data.simulated$int.plus
  int.minus      <- data.simulated$int.minus

  if(is.null(int.plus)){
    pos.plus <- matrix(1, ncol=u.len, nrow=h.len)
  } else { 
    temp <- rep(1,N)
    temp[u.left >= int.plus[2] | u.right <= int.plus[1]] <- 100 
    pos.plus <- rep(NA,length(pos.full)) 
    pos.plus[!is.na(pos.full)] <- temp
    pos.plus <- matrix(pos.plus, ncol=u.len, byrow=TRUE)
  }
  if(is.null(int.minus)){
    pos.minus <- matrix(1, ncol=u.len, nrow=h.len) 
  } else {
    temp <- rep(1,N)
    temp[u.left >= int.minus[2] | u.right <= int.minus[1]] <- 100 
    pos.minus <- rep(NA,length(pos.full)) 
    pos.minus[!is.na(pos.full)] <- temp
    pos.minus <- matrix(pos.minus, ncol=u.len, byrow=TRUE)
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
  power_spurious_matrix_temp <- matrix(NA, nrow = Nsim, ncol = 4)
  
  
  set.seed(1)
  for (i in 1:Nsim){
    #Simulating the time series
    data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = sim.design)
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
        
    if (sum(test.res$test_ms * pos.plus == 1, na.rm = TRUE) >0) {power_matrix_temp[i, 1] <- 1} else {power_matrix_temp[i, 1] <- 0}
    if (sum(test.res$test_uncor * pos.plus == 1, na.rm = TRUE) >0) {power_matrix_temp[i, 2] <- 1} else {power_matrix_temp[i, 2] <- 0}
    if (sum(test.res$test_rows * pos.plus == 1, na.rm = TRUE) >0) {power_matrix_temp[i, 3] <- 1} else {power_matrix_temp[i, 3] <- 0}
    if (sum(SiZer_results$test * pos.plus == 1, na.rm = TRUE) >0) {power_matrix_temp[i, 4] <- 1} else {power_matrix_temp[i, 4] <- 0}

    if (sum(test.res$test_ms * pos.plus == 100, na.rm = TRUE) >0) {power_spurious_matrix_temp[i, 1] <- 1} else {power_spurious_matrix_temp[i, 1] <- 0}
    if (sum(test.res$test_uncor * pos.plus == 100, na.rm = TRUE) >0) {power_spurious_matrix_temp[i, 2] <- 1} else {power_spurious_matrix_temp[i, 2] <- 0}
    if (sum(test.res$test_rows * pos.plus == 100, na.rm = TRUE) >0) {power_spurious_matrix_temp[i, 3] <- 1} else {power_spurious_matrix_temp[i, 3] <- 0}
    if (sum(SiZer_results$test * pos.plus == 100, na.rm = TRUE) >0) {power_spurious_matrix_temp[i, 4] <- 1} else {power_spurious_matrix_temp[i, 4] <- 0}
    
  }
  overall_power <- colSums(power_matrix_temp)/Nsim
  spurious_power <- colSums(power_spurious_matrix_temp)/Nsim
  return(list(overall = overall_power, spurious = spurious_power))
}

# 
# calculating_size_for_SiZer <- function(a1, different_T, alpha, sigma_eta, Nsim = 1000, SimRuns =1000){
#   size    <- c()
#   for (T in different_T){
#     set.seed(1)
#     # Construct grid
#     grid      <- grid_construction(T)
#     gset      <- grid$gset
#     u.grid    <- sort(unique(gset[,1]))
#     h.grid    <- sort(unique(gset[,2]))
#     
#     autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))
#     ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
#     deletions <- ess$del
#     grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
#     
#     
#     size_temp <- numeric(Nsim)
#     
#     for (i in 1:Nsim){
#       #Simulating the time series
#       data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
#       data           <- data.simulated$data
#       trend          <- data.simulated$trend
#       sigma_true     <- data.simulated$sigma
#       

#       
#       SiZer_results$test[SiZer_results$test == 2] <- NA
#       if (sum(abs(SiZer_results$test), na.rm = TRUE) >0) {size_temp[i] <- 1} else {size_temp[i] <- 0}
#     }
#     size <- c(size, sum(size_temp)/Nsim)
#     cat("Size of SiZer is", sum(size_temp)/Nsim, " for T = ", T, "\n")
#     rm(size_temp)
#   }
#   return(size.results.sizer = size)
# }
# 
# calculating_size_rowwise <- function(a1, T, alpha, sigma_eta, Nsim = 1000, SimRuns =1000){
#   
#   set.seed(1)
#   # Construct grid
#   grid      <- grid_construction(T)
#   gset      <- grid$gset
#   u.grid    <- sort(unique(gset[,1]))
#   h.grid    <- sort(unique(gset[,2]))
#   
#   autocov   <- (sigma_eta^2/(1-a1^2)) * (a1^seq(0,T-1,by=1))
#   ess       <- ESS.star(u.grid=u.grid, h.grid=h.grid, T=T, autocov=autocov)
#   deletions <- ess$del
#   grid      <- grid_construction(T=T, u.grid=u.grid, h.grid=h.grid, deletions=deletions)
#   
#   h.grid.new <- sort(unique(grid$gset[,2]))
#   
#   size_ms_temp     <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
#   size_uncor_temp  <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
#   size_rows_temp   <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
#   size_SiZer_temp  <- matrix(NA, nrow = Nsim, ncol = length(h.grid.new))
#   
#   # Compute kernel weights and critical value for multiscale test
#   wghts  <- kernel_weights(T=T, grid=grid)
#   quants <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
#   
#   sizer.wghts  <- SiZer_weights(T=T, grid=grid)
#   sizer.std    <- SiZer_std(weights=sizer.wghts, autocov=autocov, T)
#   sizer.quants <- SiZer_quantiles(alpha=alpha, T=T, grid=grid, autocov=autocov)
#   
#   #set.seed(1)
#   for (i in 1:Nsim){
#     #Simulating the time series
#     data.simulated <- simulating_data(T, a1, sigma_eta, sim.design = 'constant')
#     data           <- data.simulated$data
#     trend          <- data.simulated$trend
#     sigma_true     <- data.simulated$sigma
#     
#     stats    <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigma_true, grid=grid)
#     vals     <- stats$values
#     
#     test.res <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)
#     
#     test.res$test_ms[test.res$test_ms == 2] <- NA
#     test.res$test_uncor[test.res$test_uncor == 2] <- NA
#     test.res$test_rows[test.res$test_rows == 2] <- NA        
#     
#     
#     sizer.values  <- sizer.wghts %*% data
#     sizer.vals    <- as.vector(sizer.values)
#     
#     SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)
#     
#     SiZer_results$test[SiZer_results$test == 2] <- NA
#     
#     
#     for (j in 1:length(h.grid.new)){
#       bw  <- h.grid.new[j]
#       h_index <- match(bw, h.grid)
#       
#       if (sum(abs(test.res$test_ms[h_index, ]), na.rm = TRUE) >0) {size_ms_temp[i, j] <- 1} else {size_ms_temp[i, j] <- 0}
#       if (sum(abs(test.res$test_uncor[h_index, ]), na.rm = TRUE) >0) {size_uncor_temp[i, j] <- 1} else {size_uncor_temp[i, j] <- 0}
#       if (sum(abs(test.res$test_rows[h_index, ]), na.rm = TRUE) >0) {size_rows_temp[i, j] <- 1} else {size_rows_temp[i, j] <- 0}
#       if (sum(abs(SiZer_results$test[h_index, ]), na.rm = TRUE) >0) {size_SiZer_temp[i, j] <- 1} else {size_SiZer_temp[i, j] <- 0}
#     }
#   }
#   size_ms    <- colSums(size_ms_temp)/Nsim
#   size_uncor <- colSums(size_uncor_temp)/Nsim
#   size_rows  <- colSums(size_rows_temp)/Nsim
#   size_SiZer <- colSums(size_SiZer_temp)/Nsim
#   
#   return(list(size_ms = size_ms, size_uncor = size_uncor, size_rows = size_rows, size_SiZer = size_SiZer, h.grid = h.grid.new))
# }