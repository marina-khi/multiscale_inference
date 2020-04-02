multiscale_quantiles <- function(T, grid, SimRuns=1000, probs=seq(0.5,0.995,by=0.005), N_ts = 1){
  # Function that computes the quantiles of the multiscale test statistics under the null.
  #
  # Inputs:
  # T            sample size 
  # grid         grid of location-bandwidth points as produced by the function 'grid_construction',
  #              list with the elements 'gset', 'bws', 'gtype' 
  # probs        vector of probability levels (1-alpha) for which the quantiles are computed 
  #              Default is probs=seq(0.5,0.995,by=0.005).
  # N_ts         number of time series analyzed. Default is 1
  #           
  # Outputs: list with the elements  
  # quant_ms     vector of quantiles for all probability levels in the vector 'probs' 
  #              for our version of the multiscale statistic defined in Section 3

  gset    <- grid$gset
  Psi     <- rep(NA,SimRuns)

  cat("","\n")
  cat("Computing critical value of the multiscale test","\n")
  progbar <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")

  if(N_ts == 1){ 
    for(pos in 1:SimRuns){
      zeta <- rnorm(T, mean=0, sd=1)
      res  <- multiscale_statistics(data=zeta, sigmahat=1, grid=grid) 
      Psi[pos]    <- res$stat_ms
      setTxtProgressBar(progbar, pos)
    }
  }
  else {
    for (pos in 1:SimRuns){
      Psi.ij <- matrix(NA, nrow = N_ts, ncol = N_ts)      
      z_matrix <- matrix(rnorm(T * N_ts, mean=0, sd=1), nrow = T, ncol = N_ts) 
      for(i in 1:N_ts){
        for (j in i:N_ts){
          z_diff <- sqrt(sigmahat_vector_2[i]) * (z_matrix[, i] - mean (z_matrix[, i])) - sqrt(sigmahat_vector_2[j]) * (z_matrix[, j] - mean (z_matrix[, j]))
          Psi.ij[i, j] <- multiscale_statistics(data=z_diff, sigmahat=sqrt(sigmahat_vector_2[i] + sigmahat_vector_2[j]), grid=grid)
        }
      }
      Psi[pos] <- max(Psi.ij, na.rm = TRUE)
      setTxtProgressBar(progbar, pos)
    }
  }
  close(progbar)

  quant    <- as.vector(quantile(Psi, probs=probs))
  quant    <- rbind(probs, quant)

  colnames(quant) <- NULL
  rownames(quant) <- NULL

  return(quant = quant)
}
