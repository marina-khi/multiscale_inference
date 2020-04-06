multiscale_quantiles <- function(T, grid, sigma_vector = NULL, SimRuns=1000, probs=seq(0.5,0.995,by=0.005), N_ts = 1){
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
  # quant        vector of quantiles for all probability levels in the vector 'probs' 
  #              for our version of the multiscale statistic defined in Section 3

  Psi <- rep(NA,SimRuns)

  #cat("","\n")
  #cat("Computing critical value of the multiscale test","\n")
  #progbar <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")
  tic("Without C method")
  if(N_ts == 1){ 
    for(pos in 1:SimRuns){
      zeta <- rnorm(T, mean=0, sd=1)
      Psi[pos]  <- multiscale_statistics(data=zeta, sigmahat=1, grid=grid)$stat_ms
      #setTxtProgressBar(progbar, pos)
    }
  } else {
    if (is.null(sigma_vector)) {
      sigma_vector <- rep(1, N_ts) #If we did not supply the vector of estimated sigmas, we take vector of ones
    }
    
    #Now we calculate denominators for our statistic based on sqrt(sigma_i^2 + sigma_j^2)
    fun <- function(i,j) sqrt(i*i + j*j)
    sigma_matrix <- outer(sigma_vector, sigma_vector, FUN=fun)
    
    for (pos in 1:SimRuns){
      Psi.ij <- matrix(NA, nrow = N_ts, ncol = N_ts)      
      z_matrix <- matrix(rnorm(T * N_ts, mean=0, sd=1), nrow = T, ncol = N_ts) 
      for(i in 1:(N_ts - 1)){
        for (j in (i + 1):N_ts){
          z_diff <- sigma_vector[i] * (z_matrix[, i] - mean (z_matrix[, i])) - sigma_vector[j] * (z_matrix[, j] - mean (z_matrix[, j]))
          Psi.ij[i, j] <- multiscale_statistics(data=z_diff, sigmahat=sigma_matrix[i, j], grid=grid)$stat_ms
          if ((i == 1) & (j == 2)){
            Psi_tmp <- Psi.ij[i, j]
          } else if (Psi.ij[i, j] > Psi_tmp){
            Psi_tmp <- Psi.ij[i, j]          
          }
        }
      }
      Psi[pos] <- Psi_tmp
      #setTxtProgressBar(progbar, pos)
    }
  }
  #close(progbar)
  toc()
  
  quant    <- as.vector(quantile(Psi, probs=probs))
  quant    <- rbind(probs, quant)

  colnames(quant) <- NULL
  rownames(quant) <- NULL

  return(Psi)
}

# Psi2    <- rep(NA,SimRuns)
# 
# cat("","\n")
# cat("Computing critical value of the multiscale test","\n")
# progbar <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")
# tic("Old method")
# T = T_tempr
# if(N_ts == 1){ 
#   for(pos in 1:SimRuns){
#     zeta     <- rnorm(T, mean=0, sd=1)
#     res      <- multiscale_statistics(data=zeta, sigmahat=1, grid=grid) 
#     Psi[pos] <- res$stat_ms
#     setTxtProgressBar(progbar, pos)
#   }
# } else {
#   wghts <- matrix(kernel_weights(T_tempr, gset_cpp, N), ncol = T_tempr, byrow = TRUE)
#   for (pos in 1:SimRuns){
#     Psi.ij <- matrix(NA, nrow = N_ts, ncol = N_ts)      
#     z_matrix <- matrix(rnorm(T * N_ts, mean=0, sd=1), nrow = T, ncol = N_ts) 
#     for(i in 1:(N_ts-1)){
#       for (j in (i+1):N_ts){
#         z_diff <- sqrt(sigmahat_vector_2[i]) * (z_matrix[, i] - mean (z_matrix[, i])) - sqrt(sigmahat_vector_2[j]) * (z_matrix[, j] - mean (z_matrix[, j]))
#         Psi.ij[i, j] <- multiscale_statistics_old(data=z_diff, weights=wghts, sigmahat=sqrt(sigmahat_vector_2[i] + sigmahat_vector_2[j]), grid=grid)$stat_ms
#       }
#     }
#     Psi2[pos] <- max(Psi.ij, na.rm = TRUE)
#     setTxtProgressBar(progbar, pos)
#   }
# }
# close(progbar)
# toc()