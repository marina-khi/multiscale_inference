multiscale_statistics <- function(data, weights, sigmahat, grid){


  # Function that calculates the multiscale statistics  
  #
  # Arguments:
  # data        time series of length T 
  # weights     kernel weights matrix produced by the function 'kernel_weights'
  # sigmahat    estimate of sqrt(long-run variance)
  # grid        grid of location-bandwidth points as produced by the function 'grid_construction',
  #             list with the elements 'gset', 'bws'
  #
  # Outputs: list with the elements  
  # values      vector of individual test statistics (hat(psi)/sigmahat) 
  #             for each point (u,h) in the location-bandwidth grid
  # stat_ms     value of the multiscale statistic
  # stat_uncor  value of the multiscale statistic (uncorrected version)
  # stat_rows   vector with length equal to the number of bandwidth levels in the
  #             location-bandwidth grid, stat_rows[i] is the maximum of the individual 
  #             test statistics corresponding to the i-th bandwidth level 
  # stat_order  corrected version of stat_rows, needed to compute the version of the 
  #             multiscale statistics from Section ??


  gset    <- grid$gset
  h.vec   <- grid$bws   
  correct <-  sqrt(2*log(1/(2*gset[,2])))   

  values <- weights %*% data
  values <- as.vector(values/sigmahat)

  stat.ms    <- max(abs(values) - correct)
  stat.uncor <- max(abs(values))
  stat.rows  <- rep(NA,length(h.vec))
  stat.order <- rep(NA,length(h.vec))

  for(ii in 1:length(h.vec))
  {  temp <- (gset[,2] == h.vec[ii])  
     stat.rows[ii] <- max(abs(values[temp]))
     stat.order[ii] <- stat.rows[ii] - sqrt(2*log(1/(2*h.vec[ii]))) 
  }

  return(list(values=values, stat_ms=stat.ms, stat_uncor=stat.uncor, stat_rows=stat.rows, stat_order=stat.order))

}
