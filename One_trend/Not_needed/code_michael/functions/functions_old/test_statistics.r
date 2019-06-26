test_statistics <- function(data, sigmahat, grid){


  # Function that calculates the multiscale statistic  
  #
  # Arguments:
  # data       vector of length T 
  # sigmahat   estimate of sqrt(long-run variance)
  # grid       grid of location-bandwidth points as produced by the function 'grid_construction',
  #            list with the elements 'gset', 'bws', 'lens' 
  #
  # Outputs: list with the elements  
  # 
  # values       vector of individual test statistics (hat(psi)/sigmahat) 
  #              for each point (u,h) in the location-bandwidth grid
  # stat_ms      value of the multiscale statistic
  # stat_uncor   value of the multiscale statistic (uncorrected version)
  # stat_rows    vector with length equal to the number of bandwidth levels in the
  #              location-bandwidth grid, stat_rows[i] is the maximum of the individual 
  #              test statistics corresponding to the i-th bandwidth level 


  T        <- as.integer(length(data)) 
  sigmahat <- as.double(sigmahat)

  gset  <- grid$gset
  h.vec <- grid$bws   
  N     <- as.integer(dim(gset)[1])
  N.h   <- as.integer(length(h.vec))
  N.u   <- grid$lens
  gset  <- as.matrix(gset)
  gset  <- as.vector(gset) 
  correct <- sqrt(2*log(1/(2*h.vec))) 
   
  storage.mode(data)    <- "double"
  storage.mode(gset)    <- "double"
  storage.mode(h.vec)   <- "double"
  storage.mode(N.u)     <- "integer"
  storage.mode(correct) <- "double"

  values          <- vector(mode = "double", length = N)
  Psi_multiscale  <- vector(mode = "double", length = 1)
  Psi_uncorrected <- vector(mode = "double", length = 1)
  Psi_rowwise     <- vector(mode = "double", length = N.h)
  
  result <- .C("test_statistics", data, T, gset, N, h.vec, N.h, N.u, correct, sigmahat, values, Psi_multiscale, Psi_uncorrected, Psi_rowwise)

  return(list(values=result[[10]], stat_ms=result[[11]], stat_uncor=result[[12]], stat_rows=result[[13]]))

}

