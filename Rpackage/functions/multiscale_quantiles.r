multiscale_quantiles <- function(Tlen, grid, N_ts = 1, sigma_vector = NULL, SimRuns=1000, probs=seq(0.5,0.995,by=0.005)){
  # Function that computes the quantiles of the multiscale test statistics under the null.
  #
  # Inputs:
  # T            sample size 
  # grid         grid of location-bandwidth points as produced by the function 'grid_construction',
  #              list with the elements 'gset', 'bws', 'gtype' 
  # N_ts         number of time series analyzed. Default is 1
  # sigma_vector vector of estimated sqrt(long-run variance) for each time series. If not given, then
  #              the default is a vector of ones of length N_ts (1, ..., 1)
  # probs        vector of probability levels (1-alpha) for which the quantiles are computed 
  #              Default is probs=seq(0.5,0.995,by=0.005).
  #           
  # Outputs: list with the elements  
  # quant        vector of quantiles for all probability levels in the vector 'probs' 
  #              for our version of the multiscale statistic defined in Section 3


  gset                   <- grid$gset
  N                      <- as.integer(dim(grid$gset)[1])
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp) 
  storage.mode(gset_cpp) <- "double"
  

  if(is.null(sigma_vector)){
    sigma_vector <- rep(1, N_ts)
  }
  
  #Compute the quantiles for the multiscale method
  filename = paste0("quantiles/distr_T_", Tlen, "_N_", N_ts, ".RData")
  if(!file.exists(filename)) {
    Phi    <- gaussian_statistics(T = Tlen, N_ts = N_ts, SimRuns = SimRuns, gset = gset_cpp, N = N,
                                  sigma_vec = sigma_vector)
    quant  <- as.vector(quantile(Phi, probs=probs))
    quant  <- rbind(probs, quant)

    colnames(quant) <- NULL
    rownames(quant) <- NULL
    
    save(quant, file = filename)
  } else {
    load(filename)
  }
  return(quant)
}