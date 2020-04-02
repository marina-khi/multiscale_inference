multiscale_statistics <- function(data, sigmahat, grid){
  # Function that calculates the multiscale statistics  
  #
  # Arguments:
  # data        time series of length T 
  # sigmahat    estimate of sqrt(long-run variance) for the data
  # grid        grid of location-bandwidth points as produced by the function grid_construction(),
  #             list with at least an element called 'gset'
  #
  # Outputs: list with the elements  
  # values      vector of individual test statistics (hat(psi)/sigmahat) 
  #             for each point (u,h) in the location-bandwidth grid
  # stat_ms     value of the multiscale statistic

  gset    <- grid$gset
  correct <- sqrt(2*log(1/(2*gset[,2])))   
  
  T                      <- as.integer(length(data)) 
  N                      <- as.integer(dim(grid$gset)[1])
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp) 
  storage.mode(gset_cpp) <- "double"
  data                   <- as.vector(data)
  
  values <- kernel_averages(T, gset_cpp, data, sigmahat, N)
  values <- as.vector(values)

  stat.ms <- max(abs(values) - correct)
  
  return(list(values=values, stat_ms=stat.ms))
}