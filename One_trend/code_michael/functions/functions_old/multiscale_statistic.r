multiscale_statistic <- function(data, sigmahat, grid=NULL){

  # Function that calculates the multiscale statistic  
  #
  # Arguments:
  # data       vector of length T 
  # sigmahat   estimate of sqrt(long-run variance)
  # grid       list with the elements 
  #            grid[[1]]: grid of locations
  #            grid[[2]]: grid of bandwidths
  #            If grid=NULL, then a default grid is used.    
  #
  # Outputs: list with the elements  
  # grid.list  location-bandwidth grid under consideration 
  #            (with the same format as the argument 'grid')
  # psi_vals   vector of individual test statistics (hat(psi)/sigmahat) 
  #            for each point (u,h) in the location-bandwidth grid
  # test_vals  vector of individual test statistics (|hat(psi)/sigmahat| - correct) 
  #            for each point (u,h) in the location-bandwidth grid
  # max_val    overall value of the multiscale statistic

  T                  <- as.integer(length(data)) 
  sigmahat           <- as.double(sigmahat)
  storage.mode(data) <- "double"

  if(is.null(grid)) {
    GRID <- default_grid(T)
    grid.list <- list()
    grid.list[[1]] <- unique(GRID[,1])
    grid.list[[2]] <- unique(GRID[,2])    
  } else {
    GRID <- expand.grid(u = grid[[1]], h = grid[[2]])  
    grid.list <- grid
  }

  N                     <- as.integer(dim(GRID)[1]) 
  N.u                   <- as.integer(length(grid.list[[1]]))
  N.h                   <- as.integer(length(grid.list[[2]]))
  correct               <- sqrt(2*log(1/(2*GRID[,2])))
  GRID                  <- as.matrix(GRID)
  GRID                  <- as.vector(GRID) 
  storage.mode(GRID)    <- "double"
  storage.mode(correct) <- "double"
   
  psi_values  <- vector(mode = "double", length = N)
  test_values <- vector(mode = "double", length = N)
  max_value   <- vector(mode = "double", length = 1)
  
  result <- .C("multiscale_statistic", data, T, GRID, correct, N, sigmahat, psi_values, test_values, max_value)

  return(list(grid=grid.list,psi_vals=matrix(result[[7]],nrow=N.h,byrow=TRUE),test_vals=matrix(result[[8]],nrow=N.h,byrow=TRUE),max_val=result[[9]]))
}

