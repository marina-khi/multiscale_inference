kernel_weights <- function(T, grid){


  # Function that calculates the kernel weights  
  #
  # Arguments:
  # T          sample size 
  # grid       grid of location-bandwidth points as produced by the function 'grid_construction',
  #            list with the element 'gset' (and possibly others)
  #
  # Outputs: 
  # weights    matrix of kernel weights
  #            w_1(u_1,h_1), ..., w_T(u_1,h_1)
  #            w_1(u_2,h_2), ..., w_T(u_2,h_2)
  #                          ...
  #            w_1(u_N,h_N), ..., w_T(u_N,h_N)


  T     <- as.integer(T) 
  gset  <- grid$gset
  N     <- as.integer(dim(gset)[1])
  gset  <- as.matrix(gset)
  gset  <- as.vector(gset) 
   
  storage.mode(gset) <- "double"

  wghts <- vector(mode = "double", length = N*T)
  
  result <- .C("kernel_weights", T, gset, N, wghts)

  return(matrix(result[[4]],ncol=T,byrow=TRUE))
}

