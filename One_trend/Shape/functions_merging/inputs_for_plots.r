inputs_for_plots <- function(test.results, grid){


  # prepares the test results of the multiscale test to prodcue minimal intervals and 
  # SiZer plots
  #
  # Inputs:
  # test results  vector of test results as produced by the function 'multiscale_test',
  #               vector of length N with the entries -1,0 or 1, where N = number of 
  #               location-bandwidth points (u,h) in the grid
  # grid          grid of location-bandwidth points as produced by the function 'grid_construction'
  #               list with the elements 'gset', 'bws', 'lens', 'gset_full', 'pos_full'
  #
  # Outputs: list with the elements
  # mat           matrix of test results
  # cols          colour vector for SiZer plot
  # ugrid         full grid of locations
  # hgrid         full grid of bandwidths
           
  
  gset.full <- grid$gset_full
  u.grid    <- unique(gset.full[,1])
  h.grid    <- unique(gset.full[,2])  
  pos.full  <- grid$pos_full

  res.full  <- rep(2,length(pos.full))
  pos       <- 1

  for(i in 1:length(pos.full))
  {  if(!is.na(pos.full[i]))
     { res.full[i] <- test.results[pos] 
       pos <- pos+1
     }
  }  

  res.mat  <- matrix(res.full,ncol=length(u.grid),byrow=TRUE)
  col.vec <- c('red', 'purple', 'blue', 'gray') 
  temp    <- sort(unique(as.vector(res.mat))) + 2
  col.vec <- col.vec[temp]

  return(list(mat=res.mat, cols=col.vec, ugrid=u.grid, hgrid=h.grid))

}
