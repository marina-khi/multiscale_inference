multiscale_testing <- function(alpha, data, sigma_vec, grid = NULL,  SimRuns = 1000, N_ts = 1){ 
  # Function that carries out the multiscale test given that the values of the test  
  # statistics and estimatates of long-run variance have already been computed 
  #
  # Arguments:
  # alpha       significance level     
  # grid        grid of location-bandwidth points as produced by the function 'grid_construction',
  #             list with the elements 'gset', 'lens', 'gset_full', 'pos_full' 
  # sigma_vec   vector of estimated long-run variances. Length must be equal to N_ts
  # N_ts        number of time series analysed. Default is 1
  #  
  # Outputs: list with the elements  
  # test_ms     matrix of test results for our multiscale test defined in Section 3
  #             test_ms[i,j] = -1: test rejects the null for the j-th location u and the 
  #                                i-th bandwidth h and indicates a decrease in the trend
  #             test_ms[i,j] = 0:  test does not reject the null for the j-th location u  
  #                                and the i-th bandwidth h 
  #             test_ms[i,j] = 1:  test rejects the null for the j-th location u and the 
  #                                i-th bandwidth h and indicates an increase in the trend
  #             test_ms[i,j] = 2:  no test is carried out at j-th location u and i-th 
  #                                bandwidth h (because the point (u,h) is excluded from  
  #                                the grid as specified by the 'deletions'-option in the
  #                                function 'grid_construction')  
  
  
  if(N_ts == 1){
    Tlen <- length(data) 
  } else {
    Tlen <- nrow(data)
  }
  
  Tlen <- as.integer(Tlen)

  #If grid is not supplied, we construct it by default
  if(is.null(grid)){
    grid <- grid_construction(Tlen)
  }
  
  # Select (1-alpha) quantile of the multiscale statistic under the null
  quantiles <- multiscale_quantiles(Tlen, grid, N_ts, sigma_vector = sigma_vec, SimRuns = SimRuns)

  probs    <- as.vector(quantiles[1,])
  quant    <- as.vector(quantiles[2,])

  if(sum(probs == (1-alpha)) == 0)
    pos <- which.min(abs(probs-(1-alpha)))
  if(sum(probs == (1-alpha)) != 0)
    pos <- which.max(probs == (1-alpha)) 

  quant <- quant[pos]

  # Compute test results
  gset                   <- grid$gset
  N                      <- as.integer(dim(gset)[1])
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp) 
  storage.mode(gset_cpp) <- "double"

  if (N_ts == 1){
    sigma <- as.double(sigma_vec[1])
    Psi = multiscale_statistics(T = Tlen, data = data, gset = gset_cpp, N, sigma)
  } else {
    Psi_ij <- multiscale_statistics_multiple(T = Tlen, N_ts = N_ts, data = data, gset = gset_cpp,
                                 N, sigma_vec)
  }
  # results for multiscale test
  #test.results <- (vals2 > quant.ms) * sign(values)

  #gset.full   <- grid$gset_full
  #u.grid.full <- unique(gset.full[,1])
  #h.grid.full <- unique(gset.full[,2])  
  #pos.full    <- grid$pos_full

  #test.res <- matrix(2,ncol=1,nrow=length(pos.full))
  #test.res[!is.na(pos.full), 1] <- test.results[,1]

  #test.ms    <- matrix(test.res[,1], ncol=length(u.grid.full), byrow=TRUE)

  #return(list(ugrid=u.grid.full, hgrid=h.grid.full, test_ms=test.ms, quant.ms = quant.ms))
  return(list(quant = quant, Psi_ij = Psi_ij))
}