multiscale_testing <- function(alpha, data, grid, sigma_vec, N_ts = 1){ 
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

  Tlen    <- length(data) 
  gset    <- grid$gset
  correct <- sqrt(2*log(1/(2*gset[,2]))) 

  #Compute the quantiles for the multiscale method
  filename = paste0("quantiles/distr_T_", Tlen, "_N_", N_ts, ".RData")
  if(!file.exists(filename)) {
    quants <- multiscale_quantiles(T=Tlen, grid=grid, SimRuns=SimRuns, N_ts = N_ts)
    save(quants, file = filename)
  } else {
    load(filename)
  }
  
  
  # Select (1-alpha) quantile of the multiscale statistic under the null
 
  quant.ms    <- quantiles$quant_ms

  probs.ms    <- as.vector(quant.ms[1,])
  quant.ms    <- as.vector(quant.ms[2,])

  if(sum(probs.ms == (1-alpha)) == 0)
    pos.ms <- which.min(abs(probs.ms-(1-alpha)))
  if(sum(probs.ms == (1-alpha)) != 0)
    pos.ms <- which.max(probs.ms == (1-alpha)) 

  quant.ms    <- quant.ms[pos.ms]

  # Compute test results

  vals1 <- abs(values)
  vals2 <- vals1 - correct

  test.results <- matrix(NA,ncol=1,nrow=length(vals2)) 
  # results for multiscale test
  test.results[,1] <- (vals2 > quant.ms) * sign(values)

  gset.full   <- grid$gset_full
  u.grid.full <- unique(gset.full[,1])
  h.grid.full <- unique(gset.full[,2])  
  pos.full    <- grid$pos_full

  test.res <- matrix(2,ncol=1,nrow=length(pos.full))
  test.res[!is.na(pos.full), 1] <- test.results[,1]

  test.ms    <- matrix(test.res[,1], ncol=length(u.grid.full), byrow=TRUE)

  return(list(ugrid=u.grid.full, hgrid=h.grid.full, test_ms=test.ms, quant.ms = quant.ms))
}