multiscale_testing <- function(alpha, quantiles, values, grid){ 
  # Function that carries out the multiscale test given that the values of the test  
  # statistics and critical values have already been computed 
  #
  # Arguments:
  # alpha       significance level     
  # quantiles   (1-alpha) quantiles of the multiscale statistic under the null
  #             list with the elements ('quant_ms', 'quant_uncor', 'quant_order', 'quant_rowwise') 
  #             each list element specifies the quantiles for one version of the test 
  # values      values of the individual test statistics for each location-bandwidth point (u,h)
  # grid        grid of location-bandwidth points as produced by the function 'grid_construction',
  #             list with the elements 'gset', 'lens', 'gset_full', 'pos_full' 
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

  T       <- length(data) 
  gset    <- grid$gset
  N.u     <- grid$lens
  correct <- sqrt(2*log(1/(2*gset[,2]))) 


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

  test.results <- matrix(NA,ncol=4,nrow=length(vals1)) 
  # results for multiscale test
  test.results[,1] <- (vals2 > quant.ms) * sign(values)

  gset.full   <- grid$gset_full
  u.grid.full <- unique(gset.full[,1])
  h.grid.full <- unique(gset.full[,2])  
  pos.full    <- grid$pos_full

  test.res <- matrix(2,ncol=4,nrow=length(pos.full))
  for(k in 1:4)
     test.res[!is.na(pos.full), k] <- test.results[,k]

  test.ms    <- matrix(test.res[,1], ncol=length(u.grid.full), byrow=TRUE)

  return(list(ugrid=u.grid.full, hgrid=h.grid.full, test_ms=test.ms, quant.ms = quant.ms))
}