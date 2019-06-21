multiscale_test <- function(alpha, data, weights, sigmahat, grid, kappa=0.1, SimRuns=1000){ 


  # Function that carries out the multiscale test.
  #
  # Arguments:
  # alpha       significance level     
  # data        data vector 
  # weights     kernel weights matrix produced by the function 'kernel_weights'
  # sigmahat    estimate of sqrt(long-run variance)
  # grid        grid of location-bandwidth points as produced by the function 'grid_construction',
  #             list with the elements 'gset', 'bws', 'lens' 
  # kappa       specifies which order statistic to use for the version of the multiscale test
  #             defined in Section ??
  # SimRuns     number of simulations used to compute the critical value of the test, 
  #             default is SimRuns=1000. 
  #
  # Outputs: list with the elements  
  # test_ms      vector of test results for our multiscale test defined in Section 3
  #              test_ms[i] = -1: test rejects the null for the i-th point (u,h) in the 
  #                               location-bandwidth grid (grid$gset) and indicates a 
  #                               decrease in the trend
  #              test._ms[i] = 0: test does not reject the null for the i-th point (u,h)  
  #                               in the location-bandwidth grid (grid$gset) 
  #              test_ms[i] = 1:  test rejects the null for the i-th point (u,h) in the 
  #                               location-bandwidth grid (grid$gset) and indicates an 
  #                               increase in the trend
  # test_uncor   vector of test results for the uncorrected version of the multiscale test 
  # test_order   vector of test results for the version of the multiscale test from Section ??
  # test_rows    vector of test results for the rowwise version of the multiscale test 


  T       <- length(data) 
  gset    <- grid$gset
  h.vec   <- grid$bws   
  N       <- dim(gset)[1]
  N.h     <- length(h.vec)
  N.u     <- grid$lens
  correct <- sqrt(2*log(1/(2*gset[,2]))) 


  # Compute multiscale statistics

  res  <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigmahat, grid=grid) 
  vals <- res$values
  

  # Compute critical value of the test 

  res   <- critical_values(T=T, grid=grid, weights=weights, kappa=kappa, SimRuns=SimRuns)
  probs <- res$probs

  if(sum(probs == (1-alpha)) == 0)
    pos <- which.min(abs(probs-(1-alpha)))
  if(sum(probs == (1-alpha)) != 0)
    pos <- which.max(probs == (1-alpha)) 

  critval.ms    <- res$quant_ms[pos]  
  critval.uncor <- res$quant_uncor[pos]  
  critval.order <- res$quant_order[pos]  
  critval.rows  <- res$quant_rows[,pos]
  critval.rows  <- rep(critval.rows,N.u)


  # Compute test results

  vals1 <- abs(vals)
  vals2 <- vals1 - correct

  test.ms    <- (vals2 > critval.ms) * sign(vals)
  test.uncor <- (vals1 > critval.uncor) * sign(vals)
  test.order <- (vals2 > critval.order) * sign(vals)
  test.rows  <- (vals1 > critval.rows) * sign(vals)

  return(list(test_ms=test.ms, test_uncor=test.uncor, test_order=test.order, test_rows=test.rows))

}