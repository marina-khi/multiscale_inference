multiscale_test <- function(alpha, data, sigmahat, grid=NULL, kappa=0.1, SimRuns=1000){ 


  # Function that carries out the multiscale test.
  #
  # Arguments:
  # alpha       significance level     
  # data        data vector, time series of length T
  # sigmahat    estimate of sqrt(long-run variance)
  # grid        grid of location-bandwidth points as produced by the function 'grid_construction',
  #             list with the elements 'gset', 'bws', 'lens'. If grid=NULL, a default grid is used.
  # kappa       specifies which order statistic to use for the version of the multiscale test
  #             defined in Section ??
  # SimRuns     number of simulations used to compute the critical value of the test, 
  #             default is SimRuns=1000. 
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
  # test_uncor  matrix of test results for the uncorrected version of the multiscale test 
  # test_order  matrix of test results for the version of the multiscale test from Section ??
  # test_rows   matrix of test results for the rowwise version of the multiscale test 


  T <- length(data) 

  if(is.null(grid))
  { source("functions/grid_construction.r")
    grid <- grid_construction(T)
  }

  gset    <- grid$gset
  gtype   <- grid$gtype
  N.u     <- grid$lens
  correct <- sqrt(2*log(1/(2*gset[,2]))) 


  # Compute kernel weights for local linear derivative estimator

  dyn.load("functions/kernel_weights.so")
  source("functions/kernel_weights.r")
  wghts <- kernel_weights(T=T, grid=grid)


  # Compute multiscale statistics

  source("functions/multiscale_statistics.r")
  results <- multiscale_statistics(data=data, weights=wghts, sigmahat=sigmahat, grid=grid) 
  vals    <- results$values
  

  # Compute critical values of the tests

  filename.ms    <- paste0("quantiles/quant_T_", T, "_multiscale.txt")
  filename.uncor <- paste0("quantiles/quant_T_", T, "_uncorrected.txt")
  filename.order <- paste0("quantiles/quant_T_", T, "_orderstat.txt")
  filename.rows  <- paste0("quantiles/quant_T_", T, "_rowwise.txt")

  if(gtype == "default" & file.exists(filename.ms)){ 
    quant.ms    <- as.matrix(read.table(filename.ms, header=FALSE))
    quant.uncor <- as.matrix(read.table(filename.uncor, header=FALSE))
    quant.order <- as.matrix(read.table(filename.order, header=FALSE))
    quant.rows  <- as.matrix(read.table(filename.rows, header=FALSE))
  }

  if(gtype != "default" | !file.exists(filename.ms)) 
  { source("functions/multiscale_quantiles.r")
    results <- multiscale_quantiles(T=T, grid=grid, weights=wghts, kappa=kappa, SimRuns=SimRuns)
    quant.ms    <- results$quant_ms
    quant.uncor <- results$quant_uncor
    quant.order <- results$quant_order
    quant.rows  <- results$quant_rows

    colnames(quant.ms)    <- NULL
    rownames(quant.ms)    <- NULL
    colnames(quant.uncor) <- NULL
    rownames(quant.uncor) <- NULL
    colnames(quant.order) <- NULL
    rownames(quant.order) <- NULL
    colnames(quant.rows)  <- NULL
    rownames(quant.rows)  <- NULL

    if(gtype == "default" & !file.exists(filename.ms)) 
    { write.table(quant.ms, file=filename.ms, col.names=FALSE, row.names=FALSE)
      write.table(quant.uncor, file=filename.uncor, col.names=FALSE, row.names=FALSE)
      write.table(quant.order, file=filename.order, col.names=FALSE, row.names=FALSE)
      write.table(quant.rows, file=filename.rows, col.names=FALSE, row.names=FALSE)
    }
  }

  quantiles <- list(quant_ms=quant.ms, quant_uncor=quant.uncor, quant_order=quant.order, quant_rows=quant.rows)

  
  # Compute test results

  source("functions/multiscale_testing.r")
  results <- multiscale_testing(alpha=alpha, quantiles=quantiles, values=vals, grid=grid)  


  return(results)
}