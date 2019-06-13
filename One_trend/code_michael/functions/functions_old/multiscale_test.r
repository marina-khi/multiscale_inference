multiscale_test <- function(data, alpha, sigmahat, grid=NULL, SimRuns=100){ 


  # Function that carries out the multiscale test
  #
  # Arguments:
  # data        vector of length T 
  # alpha       significance level
  # sigmahat    estimate of sqrt(long-run variance)
  # grid        list with the elements 
  #             grid[[1]]: grid of locations
  #             grid[[2]]: grid of bandwidths
  #             If grid=NULL, then a default grid is used.    
  # SimRuns     number of simulations used to compute the critical value of the test, 
  #             default is SimRuns=100 
  #
  # Outputs: list with the elements  
  # grid      location-bandwidth grid under consideration 
  #           (with the same format as the argument 'grid')
  # psivals   vector of individual test statistics (hat(psi)/sigmahat) 
  #           for each point (u,h) in the location-bandwidth grid
  # testvals  vector of individual test statistics (|hat(psi)/sigmahat| - correct) 
  #           for each point (u,h) in the location-bandwidth grid
  # critval   critical value of the test for significance level alpha
  # testres   vector of the test results for each point (u,h) in grid
  #           testres[i] = -1: test rejects the null at grid point i 
  #                            and indicates a decrease in the trend
  #           testres[i] = 0:  test accepts the null at grid point i            
  #           testres[i] = 1:  test rejects the null at grid point i 
  #                            and indicates an increase in the trend


  T <- length(data)


  # Compute critical value of the test 

  filename <- paste0("distribution/distr_T_", T, ".txt")

  if(is.null(grid) & file.exists(filename)) {
    psi.null <- read.table(filename, header=FALSE)
    psi.null <- as.matrix(psi.null)
    psi.null <- as.vector(psi.null)
    critval  <- quantile(psi.null,1-alpha)
  } else {
    cat("","\n")
    cat("Computing critical value of the test","\n")
    progbar  <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")
    psi.null <- rep(NA,SimRuns)
    for(pos in 1:SimRuns)
    {  zeta <- rnorm(T, mean=0, sd=1)
       temp <- multiscale_statistic(data=zeta, sigmahat=1, grid=grid)
       psi.null[pos] <- temp$max_val
       setTxtProgressBar(progbar, pos)
    }
    close(progbar)
    write.table(psi.null, file=filename, row.names=FALSE, col.names=FALSE)
    critval <- quantile(psi.null,1-alpha)
  } 


  # Compute multiscale statistic

  res       <- multiscale_statistic(data=data, sigmahat=sigmahat, grid=grid)
  grid.list <- res$grid
  psivals   <- res$psi_vals
  testvals  <- res$test_vals 
  testres   <- (as.vector(t(testvals)) > critval)
  testres   <- testres * sign(as.vector(t(psivals)))
  testres   <- matrix(testres,nrow=length(grid.list[[2]]),byrow=TRUE)   

  return(list(grid=grid.list,psivals=psivals,testvals=testvals,critval=critval,testres=testres))

}