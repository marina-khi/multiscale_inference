multiscale_test_rowwise <- function(data, alpha, sigmahat, grid=NULL, SimRuns=100){ 


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
  # psivals   matrix of individual test statistics (hat(psi)/sigmahat), 
  #           each row of the matrix contains the test statistic for a specific bandwidth h  
  # critval   vector of critical values for significance level alpha and each bandwidth h
  # testres   matrix of the test results, 
  #           each row of the matrix contains the test results for a specific bandwidth h  
  #           testres[i,j] = -1: test rejects the null for bandwidth i and location j
  #                              and indicates a decrease in the trend
  #           testres[i,j] = 0:  test accepts the null for bandwidth i and location j
  #           testres[i,j] = 1:  test rejects the null for bandwidth i and location j
  #                              and indicates an increase in the trend


  T <- length(data)

  if(is.null(grid)) {
    GRID <- default_grid(T)
    grid.list <- list()
    grid.list[[1]] <- unique(GRID[,1])
    grid.list[[2]] <- unique(GRID[,2])      
  } else {
    GRID <- expand.grid(u = grid[[1]], h = grid[[2]])  
    grid.list <- grid
  }

  u.grid <- grid.list[[1]]
  h.grid <- grid.list[[2]]


  # Compute critical value of the test 

  psi.null <- matrix(NA,ncol=SimRuns,nrow=length(h.grid))
  critvals <- rep(NA,length(h.grid))
  filename <- paste0("distribution/distr_T_", T, "_rowwise.txt")
  
  if(is.null(grid) & file.exists(filename)) {
    psi.null <- read.table(filename, header=FALSE)
    for(h.pos in 1:length(h.grid))
       critvals[h.pos] <- quantile(psi.null[h.pos,], 1-alpha) 
  } else {
    h.pos <- 1
    for(h.pos in 1:length(h.grid))
    {  cat("","\n")
       cat(paste0("Computing critical value of the rowwise test for h = ", h.grid[h.pos]), "\n")
       progbar  <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")
     
       for(pos in 1:SimRuns)
       {  zeta <- rnorm(T, mean=0, sd=1)
          grid.temp <- list()
          grid.temp[[1]] <- u.grid
          grid.temp[[2]] <- h.grid[h.pos]
          temp <- multiscale_statistic(data=zeta, sigmahat=1, grid=grid.temp)
          psi.null[h.pos,pos] <- temp$max_val
          setTxtProgressBar(progbar, pos)
       }
       close(progbar)

       critvals[h.pos] <- quantile(psi.null[h.pos,], 1-alpha) 
       h.pos <- h.pos+1
    }
    write.table(psi.null, file=filename, row.names=FALSE, col.names=FALSE)    
  } 

   
  # Compute multiscale statistic

  psivals  <- matrix(NA,ncol=length(u.grid),nrow=length(h.grid))
  testvals <- matrix(NA,ncol=length(u.grid),nrow=length(h.grid))
  testres  <- matrix(NA,ncol=length(u.grid),nrow=length(h.grid))
  h.pos    <- 1
 
  for(h.pos in 1:length(h.grid))
  {  grid.temp <- list()
     grid.temp[[1]] <- u.grid
     grid.temp[[2]] <- h.grid[h.pos]

     res <- multiscale_statistic(data=data, sigmahat=sigmahat, grid=grid.temp)

     psivals[h.pos,]  <- as.vector(res$psi_vals)
     testvals[h.pos,] <- as.vector(res$test_vals) 
     testres[h.pos,]  <- (testvals[h.pos,] > critvals[h.pos])
     testres[h.pos,]  <- testres[h.pos,] * sign(psivals[h.pos,])  
  }

  return(list(grid=grid.list,psivals=psivals,testvals=testvals,critvals=critvals,testres=testres))

}