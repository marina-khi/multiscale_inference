estimating_sigma_for_AR1 <- function(y, L1, L2){
  #Wrapper of a C function from estimating_sigma.C
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  # y               vector of data length T 
  # L1, L2          tuning parameters
  T <- as.integer(length(y))
  L1 <- as.integer(L1)
  L2 <- as.integer(L2)

  storage.mode(y)    <- "double"
  
  sigmahat         <- vector(mode = "double", length = 1)
  
  result <- .C("estimating_sigma", y, T, L1, L2, sigmahat)
  
  #  cat("kernel = ", result[[5]], ", T=", result[[2]], ", N=", result[[4]], ", sigmahat = ", result[[6]],
  #      "maximum = ", result[[7]], "len of values = ", length(result[[8]]), sep=" ")
  
  return(result[[5]]) 

}