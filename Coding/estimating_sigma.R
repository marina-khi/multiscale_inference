estimating_sigma <- function(data){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  sigma2 = sigma_eta2 / ((1 - sum(a_hat))^2)
  return(sqrt(sigma2))
}