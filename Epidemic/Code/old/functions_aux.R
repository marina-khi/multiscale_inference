#' Calculates estimators of the standard deviations for a number of iid samples. 
#' It uses the standard first differences method.
#' @param data      Matrix of the samples. Each sample is located in one of the columns.
#' @param method    Method of estimation, either 'first', then first differences are used, or 'second',
#'                  then the two neighbours method is used, as described in Brown and Levine (2007)
#' @return lrv      Estimator of the long run variance of the error terms.
#' @return ahat     Vector of length p of estimated AR coefficients.
#' @return vareta   Estimator of the variance of the innovation term

estimate_iid_sd <- function(data, method = 'first'){
  N    = ncol(data)
  Tlen = nrow(data)
  sigmahat_vec <- c()
  if (method == 'first') {
    for (i in 1:N){
      variance_i   <- sum((data[2:Tlen, i] - data[1:(Tlen - 1), i])^2)/(2 * Tlen - 2)
      sigma_hat_i  <- sqrt(variance_i)
      sigmahat_vec <- c(sigmahat_vec, sigma_hat_i)
    }
  } else if (method == 'second'){
    for (i in 1:N){
      variance_i   <- 2 * sum((data[3:Tlen, i]/2 - data[2:(Tlen - 1), i] + data[1:(Tlen - 2), i]/2)^2)/(3 * Tlen - 6)
      sigma_hat_i  <- sqrt(variance_i)
      sigmahat_vec <- c(sigmahat_vec, sigma_hat_i)
    }
  } else {
    sigmahat_vec <- NULL
  }
  return(sigmahat_vec)
}

# correction factor for error variance

correct <- function(Y, bw=0.025)
{  Y <- as.matrix(Y)
nn <- dim(Y)[2]
TT <- dim(Y)[1]
X <- (1:TT)/TT
const <- rep(0,nn)
for(i in 1:nn)
{  lambda.fct <- nw(X,Y[,i],bw,TT)
resid <- Y[,i] - lambda.fct
pos <- (lambda.fct > 0)
resid <- resid[pos]/sqrt(lambda.fct[pos])
const[i] <- var(resid)
}   
const <- mean(const)
return(const)
}   

#' Epanechnikov kernel function.
#' @param x A number.
#' @return 3/4(1-x^2) for |x|<=1 and 0 elsewhere.
#' @example 
#' epanechnikov_kernel(1)

epanechnikov_kernel <- function(x)
{
  if (abs(x)<=1)
  {
    result = 3/4 * (1 - x*x)
  } else {
    result = 0
  }
  return(result)
}

#' Function needed for local linear smoothing
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
s_t_1 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

#' Function needed for local linear smoothing
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
s_t_2 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

#' Function needed for local linear smoothing
#' @param u      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
s_t_0 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

#Local Linear estimator using the epanechnikov kernel. 
local_linear_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    s_t_2_value = s_t_2(u, bw, T_size)
    s_t_1_value = s_t_1(u, bw, T_size) 
    for (i in 1:T_size){
      k = (s_t_2_value - s_t_1_value * ((grid_p[i] - u) / bw)) * epanechnikov_kernel((grid_p[i] - u) / bw)
      result = result + k * data_p[i]
      norm = norm + k 
    }
    return(result/norm)
  }
}
