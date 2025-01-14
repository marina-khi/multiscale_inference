#################
#SIZER FUNCTIONS#
#################

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


ESS.star <- function(u.grid, h.grid, T, autocov)
  
{ # compute the effective sample size ESS.star
  #
  # Arguments:
  # u.grid       grid of locations
  # h.grid       grid of bandwidths
  # T            time series length
  # autocov      vector of (estimated) error autocovariances 
  # 
  # Outputs:
  # ess          matrix with length(u.grid) columns and length(h.grid) rows
  #              specifying ESS for each point (u,h)
  # ess.star     matrix with length(u.grid) columns and length(h.grid) rows
  #              specifying ESS.star for each point (u,h)
  # deletions    vector of length length(u.grid)*length(h.grid) with NA elements
  #              in places where ESS.star<5
  
  N.u <- length(u.grid)
  N.h <- length(h.grid)
  
  ess      <- matrix(NA,ncol=N.u,nrow=N.h)
  ess.star <- matrix(NA,ncol=N.u,nrow=N.h)
  
  for(i in 1:N.h){
    bw <- h.grid[i]
    
    pos.int <- 1:N.u     
    temp    <- ( u.grid - bw >= 0 & u.grid + bw <= 1 )
    if(sum(temp) > 0){
      pos.int <- pos.int[temp]
      u       <- u.grid[pos.int[1]]
      arg     <- ((1:T)/T - u)/bw
      ess[i,pos.int] <- sum(epanechnikov_kernel(arg)) / 0.75
    }
    
    pos.bnd <- 1:N.u
    temp    <- ( u.grid - bw < 0 | u.grid + bw > 1 )
    if(sum(temp) > 0){
      pos.bnd <- pos.bnd[temp]
      for(j in 1:length(pos.bnd)){
        u   <- u.grid[pos.bnd[j]]
        arg <- ((1:T)/T - u)/bw
        ess[i,pos.bnd[j]] <- sum(epanechnikov_kernel(arg)) / 0.75
      }
    }
  }
  
  cov.wghts <- 1 - (1:(T-1))/T
  V.bar     <- autocov[1]/T + (2/T) * sum(cov.wghts * autocov[2:T])
  T.star    <- autocov[1] / V.bar
  ess.star  <- (T.star/T) * ess
  
  deletions <- 1:(N.u*N.h)
  temp      <- as.vector(t(ess.star))
  deletions[temp < 5] <- NA
  
  return(list(ess = ess, ess.star=ess.star,del=deletions))
}

SiZer_weights <- function(t_len, grid){
  # calculate the kernel weights for SiZer 
  #
  # Arguments:
  # t_len        sample size 
  # grid         grid of location-bandwidth points as produced by the function 'grid_construction',
  #              list with the element 'gset' (and possibly others)
  #
  # Outputs: 
  # weights      matrix of kernel weights
  #              w_1(u_1,h_1), ..., w_T(u_1,h_1)
  #              w_1(u_2,h_2), ..., w_T(u_2,h_2)
  #                          ...
  #              w_1(u_N,h_N), ..., w_T(u_N,h_N)
  
  t_len <- as.integer(t_len)
  gset  <- grid$gset
  N     <- as.integer(dim(gset)[1])
  gset  <- as.matrix(gset)
  gset  <- as.vector(gset) 
  
  storage.mode(gset) <- "double"
  
  wghts <- vector(mode = "double", length = N*T)
  
  result <- sizer_weights(t_len, gset, N)
  
  return(matrix(result, ncol=t_len, byrow=TRUE))
}


SiZer_std <- function(weights, autocov1, autocov2, t_len)
  
{ # compute local linear derivative estimator and its standard deviation on the
  # location-bandwidth grid.
  #
  # Arguments:
  # data      time series of length T
  # weights   kernel weights matrix produced by the function 'SiZer_weights'
  # autocov   vector of error autocovariances (gamma[0],...,gamma[T-1])
  #
  # Outputs:
  # std       vector of standard deviations (length = number of location-bandwidth
  #           points in the grid)
  
  autocov.mat1 <- matrix(NA, ncol=t_len, nrow=t_len)
  autocov.mat2 <- matrix(NA, ncol=t_len, nrow=t_len)
  
  for(ell in 1:(t_len-1)){
    autocov.mat1[ell,] <- c(autocov1[ell:1],autocov1[2:(t_len-ell+1)])
    autocov.mat2[ell,] <- c(autocov2[ell:1],autocov2[2:(t_len-ell+1)])
  }
  autocov.mat1[t_len,] <- autocov1[t_len:1]
  autocov.mat2[t_len,] <- autocov2[t_len:1]
  
  temp1     <- autocov.mat1 %*% t(weights)
  temp1     <- t(temp1)
  temp1     <- weights * temp1
  temp1     <- temp1 %*% rep(1,dim(temp1)[2])
  temp2     <- autocov.mat2 %*% t(weights)
  temp2     <- t(temp2)
  temp2     <- weights * temp2
  temp2     <- temp2 %*% rep(1,dim(temp2)[2])
  
  std.devs <- sqrt(temp1 + temp2)
  std.devs <- as.vector(std.devs)
  
  return(std=std.devs)
}

SiZer_quantiles <- function(alpha, t_len, grid, autocov1, autocov2)
  
{ # compute quantiles for SiZer as described in Park et al. (2009), 
  # 'Improved SiZer for time series' 
  
  gset  <- grid$gset
  u.vec <- sort(unique(gset[,1]))
  h.vec <- sort(unique(gset[,2]))
  
  Delta.tilde <- u.vec[2] - u.vec[1]
  quants      <- rep(0,length(h.vec))
  
  for(i in 1:length(h.vec)){
    gg        <- sum(gset[,2] == h.vec[i])
    
    integrand_1   <- function(s, h_, delta_, gamma_) {1000 * gamma_[floor(s * h_ / delta_) + 1] * exp(-s^2/4) * (2 - s^2)/8}
    I_gamma_num   <- 2 * integrate(integrand_1, lower = 0, upper = (t_len - 1) / (t_len * h.vec[i]),
                                   h_ = h.vec[i], delta_ = 1/t_len,
                                   gamma_ = autocov1 + autocov2,
                                   subdivisions = 500)[[1]]
    integrand_2   <- function(s, h_, delta_, gamma_) {1000 * gamma_[floor(s * h_ / delta_) + 1] * exp(-s^2/4)}
    I_gamma_denom <- 2 * integrate(integrand_2, lower = 0, upper = (t_len - 1) / (t_len * h.vec[i]),
                                   h_ = h.vec[i], delta_ = 1/t_len,
                                   gamma_ = autocov1 + autocov2,
                                   subdivisions = 500)[[1]]
    I.gamma <- I_gamma_num/I_gamma_denom
    
    # arg       <- seq(-(t_len-1), (t_len-2), by = 1)/(t_len * h.vec[i])
    # autocovs1 <- c(autocov1[t_len:2],autocov1[1:(t_len-1)])
    # autocovs2 <- c(autocov2[t_len:2],autocov2[1:(t_len-1)])
    # int1      <- sum((autocovs1 + autocovs2) * exp(-arg^2/4) * (2 - arg^2) / 8) 
    # int2      <- sum((autocovs1 + autocovs2) * exp(-arg^2/4))
    # 
    # arg       <- seq(-(t_len - 2), (t_len - 1), by = 1)/(t_len * h.vec[i])
    # autocovs1 <- c(autocov1[(t_len - 1):2], autocov1[1:t_len])
    # autocovs2 <- c(autocov2[(t_len - 1):2], autocov2[1:t_len])
    # 
    # int1      <- int1 + sum((autocovs1 + autocovs2) * exp(-arg^2/4) * (2 - arg^2) / 8) 
    # int2      <- int2 + sum((autocovs1 + autocovs2) * exp(-arg^2/4))
    #
    # I.gamma   <- int1/int2
    
    #Clustering coefficient
    theta     <- 2 * pnorm(sqrt(I.gamma) * sqrt(log(gg)) * Delta.tilde/h.vec[i]) - 1
    x         <- (1 - alpha/2)^(1/(theta * gg))
    quants[i] <- qnorm(x)
  }
  return(quants)
}


SiZer_test <- function(values1, values2, std.devs, quants, grid){ 
  
  # carry out row-wise SiZer test
  #
  # Arguments:
  # values1     vector of local linear derivative estimators of the first time series
  #             (length = number of location-bandwidth points in grid)
  # values2     vector of local linear derivative estimators of the second time series
  #             (length = number of location-bandwidth points in grid)
  # std.devs    vector of standard deviations of the local linear derivative estimators
  #             (length = number of location-bandwidth points in grid)
  # quants      vector of quantiles (length = number of bandwidth levels)
  # grid        grid of location-bandwidth points as produced by the function 'grid_construction'
  #
  # Outputs: 
  # test_sizer  matrix of SiZer test results 
  #             test_sizer[i,j] = -1: test rejects the null for the j-th location u and the 
  #                                   i-th bandwidth h and indicates a decrease in the trend
  #             test_sizer[i,j] = 0:  test does not reject the null for the j-th location u  
  #                                   and the i-th bandwidth h 
  #             test_sizer[i,j] = 1:  test rejects the null for the j-th location u and the 
  #                                   i-th bandwidth h and indicates an increase in the trend
  #             test_sizer[i,j] = 2:  no test is carried out at j-th location u and i-th 
  #                                   bandwidth h (because the point (u,h) is excluded from  
  #                                   the grid as specified by the 'deletions'-option in the
  #                                   function 'grid_construction')  
  
  gset    <- grid$gset
  h.vec   <- grid$bws   
  N       <- dim(gset)[1]
  N.h     <- length(h.vec)
  N.u     <- grid$lens
  
  quants   <- rep(quants,N.u)
  critvals <- std.devs * quants
  
  test.sizer <- rep(0,N)
  test.sizer[values1 - values2 > critvals]  <- 1
  test.sizer[values1 - values2 < -critvals] <- -1
  
  gset.full   <- grid$gset_full
  u.grid.full <- unique(gset.full[,1])
  h.grid.full <- unique(gset.full[,2])  
  pos.full    <- grid$pos_full
  
  test.full  <- rep(2,length(pos.full))  
  test.full[!is.na(pos.full)] <- test.sizer
  test.sizer <- matrix(test.full, ncol=length(u.grid.full), byrow=TRUE)
  
  return(list(ugrid=u.grid.full, hgrid=h.grid.full, test=test.sizer))
}


###############
#UCB FUNCTIONS#
###############

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
#' @param x      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
#' @param x_vec  Vector of values for the X variables, length should be T_size
s_t_0_UCB <- function(x, h, T_size, x_vec) {
  result = 0
  for (i in 1:T_size) {
    u = (x_vec[i] - x) / h
    result = result + epanechnikov_kernel(u)
  }
  return(result);
}

#' Function needed for local linear smoothing
#' @param x      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
#' @param x_vec  Vector of values for the X variables, length should be T_size
s_t_1_UCB <- function(x, h, T_size, x_vec) {
  result = 0
  for (i in 1:T_size) {
    u = x_vec[i] - x
    result = result + epanechnikov_kernel(u/h) * u
  }
  return(result);
}

#' Function needed for local linear smoothing
#' @param x      Location at which the local linear smoother is calculated.
#' @param h      Bandwidth that is used for calculating local linear smoothing function.
#' @param T_size Sample size
#' @param x_vec  Vector of values for the X variables, length should be T_size
s_t_2_UCB <- function(x, h, T_size, x_vec) {
  result = 0
  for (i in 1:T_size) {
    u = x_vec[i] - x
    result = result + epanechnikov_kernel(u/h) * u * u
  }
  return(result);
}

#Local Linear estimator using the Epanechnikov kernel. 
UCB_estimation <- function(x_, data_p, grid_p, bw){
  #  if (length(data_p) != length(grid_p)){
  #    cat("Dimensions of the grid and the data do not match, please check the arguments")
  #    return(NULL)
  #  } else {
  result       = 0
  t_len        = length(data_p)
  s_t_2_value1 = s_t_2_UCB(x = x_, h = bw, T_size = t_len, x_vec = grid_p)
  s_t_1_value1 = s_t_1_UCB(x = x_, h = bw, T_size = t_len, x_vec = grid_p)
  s_t_0_value1 = s_t_0_UCB(x = x_, h = bw, T_size = t_len, x_vec = grid_p)
  num1 = s_t_2_value1 * s_t_0_value1 - s_t_1_value1^2
  s_t_2_value2 = s_t_2_UCB(x = x_, h = bw * sqrt(2), T_size = t_len, x_vec = grid_p)
  s_t_1_value2 = s_t_1_UCB(x = x_, h = bw * sqrt(2), T_size = t_len, x_vec = grid_p)
  s_t_0_value2 = s_t_0_UCB(x = x_, h = bw * sqrt(2), T_size = t_len, x_vec = grid_p)
  num2 = s_t_2_value2 * s_t_0_value2 - s_t_1_value2^2
  for (i in 1:t_len){
    u = grid_p[i] - x_
    denom1 = (s_t_2_value1 - s_t_1_value1 * u) * epanechnikov_kernel(u / bw)
    denom2 = (s_t_2_value2 - s_t_1_value2 * u) * epanechnikov_kernel(u / (bw * sqrt(2)))
    result = result + (2 * denom1 / num1 - denom2 / num2) * data_p[i]
  }
  return(result)
}

sigma_estimation_UCB <- function(y_, x_matrix_, beta_est_, m_, k_n_){
  result <- 0
  for (i in 1:(m_ - 1)){
    tmp <- 0
    for (j in 1:k_n_){
      tmp <- tmp + as.vector((y_[j + i * k_n_] - y_[j + (i - 1) * k_n_] - (x_matrix_[j + i * k_n_, ] - x_matrix_[j + (i - 1) * k_n_, ]) %*% as.vector(beta_est_)))
    }
    result <- result + tmp^2
  }
  return(result / (2 * (m_ - 1) * k_n_))
}