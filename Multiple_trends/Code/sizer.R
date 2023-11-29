rm(list=ls())

library(MSinference)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
library(Rcpp)


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
      ess[i,pos.int] <- sum(epan(arg)) / 0.75
    }
  
    pos.bnd <- 1:N.u
    temp    <- ( u.grid - bw < 0 | u.grid + bw > 1 )
    if(sum(temp) > 0){
      pos.bnd <- pos.bnd[temp]
      for(j in 1:length(pos.bnd)){
        u   <- u.grid[pos.bnd[j]]
        arg <- ((1:T)/T - u)/bw
        ess[i,pos.bnd[j]] <- sum(epan(arg)) / 0.75
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

SiZer_weights <- function(T, grid)
  
{ # calculate the kernel weights for SiZer 
  #
  # Arguments:
  # T            sample size 
  # grid         grid of location-bandwidth points as produced by the function 'grid_construction',
  #              list with the element 'gset' (and possibly others)
  #
  # Outputs: 
  # weights      matrix of kernel weights
  #              w_1(u_1,h_1), ..., w_T(u_1,h_1)
  #              w_1(u_2,h_2), ..., w_T(u_2,h_2)
  #                          ...
  #              w_1(u_N,h_N), ..., w_T(u_N,h_N)
  
  T     <- as.integer(T) 
  gset  <- grid$gset
  N     <- as.integer(dim(gset)[1])
  gset  <- as.matrix(gset)
  gset  <- as.vector(gset) 
  
  storage.mode(gset) <- "double"
  
  wghts <- vector(mode = "double", length = N*T)
  
  result <- sizer_weights(T, gset, N)
  
  return(matrix(result,ncol=T,byrow=TRUE))
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
  # values    vector of local linear derivative estimators (length = number of 
  #           location-bandwidth points in the grid) 
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
    arg       <- seq(-(t_len-1),(t_len-2), by = 1)/(t_len * h.vec[i])
    autocovs1 <- c(autocov1[t_len:2],autocov1[1:(t_len-1)])
    autocovs2 <- c(autocov1[t_len:2],autocov1[1:(t_len-1)])
    int1      <- sum((autocovs1 + autocovs2) * exp(-arg^2/4) * (2 - arg^2) / 8) 
    int2      <- sum((autocovs1 + autocovs2) * exp(-arg^2/4))
    I.gamma   <- int1/int2
    
    #Clustering coefficient
    theta     <- 2 * pnorm( sqrt(I.gamma) * sqrt(log(gg)) * Delta.tilde/h.vec[i] ) - 1
    x         <- (1-alpha/2)^(1/(theta*gg))
    quants[i] <- qnorm(x)
  }
  return(quants)
}


SiZer_test <- function(values1, values2, std.devs, quants, grid){ 
  
  # carry out row-wise SiZer test
  #
  # Arguments:
  # values      vector of local linear derivative estimators (length = number of location-
  #             bandwidth points in grid)
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

SiZermap <- function(u.grid, h.grid, test.results, plot.title = NA)
  
{  # computes SiZer map from the test results 
  
  col.vec <- c("red", "purple", "blue", "gray") 
  #col.vec <- c("#F7F7F7", "#969696", "#525252", "#636363") 
  temp    <- sort(unique(as.vector(test.results))) + 2
  temp    <- seq(min(temp),max(temp),by=1)
  col.vec <- col.vec[temp]
  
  image(x=u.grid, y=log(h.grid,10), z=t(test.results), col=col.vec, xlab = '', ylab=expression(log[10](h)), main = plot.title, xaxt = 'n', mgp=c(1,0.5,0))
}


#Load necessary functions  
# source("functions/ConstructGrid.r")
# source("functions/multiscale_statistics.r")
# source("functions/multiscale_quantiles.r")
# source("functions/multiscale_testing.r")
# source("functions/long_run_variance.r")
# source("functions/sim.r")

sourceCpp("functions/SiZer_functions.cpp")


##############################
#Defining necessary constants#
##############################
#set.seed(12321)

n_ts  <- 2 #Number of time series
t_len <- 500

#n_rep    <- 5000 #number of simulations for calculating size and power
#sim_runs <- 1000 #number of simulations to calculate the Gaussian quantiles

#different_T     <- c(100, 200, 350, 500, 1000) #Different lengths of time series
different_alpha <- c(0.01, 0.05, 0.1) #Different confidence levels
#different_b     <- c(0) #Zero is for calculating the size

#For the covariate process
#beta    <- c(1, 1, 1) 
#a_x_vec <- c(0.25, 0.25, 0.25) #VAR(1) coefficients
#phi     <- 0              #dependence between the innovations

#For the error process
a     <- 0.25
sigma <- 0.25

#For the fixed effects
#rho <- 0.25 #covariance between the fixed effects

#Parameters for the estimation of long-run-variance
q <- 25 
r <- 15

#For parallel computation
#numCores  <- round(parallel::detectCores() * .70)

#Construct grid
grid      <- construct_grid(t_len)
gset      <- grid$gset
u.grid    <- sort(unique(gset[,1]))
h.grid    <- sort(unique(gset[,2]))

errors1 <- arima.sim(model = list(ar = a),
                     innov = rnorm(t_len, 0, sigma),
                     n = t_len)

errors2 <- arima.sim(model = list(ar = a),
                     innov = rnorm(t_len, 0, sigma),
                     n = t_len)

b = 3
m_matrix      <- matrix(0, nrow = t_len, ncol = n_ts)
m_matrix[, 1] <- (1:t_len - 0.5 * t_len) * (b / t_len)

timeseries1 <- m_matrix[, 1] + errors1
timeseries2 <- m_matrix[, 2] + errors2


AR.struc1   <- estimate_lrv(data = timeseries1, q = q, r_bar = r, p = 1)
sigma_hat_1 <- sqrt(AR.struc1$lrv)
a_1         <- AR.struc1$ahat
autocov1    <- (sigma_hat_1^2/(1-a_1^2)) * (a_1^seq(0,t_len-1,by=1))  

AR.struc2   <- estimate_lrv(data = timeseries2, q = q, r_bar = r, p = 1)
sigma_hat_2 <- sqrt(AR.struc2$lrv)
a_2         <- AR.struc2$ahat
autocov2    <- (sigma_hat_2^2/(1-a_2^2)) * (a_2^seq(0,t_len-1,by=1))  

gset        <- grid$gset
N           <- as.integer(dim(gset)[1])
h.grid.new  <- sort(unique(grid$gset[,2]))

sizer.wghts  <- SiZer_weights(T=t_len, grid=grid)
sizer.std    <- SiZer_std(weights = sizer.wghts, autocov1 = autocov1,
                          autocov2 = autocov2, t_len = t_len)

sizer.quants <- vector("list", length(different_alpha))
for (k in 1:length(different_alpha)){
  sizer.quants[[k]] <- SiZer_quantiles(alpha=different_alpha[k], t_len =t_len,
                                       grid=grid, autocov1=autocov1,
                                       autocov2=autocov2)
}

#Values of SiZer
values1     <- sizer.wghts %*% timeseries1
sizer.vals1 <- as.vector(values1)
values2     <- sizer.wghts %*% timeseries2
sizer.vals2 <- as.vector(values2)

#Based on the same values of the test statistic, perform the test at different significance levels
for (j in 1:length(different_alpha)){
  alpha         <- different_alpha[j]

  SiZer_results <- SiZer_test(values1=sizer.vals1, values2 = sizer.vals2,
                              std.devs=sizer.std, quants=sizer.quants[[j]],
                              grid=grid)
  SiZermap(u.grid, h.grid, SiZer_results$test, plot.title = NA)
}

