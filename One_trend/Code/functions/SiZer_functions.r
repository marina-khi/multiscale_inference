
epan <- function(x)
{  return(0.75*(1-x^2)*((sign(1-x^2)+1)/2))}


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
  # ess.star     matrix with length(u.grid) columns and length(h.grid) rows
  #              specifying ESS.star for each point (u,h)

  N.u <- length(u.grid)
  N.h <- length(h.grid)

  ess.star <- matrix(NA,ncol=N.u,nrow=N.h)
  
  for(i in 1:N.h)
  {  bw <- h.grid[i]

     pos.int <- 1:N.u     
     temp    <- ( u.grid - bw >= 0 & u.grid + bw <= 1 )
     if(sum(temp) > 0)
     { pos.int <- pos.int[temp]
       u       <- u.grid[pos.int[1]]
       arg     <- ((1:T)/T - u)/bw
       ess.star[i,pos.int] <- sum(epan(arg)) / 0.75
     }

     pos.bnd <- 1:N.u
     temp    <- ( u.grid - bw < 0 | u.grid + bw > 1 )
     if(sum(temp) > 0)
     { pos.bnd <- pos.bnd[temp]
       for(j in 1:length(pos.bnd))
       {  u   <- u.grid[pos.bnd[j]]
          arg <- ((1:T)/T - u)/bw
          ess.star[i,pos.bnd[j]] <- sum(epan(arg)) / 0.75
       }
     }
  }

  cov.wghts <- 1 - (1:(T-1))/T
  V.bar     <- autocov[1]/T + (2/T) * sum(cov.wghts * autocov[2:T])
  T.star    <- autocov[1] / V.bar
  ess.star  <- (T.star/T) * ess.star

  deletions <- 1:(N.u*N.h)
  temp      <- as.vector(t(ess.star))
  deletions[temp < 5] <- NA

  return(list(ess=ess.star,del=deletions))
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


SiZer_std <- function(weights, autocov, T)

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

  autocov.mat <- matrix(NA,ncol=T,nrow=T)
  for(ell in 1:(T-1))
     autocov.mat[ell,] <- c(autocov[ell:1],autocov[2:(T-ell+1)])
  autocov.mat[T,] <- autocov[T:1]

  temp     <- autocov.mat %*% t(weights)
  temp     <- t(temp)
  temp     <- weights * temp
  temp     <- temp %*% rep(1,dim(temp)[2])
  std.devs <- sqrt(temp)
  std.devs <- as.vector(std.devs)

  return(std=std.devs)
}


SiZer_quantiles <- function(alpha, T, grid, autocov=NULL)

{ # compute quantiles for SiZer as described in Park et al. (2009), 
  # 'Improved SiZer for time series' 

  gset <- grid$gset
  u.vec <- sort(unique(gset[,1]))
  h.vec <- sort(unique(gset[,2]))

  Delta.tilde <- u.vec[2] - u.vec[1]
  quants      <- rep(0,length(h.vec))

  if(is.null(autocov))
  { for(i in 1:length(h.vec))
    {  gg        <- sum(gset[,2] == h.vec[i])
       theta     <- 2 * pnorm( (sqrt(3)/2) * sqrt(log(gg)) * Delta.tilde/h.vec[i] ) - 1
       arg       <- (1-alpha/2)^(1/(theta*gg))
       quants[i] <- qnorm(arg)
    }
  } 

  
  if(!is.null(autocov))
  { for(i in 1:length(h.vec))
    { gg       <- sum(gset[,2] == h.vec[i])
      arg      <- seq(-(T-1),(T-2),by=1)/(T*h.vec[i])
      autocovs <- c(autocov[T:2],autocov[1:(T-1)])
      int1     <- sum(autocovs * exp(-arg^2/4) * (12 - 12*arg^2 + arg^4) / 16) 
      int2     <- sum(autocovs * exp(-arg^2/4) * (1 - arg^2/2))

      arg      <- seq(-(T-2),(T-1),by=1)/(T*h.vec[i])
      autocovs <- c(autocov[(T-1):2],autocov[1:T])
      int1     <- int1 + sum(autocovs * exp(-arg^2/4) * (12 - 12*arg^2 + arg^4) / 16) 
      int2     <- int2 + sum(autocovs * exp(-arg^2/4) * (1 - arg^2/2))
      I.gamma  <- int1/int2
      
      #integrand_1 <- function(s, h_, delta_, gamma_) {1000* gamma_[floor(s * h_ / delta_) + 1] * exp(-s^2/4) * (12 - 12 * s^2+ s^4)/16}
        
      #I_gamma_num <- 2 * integrate(integrand_1, lower = 0, upper = (T - 1) / (T * h.vec[i]), h_ = h.vec[i], delta_ = 1/T, gamma_ = autocov, subdivisions = 500)[[1]]
        
      #integrand_2 <- function(s, h_, delta_, gamma_) {1000 *gamma_[floor(s * h_ / delta_) + 1] * exp(-s^2/4) * (1 - s^2/2)}
        
      #I_gamma_denom <- 2 * integrate(integrand_2, lower = 0, upper = (T - 1) / (T * h.vec[i]), h_ = h.vec[i], delta_ = 1/T, gamma_ = autocov, subdivisions = 500)[[1]]
      #I.gamma <- I_gamma_num/I_gamma_denom
 
      theta   <- 2 * pnorm( sqrt(I.gamma) * sqrt(log(gg)) * Delta.tilde/h.vec[i] ) - 1
 
      arg       <- (1-alpha/2)^(1/(theta*gg))
      quants[i] <- qnorm(arg)
    }
  }

  return(quants)
}


SiZer_test <- function(values, std.devs, quants, grid){ 

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
  test.sizer[values > critvals]  <- 1
  test.sizer[values < -critvals] <- -1
  
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
  temp    <- sort(unique(as.vector(test.results))) + 2
  temp    <- seq(min(temp),max(temp),by=1)
  col.vec <- col.vec[temp]
  
  image(x=u.grid, y=log(h.grid,10), z=t(test.results), col=col.vec, xlab = '', ylab=expression(log[10](h)), main = plot.title, xaxt = 'n', mgp=c(1,0.5,0))
}