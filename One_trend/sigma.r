rm(list=ls())

p      <- 1
a_star <- 0.95
T      <- 200
Nsim   <- 500


# Functions

#dyn.load("sigma_HvK.so")
#source("sigma_HvK.r")

emp_acf <- function(y_data,ell,p)
{   # computes autocovariances at lags 0 to p 
    # for the ell-th differences of y_data
    y_diff   <- diff(y_data,ell)
    len      <- length(y_diff)
    autocovs <- rep(0,p+1)
    for(k in 1:(p+1))
       autocovs[k] <- sum(y_diff[k:len]*y_diff[1:(len-k+1)])/len
    return(autocovs)
}

sigma_eta <- function(y_data,coefs,p)
{   y_diff <- diff(y_data)
    len    <- length(y_diff)
    res    <- rep(0,len)
    for(i in (p+1):len)
       res[i] <- y_diff[i] - sum(coefs * y_diff[(i-1):(i-p)])
    var_eta <- mean(res^2)/2
    return(var_eta)
}

corrections <- function(coefs,var_eta,len)
{   p <- length(coefs)  
    c_vec <- rep(0,len)
    c_vec[1] <- 1
    for(j in 2:len)
    { lags <- (j-1):max(j-p,1)
      c_vec[j] <- sum(coefs[1:length(lags)] * c_vec[lags])  
    }
    c_vec <- c_vec * var_eta
    return(c_vec)
}

AR_coef <- function(y_data,L1,L2,correct,p)
{   pos <- 0
    a_mat <- matrix(0,nrow=L2-L1+1,ncol=p)
    for(ell in L1:L2)
    {  pos <- pos + 1
       autocovs <- emp_acf(y_data,ell,p)
       cov_mat <- matrix(0,ncol=p,nrow=p)
       for(i in 1:p)
       {  for(j in 1:p)
             cov_mat[i,j] <- autocovs[abs(i-j)+1]
       }
       cov_vec <- autocovs[2:(p+1)] + correct[(ell-1+1):(ell-p+1)]
       a_mat[pos,] <- solve(cov_mat) %*% cov_vec 
    }    
    a_hat <- colMeans(a_mat)
    a_hat <- as.vector(a_hat)
    return(a_hat)
}

AR_coef_HvK <- function(y_data,L1,L2,p)
{   pos <- 0
    var_vec <- rep(0,L2-L1+1)
    for(ell in L1:L2)
    {  pos <- pos + 1
       y_diff <- diff(y_data,lag=ell)
       var_vec[pos] <- mean(y_diff^2)
    }
    g0 <- mean(var_vec)/2

    g_vec <- rep(0,p)
    for(ell in 1:p)
    {  y_diff <- diff(y_data,lag=ell)
       g_vec[ell] <- g0 - mean(y_diff^2)/2
    }
    g_vec <- c(g0,g_vec)

    cov_mat <- matrix(0,ncol=p,nrow=p)
    for(i in 1:p)
    {  for(j in 1:p)
          cov_mat[i,j] <- g_vec[abs(i-j)+1]
    }
    cov_vec <- g_vec[2:(p+1)] 
  
    a_hat_HvK <- solve(cov_mat) %*% cov_vec 
    a_hat_HvK <- as.vector(a_hat_HvK)
    return(a_hat_HvK)
}


# Simulations

a_vec1 <- a_vec2 <- a_vec_HvK <- a_vec_oracle <- rep(NA,Nsim)
sig_vec1 <- sig_vec2 <- sig_vec_HvK <- sig_vec_oracle <- rep(NA,Nsim)
a_vec1 <- a_vec2 <- a_vec_HvK <- a_vec_oracle <- rep(NA,Nsim)
lrv_vec1 <- lrv_vec2 <- lrv_vec_HvK <- lrv_vec_oracle <- rep(NA,Nsim)

for(nb in 1:Nsim){
  set.seed(nb)

  # simulate data

  burn  <- 50
  eta   <- rnorm(T+burn)
  eps   <- rep(0,T + burn)
  for(i in (p+1):(T+burn))
     eps[i] <- eps[(i-1):(i-p)] %*% a_star + eta[i]
  eps <- eps[-(1:burn)]

  slope <- 4
  trend <- slope * (1:T)/T

  y_data <- trend + eps

  # compute estimators

  K1 <- p+1
  K2 <- 10
  L1 <- 20
  L2 <- 30

  a_hat1       <- AR_coef(y_data,L1,L2,rep(0,L2),p)
  sig_eta_hat1 <- sigma_eta(y_data,a_hat1,p)
  lrv_hat1     <- sig_eta_hat1/(1-sum(a_hat1))^2 
  correct      <- corrections(a_hat1,sig_eta_hat1,K2+1)
  a_hat2       <- AR_coef(y_data,K1,K2,correct,p)
  sig_eta_hat2 <- sigma_eta(y_data,a_hat2,p)
  lrv_hat2     <- sig_eta_hat2/(1-sum(a_hat2))^2 
  a_vec1[nb]   <- a_hat1
  a_vec2[nb]   <- a_hat2
  sig_vec1[nb] <- sig_eta_hat1
  sig_vec2[nb] <- sig_eta_hat2
  lrv_vec1[nb] <- lrv_hat1
  lrv_vec2[nb] <- lrv_hat2
 
  a_vec_HvK[nb]   <- AR_coef_HvK(y_data,L1,L2,p)
  sig_vec_HvK[nb] <- sigma_eta(y_data,a_vec_HvK[nb],p)
  lrv_vec_HvK[nb] <- sig_vec_HvK[nb]/(1-sum(a_vec_HvK[nb]))^2
  #res_HvK <- sigma_HvK(y_data,L1,L2)
  #a_vec_HvK[nb] <- res_HvK[[2]] 

  res_oracle         <- lm(eps[2:T] ~ eps[1:(T-1)] - 1)
  a_vec_oracle[nb]   <- as.vector(res_oracle$coef)
  sig_vec_oracle[nb] <- mean(res_oracle$residuals^2)
  lrv_vec_oracle[nb] <- sig_vec_oracle[nb]/(1-sum(a_vec_oracle[nb]))^2 

  print(nb)
}


# Plots

bounds <- "no"

a_min <- min(c(a_vec1,a_vec2,a_vec_HvK,a_vec_oracle))
a_max <- max(c(a_vec1,a_vec2,a_vec_HvK,a_vec_oracle))
sig_min <- min(c(sig_vec1,sig_vec2,sig_vec_HvK,sig_vec_oracle))
sig_max <- max(c(sig_vec1,sig_vec2,sig_vec_HvK,sig_vec_oracle))
lrv_min <- min(c(lrv_vec1,lrv_vec2,lrv_vec_HvK,lrv_vec_oracle))
lrv_max <- max(c(lrv_vec1,lrv_vec2,lrv_vec_HvK,lrv_vec_oracle))

if(bounds == "yes")
{ a_min <- max(c(a_min,-1))
  a_max <- min(c(a_max,1))
  sig_min <- min(c(sig_vec1,sig_vec2,sig_vec_oracle))
  sig_max <- max(c(sig_vec1,sig_vec2,sig_vec_oracle))
  lrv_min <- min(c(lrv_vec1,lrv_vec2,lrv_vec_oracle))
  lrv_max <- max(c(lrv_vec1,lrv_vec2,lrv_vec_oracle))
}

dev.new()
par(mfrow=c(3,4))
hist(a_vec_oracle,breaks=100,xlim=c(a_min,a_max),main="oracle")
hist(a_vec1,breaks=100,xlim=c(a_min,a_max),main="our method 1")
hist(a_vec2,breaks=100,xlim=c(a_min,a_max),main="our method 2")
hist(a_vec_HvK,breaks=100,xlim=c(a_min,a_max),main="HvK method")
hist(sig_vec_oracle,breaks=100,xlim=c(sig_min,sig_max),main="oracle")
hist(sig_vec1,breaks=100,xlim=c(sig_min,sig_max),main="our method 1")
hist(sig_vec2,breaks=100,xlim=c(sig_min,sig_max),main="our method 2")
hist(sig_vec_HvK,breaks=100,xlim=c(sig_min,sig_max),main="HvK method")
hist(lrv_vec_oracle,breaks=100,xlim=c(lrv_min,lrv_max),main="oracle")
hist(lrv_vec1,breaks=100,xlim=c(lrv_min,lrv_max),main="our method 1")
hist(lrv_vec2,breaks=100,xlim=c(lrv_min,lrv_max),main="our method 2")
hist(lrv_vec_HvK,breaks=100,xlim=c(lrv_min,lrv_max),main="HvK method")


