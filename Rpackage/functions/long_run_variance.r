
emp_acf <- function(data,ell,p)

{   # computes autocovariances at lags 0 to p 
    # for the ell-th differences of data

    y.diff   <- diff(data,ell)
    len      <- length(y.diff)
    autocovs <- rep(0,p+1)
    for(k in 1:(p+1))
       autocovs[k] <- sum(y.diff[k:len]*y.diff[1:(len-k+1)])/len
    return(autocovs)
}


variance_eta <- function(data,coefs,p){
    # computes variance of AR innovation terms eta   
 
    y.diff <- diff(data)
    len    <- length(y.diff)
    x.diff <- matrix(0,nrow=len-p,ncol=p)
    for(i in 1:p)
       x.diff[,i] <- y.diff[(p+1-i):(len-i)]
    y.diff <- y.diff[(p+1):len]      

    resid   <- y.diff - as.vector(x.diff %*% coefs)
    var.eta <- mean(resid^2)/2
    return(var.eta)
}


corrections <- function(coefs,var.eta,len)

{   # computes vector of correction terms for second-stage estimator
    # of AR parameters

    p <- length(coefs)  
    c.vec <- rep(0,len)
    c.vec[1] <- 1
    for(j in 2:len)
    { lags <- (j-1):max(j-p,1)
      c.vec[j] <- sum(coefs[1:length(lags)] * c.vec[lags])  
    }
    c.vec <- c.vec * var.eta
    return(c.vec)
}


AR_coef <- function(data,L1,L2,correct,p){
    # computes estimator of the AR coefficients  
    pos <- 0
    a.mat <- matrix(0,nrow=L2-L1+1,ncol=p)
    for(ell in L1:L2)
    {  pos <- pos + 1
       autocovs <- emp_acf(data,ell,p)
       cov.mat  <- matrix(0,ncol=p,nrow=p)
       for(i in 1:p)
       {  for(j in 1:p)
             cov.mat[i,j] <- autocovs[abs(i-j)+1]
       }
       if(ell >= p)
         correct.vec <- correct[(ell-1+1):(ell-p+1)]
       if(ell < p)
         correct.vec <- c(correct[(ell-1+1):1],rep(0,p-ell))  
       cov.vec <- autocovs[2:(p+1)] + correct.vec
       a.mat[pos,] <- solve(cov.mat) %*% cov.vec 
    }    
    a.hat <- colMeans(a.mat)
    a.hat <- as.vector(a.hat)
    return(a.hat)
}


AR_lrv <- function(data, q, r.bar, p){
    # computes long-run variance of error terms
    a.tilde       <- AR_coef(data=data, L1=q, L2=q, correct=rep(0,max(p,q)+1), p=p)
    sig.eta.tilde <- variance_eta(data=data, coefs=a.tilde, p=p)
    correct       <- corrections(coefs=a.tilde, var.eta=sig.eta.tilde, len=r.bar+1)
    a.hat         <- AR_coef(data=data, L1=1, L2=r.bar, correct=correct, p=p)
    sig.eta.hat   <- variance_eta(data=data, coefs=a.hat, p=p)
    lrv.hat       <- sig.eta.hat/(1-sum(a.hat))^2 

    return(list(lrv=lrv.hat, ahat=a.hat, vareta=sig.eta.hat, atilde=a.tilde))
}


AR_acf <- function(coefs, var.eta, len){
    # computes autocovariance function based on AR coefficients
    p <- length(coefs)  
    len <- max(len,50)
    c.vec <- rep(0,len)
    c.vec[1] <- 1
    for(j in 2:len)
    { lags <- (j-1):max(j-p,1)
      c.vec[j] <- sum(coefs[1:length(lags)] * c.vec[lags])  
    }    

    gamma.AR <- rep(0,len)
    gamma.AR[1] <- var.eta * sum(c.vec^2)
    for(j in 1:(len-1))
       gamma.AR[j+1] <- var.eta * sum(c.vec[1:(len-j)] * c.vec[(j+1):len])
    
    return(gamma.AR)
}


AR_coef_HvK <- function(data,L1,L2,p)

{   # computes estimator of AR coefficients from Hall & vanKeilegom (2003)

    pos <- 0
    var.vec <- rep(0,L2-L1+1)
    for(ell in L1:L2)
    {  pos <- pos + 1
       y.diff <- diff(data,lag=ell)
       var.vec[pos] <- mean(y.diff^2)
    }
    g0 <- mean(var.vec)/2

    g.vec <- rep(0,p)
    for(ell in 1:p)
    {  y.diff <- diff(data,lag=ell)
       g.vec[ell] <- g0 - mean(y.diff^2)/2
    }
    g.vec <- c(g0,g.vec)

    cov.mat <- matrix(0,ncol=p,nrow=p)
    for(i in 1:p)
    {  for(j in 1:p)
          cov.mat[i,j] <- g.vec[abs(i-j)+1]
    }
    cov.vec <- g.vec[2:(p+1)] 
  
    a.hat.HvK <- solve(cov.mat) %*% cov.vec 
    a.hat.HvK <- as.vector(a.hat.HvK)
    return(a.hat.HvK)
}

sigmahat_vec_iid <- function(data){
  #function that calculates the vector of sqrt(variance) estimators for a number of time series.
  #It uses the standard first differences method.
  N    = ncol(data)
  Tlen = nrow(data)
  sigmahat_vec <- c()
  for (i in 1:N){
    variance_i   <- sum((data[2:Tlen, i] - data[1:(Tlen - 1), i])^2)/(2 * Tlen - 2)
    sigma_hat_i  <- sqrt(variance_i)
    sigmahat_vec <- c(sigmahat_vec, sigma_hat_i)
  }
  return(sigmahat_vec)
}

sigmahat_vec_iid2 <- function(data){
  #function that calculates the vector of sqrt(variance) estimators for a number of time series.
  #It uses the standard "three neigbours" method.
  N    = ncol(data)
  Tlen = nrow(data)
  sigmahat_vec <- c()
  for (i in 1:N){
    variance_i   <- 2 * sum((data[3:Tlen, i]/2 - data[2:(Tlen - 1), i] + data[1:(Tlen - 2), i]/2)^2)/(3 * Tlen - 6)
    sigma_hat_i  <- sqrt(variance_i)
    sigmahat_vec <- c(sigmahat_vec, sigma_hat_i)
  }
  return(sigmahat_vec)
}