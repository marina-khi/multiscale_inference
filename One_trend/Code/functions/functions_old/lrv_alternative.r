
lrv_hat <- function(data, p, q, r.bar)
 
{  T <- length(data)
 

   # Step 1: compute pilot estimator of AR parameters

   y.diff  <- data[(q+1):T] - data[1:(T-q)]
   acf.hat <- acf(y.diff, lag.max=p, type="covariance", plot=FALSE, demean = FALSE)
   acf.hat <- as.vector(acf.hat$acf)

   gamma.hat <- acf.hat[2:(p+1)]
   Gamma.hat <- matrix(0,ncol=p,nrow=p)
   if(p == 1)
     Gamma.hat <- acf.hat[1]
   if(p > 1)
   { for(j in 1:(p-1))
        Gamma.hat[j,] <- c(acf.hat[j:1],acf.hat[2:(p-j+1)])
     Gamma.hat[p,] <- acf.hat[p:1]
   }

   if(p == 1) 
     a.tilde <- gamma.hat/Gamma.hat
   if(p > 1)
     a.tilde <- as.vector(solve(Gamma.hat) %*% gamma.hat) 
   

   # Step 2: compute pilot estimator of residual variance

   y.diff <- data[2:T] - data[1:(T-1)]
   T.diff <- length(y.diff)
   x.diff <- matrix(0,nrow=T.diff-p,ncol=p)
   for(i in 1:p)
      x.diff[,i] <- y.diff[(p+1-i):(T.diff-i)]
   y.diff <- y.diff[(p+1):T.diff]      

   resid <- y.diff - as.vector(x.diff %*% a.tilde)
   nu2.tilde <- sum(resid^2)/(2*T)


   # Step 3: Compute coefficients of MA expansion (c_0, c_1, ...)

   c.hat <- rep(NA,r.bar+1)
   c.hat[1] <- 1
   a.tilde.long <- c(a.tilde,rep(0,r.bar+1))

   for(j in 2:(r.bar+1))
      c.hat[j] <- sum(c.hat[(j-1):1] * a.tilde.long[1:(j-1)])
 

   # Step 4: Compute updated estimator of AR parameters

   for(r in 1:r.bar)

   {  if(r <= (p+1))
        d.hat <- - nu2.tilde * c(c.hat[(r-2):1],rep(0,(p-(r-2)))) 
      if(r > (p+1))
        d.hat <- - nu2.tilde * c.hat[(r-2):(r-p-1)]
 
      y.diff  <- data[(r+1):T] - data[1:(T-r)]
      acf.hat <- acf(y.diff, lag.max=p, type="covariance", plot=FALSE, demean = FALSE)
      acf.hat <- as.vector(acf.hat$acf)

      gamma.hat <- acf.hat[2:(p+1)]
      Gamma.hat <- matrix(0,ncol=p,nrow=p)
      if(p == 1)
        Gamma.hat <- acf.hat[1]
      if(p > 1)
      { for(j in 1:(p-1))
           Gamma.hat[j,] <- c(acf.hat[j:1],acf.hat[2:(p-j+1)])
        Gamma.hat[p,] <- acf.hat[p:1]
      }

      
   if(p == 1) 
     a.tilde <- gamma.hat/Gamma.hat
   if(p > 1)
     a.tilde <- as.vector(solve(Gamma.hat) %*% gamma.hat) 
   
      
   }
   

}

