
# NW estimator

epan <- function(x)
{  
   return(0.75*(1-x^2)*((sign(1-x^2)+1)/2))
}

nw <- function(X,Y,bw,N) 
{  
   res <- rep(0,N)
   for(j in 1:N)
   {  rh <- sum(epan((X-j/N)/bw) * Y) 
      fh <- sum(epan((X-j/N)/bw))
      res[j] <- rh/fh
   }
   return(res)
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
   

# intervals for multiscale statistic

intervals <- function(TT, hmin=7, K=4)
{  
   ivals.mat <- numeric(0)
   for(k in 1:K)
   {  h <- k*hmin  
      vec <- rep(c(rep(1,h),rep(0,TT-h),rep(0,hmin)), ceiling(TT/hmin))
      vec <- vec[1:(TT * (floor(TT/hmin) - (k-1)))]              
      mat <- matrix(vec,ncol=TT,nrow=floor(TT/hmin)-(k-1),byrow=TRUE)
      ivals.mat <- rbind(ivals.mat,mat)  
   }
   return(ivals.mat)
}


# multiscale statistic

statistics <- function(Y, correction)
{  
   # Inputs:
   # Y            matrix of time series (nn time series of length TT)
   # correction   correction factor for noise variance
   # Outputs:
   # stat.ms      nn x nn matrix containing the value of the test statistic for each pair (i,j)
   # stat.list    list of indvidual statistics
   
   TT <- dim(Y)[1]
   kernel.mat <- intervals(TT)
   h.vec <- rowSums(kernel.mat)/TT
   stat.list <- list()
   stat.ms <- matrix(0,ncol=nn,nrow=nn)
   pos <- 1
   for(i in 1:(nn-1))
   {  for(j in (i+1):nn)
      {  nom <- as.vector(kernel.mat %*% (Y[,i] - Y[,j]))
         denom <- as.vector(sqrt(correction) * sqrt(kernel.mat %*% (Y[,i] + Y[,j])))
         stat.list[[pos]] <- abs(nom/denom) - sqrt(2*log(1/h.vec))
         stat.ms[i,j] <- max(abs(nom/denom) - sqrt(2*log(1/h.vec)))
         pos <- pos+1
      }
   }
   return(list(ms=stat.ms, stats=stat.list))
}


# critical value

statistics.sim <- function(Z)
{  
   # Inputs:
   # Z            TT x nn matrix of simulated N(0,1) variables
   # Outputs:
   # stat.sim     nn x nn matrix containing the value of the test statistic for each pair (i,j)
   
   TT <- dim(Z)[1]
   kernel.mat <- intervals(TT)
   h.vec <- rowSums(kernel.mat)/TT
   stat.sim <- matrix(0,ncol=nn,nrow=nn)
   for(i in 1:(nn-1))
   {  for(j in (i+1):nn)
      {  nom <- as.vector(kernel.mat %*% (Z[,i] - Z[,j])) / sqrt(2*TT*h.vec)
         stat.sim[i,j] <- max(abs(nom) - sqrt(2*log(1/h.vec))) 
      }
   }
   return(stat.sim)
}


critical.value <- function(nn, TT, level, Nsim=1000)
{  
   mstat.vec <- rep(NA, Nsim)
   for(isim in 1:Nsim)
   {  Z <- matrix(rnorm(nn*TT),ncol=nn,nrow=TT)
      mstat.vec[isim] <- max(statistics.sim(Z))
   }
   return(as.vector(quantile(mstat.vec, probs=1-level)))  
}



