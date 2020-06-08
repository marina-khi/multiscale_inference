
# NW estimator

epan <- function(x)
{  
   return(0.75*(1-x^2)*((sign(1-x^2)+1)/2))
}

rect <- function(x)
{
   return(0.5*((sign(1-x^2)+1)/2))
}   

nw <- function(X,Y,bw,TT) 
{  
   res <- rep(0,TT)
   for(jj in 1:TT)
   {  rh <- sum(rect((X-jj/TT)/bw) * Y) 
      fh <- sum(rect((X-jj/TT)/bw))
      res[jj] <- rh/fh
   }
   return(res)
}


# scaling factor for error variance

sigma.square <- function(Y)
{  
   Y <- as.matrix(Y)
   nn <- dim(Y)[2]
   TT <- dim(Y)[1]
   var.vec <- rep(NA,nn)
   for(ii in 1:nn)
   {  y <- Y[,ii]
      y.diff <- diff(y)
      var.vec[ii] <- sum(y.diff^2)/(2*sum(y))
   }   
   var.mean <- mean(var.vec)
   return(list(homo=var.mean, hetero=var.vec))
}   

# scaling <- function(Y, bw=0.05)
# {  Y <- as.matrix(Y)
#    nn <- dim(Y)[2]
#    TT <- dim(Y)[1]
#    X <- (1:TT)/TT
#    const <- rep(0,nn)
#    for(i in 1:nn)
#    {  lambda.fct <- nw(X,Y[,i],bw,TT)
#       resid <- Y[,i] - lambda.fct
#       pos <- (lambda.fct > 0)
#       resid <- resid[pos]/sqrt(lambda.fct[pos])
#       const[i] <- var(resid)
#    }   
#    const <- mean(const)
#    return(const)
# }   


# family of intervals for multiscale statistic

intervals <- function(TT, hmin=7, K=4)
{  
   ints.mat <- numeric(0)
   step <- floor(hmin/2)
   for(k in 1:K)
   {  nbr <- floor(TT/hmin)
      mat.base <- matrix(c(rep(1,k*hmin),rep(0,TT-k*hmin),rep(0,step),rep(1,k*hmin),rep(0,TT-k*hmin-step)),nrow=2,byrow=TRUE)
      mat.temp <- mat.base
      if(nbr>1)
      { for(jj in 1:(nbr-1))
           mat.temp <- rbind(mat.temp, cbind(matrix(0,ncol=jj*hmin,nrow=2),mat.base[,1:(TT-jj*hmin)]))
      }  
      mat.temp <- mat.temp[rowSums(mat.temp)==k*hmin, ]
      ints.mat <- rbind(ints.mat, mat.temp)
   }
   return(ints.mat)
}   
   
# ints <- intervals(TT)
# nb.ints <- dim(ints)[1]
# ints <- as.vector(ints)
# ints[ints==0] <- NA
# ints <- matrix(ints,nrow=nb.ints)
# plot(ints[1,], type="l", ylim=c(0,dim(ints)[1]+1), ylab="", xlab="time")
# for(jj in 2:dim(ints)[1])
#    lines(jj*ints[jj,])


# multiscale statistic

statistics <- function(Y, scaling)
{  
   # Inputs:
   # Y            matrix of time series (nn time series of length TT)
   # scaling      scaling factor for noise variance
   # Outputs:
   # stat.ms      nn x nn matrix containing the value of the test statistic for each pair (i,j)
   # stat.list    list of indvidual statistics
   
   TT <- dim(Y)[1]
   kernel.mat <- intervals(TT)
   h.vec <- rowSums(kernel.mat)/TT
   stat.list <- list()
   stat.ms <- matrix(NA,ncol=nn,nrow=nn)
   pos <- 1
   for(i in 1:(nn-1))
   {  for(j in (i+1):nn)
      {  nom <- as.vector(kernel.mat %*% (Y[,i] - Y[,j]))
         denom <- as.vector(sqrt(scaling) * sqrt(kernel.mat %*% (Y[,i] + Y[,j])))
         const1 <- sqrt(log(exp(1)/h.vec)) / log(log(exp(1)^exp(1)/h.vec))
         const2 <- sqrt(2*log(1/h.vec))
         stat.list[[pos]] <- const1 * (abs(nom/denom) - const2)
         stat.ms[i,j] <- max(const1 * (abs(nom/denom) - const2))
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
   nn <- dim(Z)[2]
   kernel.mat <- intervals(TT)
   h.vec <- rowSums(kernel.mat)/TT
   stat.sim <- matrix(NA,ncol=nn,nrow=nn)
   for(i in 1:(nn-1))
   {  for(j in (i+1):nn)
      {  nom <- as.vector(kernel.mat %*% (Z[,i] - Z[,j])) / sqrt(2*TT*h.vec)
         const1 <- sqrt(log(exp(1)/h.vec)) / log(log(exp(1)^exp(1)/h.vec))
         const2 <- sqrt(2*log(1/h.vec))
         stat.sim[i,j] <- max(const1 * (abs(nom) - const2)) 
      }
   }
   return(stat.sim)
}

critical.value <- function(nn, TT, level, Nsim=1000)
{  
   mstat.vec <- rep(NA, Nsim)
   for(isim in 1:Nsim)
   {  Z <- matrix(rnorm(nn*TT),ncol=nn,nrow=TT)
      mstat.vec[isim] <- max(statistics.sim(Z), na.rm=TRUE)
   }
   return(as.vector(quantile(mstat.vec, probs=1-level)))  
}



