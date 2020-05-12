rm(list=ls())
source("functions.r")
source("data_read.r")


Rsim <- 1000  # number of simulation runs
alpha <- 0.05  # significance level
nn <- 2  # number of time series
TT <- dim(covid_mat)[1]  # length of time series
correction <- 1


# function to simulate data

simulate.data <- function(nn, TT, bw=0.025)
{   
   cases <- covid_mat[,1]
   X <- (1:TT)/TT
   lambda.hat <- nw(X,cases,bw,TT)
   data <- matrix(0,ncol=nn,nrow=TT)
   for(ii in 1:nn)
   {  for(tt in 1:TT)
         data[tt,ii] <- rpois(n=1, lambda=lambda.hat[tt])
   }
   return(data)
}

# compute critical value

crit.val <- critical.value(nn, TT, alpha)


# carry out multiscale test

test.res <- rep(NA,Rsim)
for(rsim in 1:Rsim)
{
   Y <- simulate.data(nn=nn, TT=TT)
   res <- statistics(Y, correction)
   test.stat <- max(res$ms)
   test.res[rsim] <- as.numeric(test.stat > crit.val)
   print(rsim)
}

print(paste("Empirical size: ",  sum(test.res)/Rsim, sep=""))
