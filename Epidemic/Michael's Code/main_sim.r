rm(list=ls())
source("functions.r")
source("data_read.r")
library("tictoc")

Rsim  <- 5000  # number of simulation runs
alpha <- 0.05  # significance level
nn    <- 5     # number of time series
TT    <- 100   # time series length
sigma <- 15    # overdispersion parameter


# functions to simulate data

lambda.fct <- function(u){5000*exp(-(8*u-3)^2/2)+50}
lambda.vec <- lambda.fct((1:TT)/TT)

r.doublepois <- function(n, mu, theta) {
   rnbinom(n = n, mu = mu, size = mu/(theta-1))
}

simulate.data <- function(nn, TT)
{   
   data <- matrix(0,ncol=nn,nrow=TT)
   for(tt in 1:TT)
      data[tt,] <- r.doublepois(n=nn, lambda.vec[tt], sigma^2)   
   if(Rsim == 1)
   { plot(covid_mat[,1],type="l",ylim=c(0,max(max(covid_mat[,1]),max(data[,1]))))   
     lines(data[,1],col="blue")   
   }   
   return(data)
}


# compute critical value

crit.val <- critical.value(nn, TT, alpha)


# carry out multiscale test
tic("Simulations")
test.res <- rep(NA,Rsim)
for(rsim in 1:Rsim)
{
   Y <- simulate.data(nn=nn, TT=TT)
   scaling <- sigma.square(Y)$homo
   res <- statistics(Y, scaling)
   test.stat <- max(res$ms, na.rm=TRUE)
   test.res[rsim] <- as.numeric(test.stat > crit.val)
   #print(rsim)
}
toc(log = TRUE)
print(paste("Empirical size: ",  sum(test.res)/Rsim, sep=""))
