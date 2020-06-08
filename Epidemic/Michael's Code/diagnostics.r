rm(list=ls())
source("data_read.r")
source("functions.r")


# fit trend to data

cases <- covid_mat[,1]
TT <- length(cases)
Y <- cases
X <- (1:TT)/TT
bw <- 0.025
lambda.hat <- nw(X,Y,bw,TT)

plot(lambda.hat,type="l",ylim=c(min(Y),max(Y)))
lines(Y,lty="dashed")


# compute correction factor for variance 

resid <- Y - lambda.hat
acf(resid,type="correlation")

pos <- (lambda.hat > 0)
resid <- resid[pos]/sqrt(lambda.hat[pos])
c.var <- sqrt(var(resid))


# compare actual time series with a Poisson simulation

Y.sim <- rep(0,TT)
for(i in 1:TT)
  Y.sim[i] <- rpois(n=1, lambda=lambda.hat[i])

plot(Y,type="l",ylim=c(0,max(max(Y),max(Y.sim))))
lines(Y.sim,col="red")
lines(lambda.hat,col="blue")

