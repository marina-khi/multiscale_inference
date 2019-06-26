rm(list=ls())

source("functions/long_run_variance.r")


# Parameters


q     <- 25  # tuning parameter for first-step estimator
r.bar <- 10  # tuning parameter for second-step estimator

p          <- 1         # AR order
sigma_eta  <- 1         # standard deviation of the innovation term in the AR model
sim.design <- "linear"  # trend specification: "constant", "blocks", "spike", ...
a1         <- 0.95      # AR parameter       
slope.fac  <- 10        # slope of linear trend = slope_fac * sqrt(sigma_eta^2/(1-a1^2))

T     <- 500   # sample sizes 
Nsim  <- 1000  # number of simulation runs


# Simulations

a.hat <- a.oracle <- rep(0,Nsim)
lrv.hat <- lrv.oracle <- rep(0,Nsim)

for(nb in 1:Nsim)

{  set.seed(nb)

   # simulate data
   source("simulations/sim.r")  

   # compute estimator
   AR.struc    <- AR_lrv(data=data, q=25, r.bar=10, p=1)    
   a.hat[nb]   <- AR.struc$ahat 
   lrv.hat[nb] <- AR.struc$lrv

   # compute oracle estimators
   res.oracle    <- lm(eps[2:T] ~ eps[1:(T-1)] - 1)
   a.oracle[nb]  <- as.vector(res.oracle$coef)
   sig.oracle     <- mean(res.oracle$residuals^2)
   lrv.oracle[nb] <- sig.oracle/(1-sum(a.oracle[nb]))^2 

   print(nb)
} 


dev.new()
par(mfrow=c(2,2))
hist(a.hat)
abline(v=a1,col="red",lwd=2)
hist(a.oracle)
abline(v=a1,col="red",lwd=2)
hist(lrv.hat)
abline(v=sigma^2,col="red",lwd=2)
hist(lrv.oracle)
abline(v=sigma^2,col="red",lwd=2)