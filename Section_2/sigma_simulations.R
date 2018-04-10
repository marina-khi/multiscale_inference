dyn.load("estimating_sigma.dll")
source("estimating_sigma.R")
source("auxiliary.R")
library("robts")
library("microbenchmark")

L1<-10
L2<-20

a <- arima.sim(model = list(ar = 0.5), innov = rnorm(1000, 0, 1), n = 1000)
microbenchmark(estimating_sigma_for_AR1(a, L1, L2))
microbenchmark(estimating_sigma_for_AR1_old_version(a, L1, L2))

#sigmahat_for_iid <- replicate(1000, {
#  x <- rnorm(400, 0, 1)
#  estimating_sigma_for_AR1(x, L1, L2)
#})

# sigmahat_for_AR_075 <- replicate(1000, {
#   a <- arima.sim(model = list(ar = 0.75), innov = rnorm(1000, 0, 1), n = 1000)
#   estimating_sigma_for_AR1(a, L1, L2)
# }) #True sigma is 16
# 
# sigmahat_for_AR_075_with_param <- replicate(1000, {
#   a <- arima.sim(model = list(ar = 0.75), innov = rnorm(1000, 0, 1), n = 1000)
#   estimating_sigma_for_AR1_with_param(a, L1, L2)
# }) #True sigma is 16
# 
# sigmahat_for_AR_05 <- replicate(1000, {
#   b <- arima.sim(model = list(ar = 0.5), innov = rnorm(1000, 0, 1), n = 1000)
#   estimating_sigma_for_AR1(b, L1, L2)
# }) #True sigma is 4
# 
# sigmahat_for_AR_05_with_param <- replicate(1000, {
#   b <- arima.sim(model = list(ar = 0.5), innov = rnorm(1000, 0, 1), n = 1000)
#   estimating_sigma_for_AR1_with_param(b, L1, L2)
# }) #True sigma is 4
# 
# 
# hist(sigmahat_for_AR_05, breaks=seq(1.0,7.0,by=0.3))
# hist(sigmahat_for_AR_05_with_param, breaks=seq(1.0,7.0,by=0.3))
# hist(sigmahat_for_AR_075, breaks=seq(0,32, by = 2))
# hist(sigmahat_for_AR_075_with_param, breaks=seq(0,32, by = 2))

#estimating_sigma_for_AR1(yearly_temp, L1, L2)
#asymvar.acf(yearly_temp)$lrv