library(mvtnorm)

set.seed(11)
x = runif(100)
d = abs(outer(x,x,"-")) # compute distance matrix, d_{ij} = |x_i - x_j|
l = 1 # length scale
Sigma_SE = exp(-d^2/(2*l^2)) # squared exponential kernel
y = mvtnorm::rmvnorm(1,sigma=Sigma_SE)
plot(x,y)