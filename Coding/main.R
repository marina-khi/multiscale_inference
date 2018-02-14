source("functions.R")
library(microbenchmark)

set.seed(200)
sigma <- 2
T<-1000
y_data <- rnorm(T, 0, sigma)
sigmahat <- sigma

#Creating g_t_set like in Section 2.1
u <- seq(1/T, 1, length.out = T)
h <- seq(1/T, 1, length.out = T)

microbenchmark(g_t_set_temp <- expand.grid(u = u, h = h)) #Creating a dataframe with all possible combination of u and h
g_t_set_temp$values <-numeric(T^2) # Setting the values of the statistic to be zero

#Subsetting such that [u-h, u+h] \in [0,1]
microbenchmark(g_t_set <- subset(g_t_set_temp, u - h >= 0 & u + h <= 1, select = c(u, h, values)))



microbenchmark(psihat_statistic(y_data, g_t_set, epanechnikov_kernel, sigmahat))