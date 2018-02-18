source("functions.R")
library(microbenchmark)

set.seed(1)
sigma <- 2
T<-100
y_data <- 1 + rnorm(T, 0, sigma)
sigmahat <- sigma
alpha <-0.05

#Creating g_t_set like in Section 2.1
u <- seq(4/T, 1, length.out = T/4)
h <- seq(3/T, 1/4+3/T, length.out = T/20)

g_t_set_temp <- expand.grid(u = u, h = h) #Creating a dataframe with all possible combination of u and h
g_t_set_temp$values <-numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero

#Subsetting u and h such that [u-h, u+h] lies in [0,1]
g_t_set <- subset(g_t_set_temp, u - h >= 0 & u + h <= 1, select = c(u, h, values))
g_t_set$lambda <- lambda(g_t_set[['h']])

psihat_statistic_value<- psihat_statistic(y_data, g_t_set, epanechnikov_kernel, sigmahat)

microbenchmark(psihat_statistic(y_data, g_t_set, epanechnikov_kernel, sigmahat))

gaussian_statistic_distribution <- replicate(1000, {
  z = rnorm(T, 0, 1)
  z_temp = sigma * z
  psihat_statistic(z_temp, g_t_set, epanechnikov_kernel, sigma)
})
save(gaussian_statistic_distribution, file = 'distribution.RData')
gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)

if (psihat_statistic_value>=gaussian_quantile)
{
  cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
      "Gaussian quantile value = ", gaussian_quantile)
} else {
  cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
      "Gaussian quantile value = ", gaussian_quantile)
}
g_t_set_rejected <- subset(g_t_set, values >= gaussian_quantile, select = c(u, h, values))