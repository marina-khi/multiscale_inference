########################
#Analysis of covid data#
########################
#rm(list=ls())

library(tidyr)
library(multiscale)
library(xtable)
library(aweek)

options(xtable.floating = FALSE)
options(xtable.timestamp = "")

library(dendextend)
library(tictoc)

source("functions.R")

#Defining necessary constants
b_bar  <- 1.02
bw_abs <- 7
t_len  <- 200
n_ts   <- 6
sigma  <- 10

# functions for data simulations
lambda_fct <- function(u, c = 1000, height = 5000, position = 10) {
  return (height * exp(-(position * u - 3) ^ 2 / 2) + c)
}

r_doublepois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta - 1))
}

simulate_data <- function(n_ts, t_len, lambda_vec, sigma) {
  data <- matrix(0, ncol = n_ts, nrow = t_len)
  for(t in 1:t_len) {
    data[t, ] <- r_doublepois(n = n_ts, lambda_vec[t], sigma^2)
  }
  return(data)
}

lambda_vec_1 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10)
lambda_vec_2 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 11)

Y1 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_1, sigma = sigma)
Y2 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_2, sigma = sigma)
#Y2 <- simulate_data(n_ts = n_ts - 1, t_len = t_len, lambda_vec = lambda_vec, sigma = sigma)
#Y  <- cbind(Y1, Y2)

plot((1:t_len) / t_len, Y1,  ylim = c(0, max(Y1, Y2) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_1,  lty = "dashed")
lines((1:t_len) / t_len, Y2, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_2, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 1  ~ and  ~ 2), line = 1)

lambda_vec_3 <- (1:t_len) * 20 + 1000
lambda_vec_4 <- lambda_vec_3 * 1.3

Y3 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_3, sigma = sigma)
Y4 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_4, sigma = sigma)

plot((1:t_len) / t_len, Y3,  ylim = c(0, max(Y3, Y4) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_3,  lty = "dashed")
lines((1:t_len) / t_len, Y4, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_4, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 3  ~ and  ~ 4), line = 1)

lambda_vec_5 <- c(rep(1000, t_len/2), rep(4000, t_len/2))
tmp <- rep(1000, t_len/3)
lambda_vec_6 <- c(tmp, rep(4000, t_len - length(tmp)))

Y5 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_5, sigma = sigma)
Y6 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_6, sigma = sigma)

plot((1:t_len) / t_len, Y5,  ylim = c(0, max(Y5, Y6) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_5,  lty = "dashed")
lines((1:t_len) / t_len, Y6, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_6, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 5  ~ and  ~ 6), line = 1)

Y <- cbind(Y1, Y2, Y3, Y4, Y5, Y6)

b_grid <- seq(1, b_bar, by = 0.01)

m_hat <- function(vect_u, b, data_p, grid_p, bw){
  t_len <- length(data_p)
  m_hat_vec <- c()
  for (u in vect_u){
    result = sum((abs((grid_p - u * b) / bw) <= 1) * data_p)
    #norm = sum((abs((grid_p - u * b) / bw) <= 1))
    norm = min(floor((u * b + bw) * t_len), t_len) - max(ceiling((u * b - bw) * t_len), 1) + 1
    m_hat_vec <- c(m_hat_vec, result/norm)
  }
  return(m_hat_vec)
}

grid_points <- seq(1/t_len, 1, by = 1/t_len)
integral_points <- seq(1/t_len, 1, by = 0.01/t_len)

m_hat_vec <- m_hat(grid_points, b = 1, Y[, 1], grid_points, bw = bw_abs/t_len)
plot(grid_points, Y[, 1], type = "l")
lines(grid_points, m_hat_vec, col = "red")

m_hat_vec <- m_hat(grid_points, b = 1, Y[, 3], grid_points, bw = bw_abs/t_len)
plot(grid_points, Y[, 3], type = "l")
lines(grid_points, m_hat_vec, col = "blue")

m_hat_vec <- m_hat(grid_points, b = 1, Y[, 5], grid_points, bw = bw_abs/t_len)
plot(grid_points, Y[, 5], type = "l")
lines(grid_points, m_hat_vec, col = "green")


integrand <- function(vect_u, b, data_points_i, data_points_j, 
                      norm_b, norm, grid_points, bw) {
  tmp <- m_hat(vect_u, b, data_points_i, grid_points, bw)/(norm_b * 1/b) - m_hat(vect_u, b = 1, data_points_j, grid_points, bw) / (norm * 1/b)
  return(tmp^2)
}

Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

tic("f-c version")
for (b in b_grid){
  norm_b <- c()
  norm   <- c()
  for (k in 1:n_ts){
    norm_b <- c(norm_b, integrate(m_hat, lower = 0, upper = 1/b, b = b,
                                  data_p = Y[, k], grid_p = grid_points,
                                  bw = bw_abs/t_len, subdivisions=2000)$value)
    norm <- c(norm, integrate(m_hat, lower = 0, upper = 1/b, b = 1,
                              data_p = Y[, k], grid_p = grid_points, 
                              bw = bw_abs/t_len, subdivisions=2000)$value)
  }
  for (i in 1:(n_ts - 1)){
    for (j in (i + 1):n_ts){
      delta_ij <- 1/b * integrate(integrand, lower = 0, upper = 1/b, b = b,
                                  data_points_i = Y[, i],
                                  data_points_j = Y[, j],
                                  norm_b = norm_b[i], norm = norm[j],
                                  grid_points = grid_points, bw = bw_abs/t_len,
                                  subdivisions=2000)$value
      delta_ji <- 1/b * integrate(integrand, lower = 0, upper = 1/b, b = b,
                                  data_points_i = Y[, j],
                                  data_points_j = Y[, i],
                                  norm_b = norm_b[j], norm = norm[i],
                                  grid_points = grid_points, bw = bw_abs/t_len,
                                  subdivisions=2000)$value
      if (b == 1) {
        Delta_hat[i, j]   <- min(delta_ij, delta_ji)
      } else {
        if (min(delta_ij, delta_ji) < Delta_hat[i, j]) {
          Delta_hat[i, j] <- min(delta_ij, delta_ji)
        }
      }
      Delta_hat[j, i] <- Delta_hat[i, j]
    }
  }
  cat("b = ", b, ": done. \n")
}

toc()

#colnames(Delta_hat) <- countries
#rownames(Delta_hat) <- countries

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
plot(res)