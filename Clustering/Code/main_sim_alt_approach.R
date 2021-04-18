rm(list=ls())

library(tidyr)
library(aweek)
library(dendextend)

#source("functions.R")

#Defining necessary constants
bw_abs <- 3.5
t_len  <- 200
n_ts   <- 8
sigma  <- 15

# functions for data simulations
lambda_fct <- function(u, c = 1000, height = 5000, position = 10, a = 3) {
  return (height * exp(-(position * u - a) ^ 2 / 2) + c)
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

lambda_vec_1 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10, a = 2)
lambda_vec_2 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 11, a = 2)

Y1 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_1, sigma = sigma)
Y2 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_2, sigma = sigma)

plot((1:t_len) / t_len, Y1,  ylim = c(0, max(Y1, Y2) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_1,  lty = "dashed")
lines((1:t_len) / t_len, Y2, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_2, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 1  ~ and  ~ 2), line = 1)

lambda_vec_3 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 2000, position = 10, a = 8)
lambda_vec_4 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 2000, position = 11, a = 8)

Y3 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_3, sigma = sigma)
Y4 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_4, sigma = sigma)

plot((1:t_len) / t_len, Y3,  ylim = c(0, max(Y3, Y4) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_3,  lty = "dashed")
lines((1:t_len) / t_len, Y4, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_4, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 3  ~ and  ~ 4), line = 1)

lambda_vec_5 <- c(lambda_vec_1[1:(t_len / 2)], lambda_vec_3[(t_len / 2 + 1):t_len])
lambda_vec_6 <- c(lambda_vec_2[1:(t_len / 2)], lambda_vec_4[(t_len / 2 + 1):t_len])

Y5 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_5, sigma = sigma)
Y6 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_6, sigma = sigma)

plot((1:t_len) / t_len, Y5,  ylim = c(0, max(Y5, Y6) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_5,  lty = "dashed")
lines((1:t_len) / t_len, Y6, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_6, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 5  ~ and  ~ 6), line = 1)

lambda_vec_1_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 2000, position = 10, a = 2)
lambda_vec_2_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 2000, position = 11, a = 2)
lambda_vec_3_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10, a = 8)
lambda_vec_4_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 11, a = 8)

lambda_vec_7 <- c(lambda_vec_1_tmp[1:(t_len / 2)], lambda_vec_3_tmp[(t_len / 2 + 1):t_len])
lambda_vec_8 <- c(lambda_vec_2_tmp[1:(t_len / 2)], lambda_vec_4_tmp[(t_len / 2 + 1):t_len])

Y7 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_7, sigma = sigma)
Y8 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_8, sigma = sigma)

plot((1:t_len) / t_len, Y7,  ylim = c(0, max(Y7, Y8) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_7,  lty = "dashed")
lines((1:t_len) / t_len, Y8, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_8, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 7  ~ and  ~ 8), line = 1)

Y <- cbind(Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8)
colnames(Y) <- c("1", "2", "3", "4", "5", "6", "7", "8")

m_hat <- function(vect_u, data_p, grid_p, bw){
  m_hat_vec <- c()
  for (u in vect_u){
    result <- sum((((u - grid_p) / bw <= 1) & ((u - grid_p) / bw >= -1)) * data_p)
    norm   <- sum((((u - grid_p) / bw <= 1) & ((u - grid_p) / bw >= -1)))
    if (norm == 0){
      m_hat_vec <- c(m_hat_vec, 0)
    } else {
      m_hat_vec <- c(m_hat_vec, result/norm)
    }
  }
  return(m_hat_vec)
}

grid_points <- (1:t_len)/sqrt(t_len)

#Step 2
norm   <- c()
a_vec  <- c()
b_vec  <- c()
c_vec  <- c()
norm_p <- c()
for (k in 1:n_ts) {
  norm <- c(norm, integrate(m_hat, lower = - Inf, upper = Inf,
                            data_p = Y[, k], grid_p = grid_points,
                            bw = bw_abs/sqrt(t_len), subdivisions = 2000)$value)
  
  integrand1 <- function(x) {x * (m_hat(x, data_p = Y[, k],
                                        grid_p = grid_points,
                                        bw = bw_abs/sqrt(t_len)) / norm[k])}
  a_vec      <- c(a_vec, integrate(integrand1, lower = - Inf, upper = Inf,
                                   subdivisions = 2000)$value)
  
  integrand2 <- function(x) {x * x * (m_hat(x, data_p = Y[, k],
                                            grid_p = grid_points,
                                            bw = bw_abs/sqrt(t_len)) / norm[k])}
  tmp        <- integrate(integrand2, lower = - Inf, upper = Inf,
                          subdivisions = 2000)$value
  b_vec      <- c(b_vec, sqrt(tmp - a_vec[k]^2))
  c_vec      <- c(c_vec, norm[k] / b_vec[k])
  integrand3 <- function(x) {m_hat(a_vec[k] + b_vec[k] * x,
                                   data_p = Y[, k], grid_p = grid_points,
                                   bw = bw_abs/sqrt(t_len)) / c_vec[k]}
  norm_p     <- c(norm_p, integrate(integrand3, lower = - Inf, upper = Inf,
                                    subdivisions = 2000)$value)
}

#Matrix with the distances: Step 3
Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

for (i in 1:(n_ts - 1)){
  p_i_star <- function(x) {(m_hat(a_vec[i] + b_vec[i] * x, data_p = Y[, i],
                                  grid_p = grid_points,
                                  bw = bw_abs/sqrt(t_len)) / c_vec[i]) / norm_p[i]}
  for (j in (i + 1):n_ts){
    p_j_star <- function(x) {(m_hat(a_vec[j] + b_vec[j] * x, data_p = Y[, j],
                                   grid_p = grid_points,
                                   bw = bw_abs/sqrt(t_len)) / c_vec[j]) / norm_p[j]}
    integrand <- function(x) {(sqrt(p_i_star(x)) - sqrt(p_j_star(x)))^2}
    tmp <- integrate(integrand, lower = -Inf, upper = Inf,
                     subdivisions=2000)$value
    Delta_hat[i, j] <- tmp
    Delta_hat[j, i] <- tmp
    cat("i = ", i, ", j = ", j, " - success\n")
  }
}


colnames(Delta_hat) <- c("1", "2", "3", "4", "5", "6", "7", "8")
rownames(Delta_hat) <- c("1", "2", "3", "4", "5", "6", "7", "8")

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
plot(res)

subgroups <- cutree(res, 3)

for (cl in 1:3){
  countries_cluster <- colnames(Delta_hat)[subgroups == cl]
  if (length(countries_cluster) == 1){
    m_hat_vec <- m_hat(grid_points, b = 1, Y[, countries_cluster],
                       grid_points, bw = bw_abs/t_len)
    plot((1:t_len) / t_len, m_hat_vec,
         ylim = c(0, max(m_hat_vec) + 100), xlab="u",
         ylab = "", mgp = c(2,0.5,0), type = "l")
    title(main = paste("Representative of cluster", cl), line = 1)
  } else {
    b_res_cl <- b_res[subgroups == cl, subgroups == cl]
    colnames(b_res_cl) <- countries_cluster
    rownames(b_res_cl) <- countries_cluster
    inds               <- arrayInd(which.min(b_res_cl), dim(b_res_cl))
    repr_country       <- rownames(b_res_cl)[inds[,1]]
    repr_b             <- min(b_res_cl, na.rm = TRUE)
    m_hat_vec <- m_hat(grid_points, b = repr_b, Y[, repr_country],
                       grid_points, bw = bw_abs/t_len)
    tmp <- rep(NA, 200 - length(m_hat_vec))
    plot((1:t_len) / t_len, c(m_hat_vec, tmp),
         ylim = c(0, max(m_hat_vec) + 100), xlab="u",
         ylab = "m_hat(b * u)", mgp = c(2,0.5,0), type = "l")
    countries_cluster <- countries_cluster[countries_cluster != repr_country]
    for (country in countries_cluster){
      repr_b <- min(b_res_cl[country, ], na.rm = TRUE)
      m_hat_vec_1 <- m_hat(grid_points, b = repr_b, Y[, country],
                           grid_points, bw = bw_abs/t_len)
      tmp <- rep(NA, 200 - length(m_hat_vec_1))
      lines((1:t_len) / t_len, c(m_hat_vec_1, tmp), col = "red")
    }
    title(main = paste("Representative of cluster", cl), line = 1)
  }
}
