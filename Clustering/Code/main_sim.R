#####################
#Simulation analysis#
#####################
rm(list=ls())

library(tidyr)
library(aweek)
library(dendextend)
library(Rcpp)

##### Defining necessary constants
b_bar  <- 2
b_grid <- seq(1, b_bar, by = 0.1) # The grid for different b
bw_abs <- 3.5
t_len  <- 200
grid_points <- seq(1/t_len, 1, by = 1/t_len)
n_ts   <- 10
sigma  <- 15

#### Auxiliary functions
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

m_hat <- function(vect_u, b, data_p, grid_p, bw){
  m_hat_vec <- c()
  for (u in vect_u){
    result <- sum((((u * b - grid_p) / bw <= 1) & ((u * b - grid_p) / bw >= -1)) * data_p)
    norm   <- sum((((u * b- grid_p) / bw <= 1) & ((u * b - grid_p) / bw >= -1)))
    if (norm == 0){
      m_hat_vec <- c(m_hat_vec, 0)
    } else {
      m_hat_vec <- c(m_hat_vec, result/norm)
    }
  }
  return(m_hat_vec)
}

##### Simulating the data

lambda_vec_1 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10, a = 2)
lambda_vec_2 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 11, a = 2)

Y1 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_1, sigma = sigma)
Y2 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_2, sigma = sigma)

pdf("plots/simulations/plots_1_and_2.pdf", width = 6, height = 4, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(2, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot((1:t_len) / t_len, Y1,  ylim = c(0, max(Y1, Y2) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_1,  lty = "dashed")
lines((1:t_len) / t_len, Y2, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_2, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 1  ~ and  ~ 2), line = 1)
dev.off()

lambda_vec_3 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10, a = 8)
lambda_vec_4 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 11, a = 8)

Y3 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_3, sigma = sigma)
Y4 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_4, sigma = sigma)

pdf("plots/simulations/plots_3_and_4.pdf", width = 6, height = 4, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(2, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot((1:t_len) / t_len, Y3,  ylim = c(0, max(Y3, Y4) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_3,  lty = "dashed")
lines((1:t_len) / t_len, Y4, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_4, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 3  ~ and  ~ 4), line = 1)
dev.off()

lambda_vec_3_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 2000, position = 10, a = 8)
lambda_vec_4_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 2000, position = 11, a = 8)

lambda_vec_5 <- c(lambda_vec_1[1:(t_len / 2)], lambda_vec_3_tmp[(t_len / 2 + 1):t_len])
lambda_vec_6 <- c(lambda_vec_2[1:(t_len / 2)], lambda_vec_4_tmp[(t_len / 2 + 1):t_len])

Y5 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_5, sigma = sigma)
Y6 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_6, sigma = sigma)

pdf("plots/simulations/plots_5_and_6.pdf", width = 6, height = 4, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(2, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot((1:t_len) / t_len, Y5,  ylim = c(0, max(Y5, Y6) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_5,  lty = "dashed")
lines((1:t_len) / t_len, Y6, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_6, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 5  ~ and  ~ 6), line = 1)
dev.off()

lambda_vec_1_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 2000, position = 10, a = 2)
lambda_vec_2_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 2000, position = 11, a = 2)
lambda_vec_3_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10, a = 8)
lambda_vec_4_tmp <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 11, a = 8)

lambda_vec_7 <- c(lambda_vec_1_tmp[1:(t_len / 2)], lambda_vec_3_tmp[(t_len / 2 + 1):t_len])
lambda_vec_8 <- c(lambda_vec_2_tmp[1:(t_len / 2)], lambda_vec_4_tmp[(t_len / 2 + 1):t_len])

Y7 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_7, sigma = sigma)
Y8 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_8, sigma = sigma)

pdf("plots/simulations/plots_7_and_8.pdf", width = 6, height = 4, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(2, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot((1:t_len) / t_len, Y7,  ylim = c(0, max(Y7, Y8) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_7,  lty = "dashed")
lines((1:t_len) / t_len, Y8, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_8, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 7  ~ and  ~ 8), line = 1)
dev.off()

lambda_vec_9 <- c(lambda_vec_1[1:(t_len * 2 / 10)], rep(lambda_vec_1[t_len * 2 / 10], t_len / 2),
                  lambda_vec_1[(t_len * 2 / 10 + 1):(t_len / 2)])

Y9 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_9, sigma = sigma)
Y10 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_9, sigma = sigma)

pdf("plots/simulations/plots_9_and_10.pdf", width = 6, height = 4, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(2, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot((1:t_len) / t_len, Y9,  ylim = c(0, max(Y9, Y10) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
lines((1:t_len) / t_len, lambda_vec_9,  lty = "dashed")
lines((1:t_len) / t_len, Y10, type = "l", col = "blue")
lines((1:t_len) / t_len, lambda_vec_9, col = "blue", lty = "dashed")
title(main = expression(Plot ~ of ~ the ~ time ~ series ~ 9  ~ and  ~ 10), line = 1)
dev.off()

Y <- cbind(Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y8, Y9, Y10)
colnames(Y) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

##### Old approach (with multiple b)

Rcpp::sourceCpp("integration.cpp")

Delta_hat_tmp <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
b_res         <- matrix(data = rep(NA, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

for (b in b_grid){
  norm_b <- c()
  norm   <- c()
  for (k in 1:n_ts){
    norm_b <- c(norm_b, integrate1_cpp(b = b, data_points = Y[, k],
                                       grid_points = grid_points,
                                       bw = bw_abs/t_len, subdiv = 2000)$res)
    norm <- c(norm, integrate1_cpp(b = 1.0, data_points = Y[, k],
                                   grid_points = grid_points,
                                   bw = bw_abs/t_len, subdiv = 2000)$res)
  }
  for (i in 1:(n_ts - 1)){
    for (j in (i + 1):n_ts){
      delta_ij <- 1/b * integrate2_cpp(b = b, data_points_1 = Y[, i],
                                       data_points_2 = Y[, j],
                                       norm_1 = norm_b[i], norm_2 = norm[j],
                                       grid_points = grid_points, bw = bw_abs/t_len,
                                       subdiv=2000)$res
      delta_ji <- 1/b * integrate2_cpp(b = b, data_points_1 = Y[, j],
                                       data_points_2 = Y[, i],
                                       norm_1 = norm_b[j], norm_2 = norm[i],
                                       grid_points = grid_points, bw = bw_abs/t_len,
                                       subdiv=2000)$res
      if (b == 1) {
        Delta_hat_tmp[i, j] <- delta_ij
        Delta_hat_tmp[j, i] <- delta_ji
        b_res[i, j] <- 1
        b_res[j, i] <- 1
      } else {
        if (delta_ij < Delta_hat_tmp[i, j]) {
          Delta_hat_tmp[i, j] <- delta_ij
          b_res[i, j] <- b
          b_res[j, i] <- 1
        } 
        if (delta_ji < Delta_hat_tmp[j, i]) {
          Delta_hat_tmp[j, i] <- delta_ji
          b_res[j, i] <- b
          b_res[i, j] <- 1          
        }
      }
    }
  }
  cat("b = ", b, " - success\n")
}

#Delta_hat_tmp was a temporary non-symmetrical matrix,
#for the distance matrix we need a symmetrical one

Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
for (i in 1:(n_ts - 1)){
  for (j in (i + 1):n_ts){
    Delta_hat[i, j] <- min(Delta_hat_tmp[i, j], Delta_hat_tmp[j, i])
    Delta_hat[j, i] <- Delta_hat[i, j]
  }
}

colnames(Delta_hat) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
rownames(Delta_hat) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

colnames(b_res) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
rownames(b_res) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)

pdf("plots/simulations/dendrogram.pdf", width = 15, height = 6, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot(res, cex = 0.8, xlab = "", ylab = "")
rect.hclust(res, k = 5, border = 2:6)
dev.off()

###########################
##### Alternative approach#
###########################
grid_points_alt <- (1:t_len)/sqrt(t_len)

#Step 2
norm   <- c()
a_vec  <- c()
b_vec  <- c()
c_vec  <- c()
norm_p <- c()
for (k in 1:n_ts) {
  norm <- c(norm, integrate(m_hat, lower = - Inf, upper = Inf,
                            b = 1, data_p = Y[, k], grid_p = grid_points_alt,
                            bw = bw_abs/sqrt(t_len), subdivisions = 2000)$value)
  
  integrand1 <- function(x) {x * (m_hat(x, b = 1, data_p = Y[, k],
                                        grid_p = grid_points_alt,
                                        bw = bw_abs/sqrt(t_len)) / norm[k])}
  a_vec      <- c(a_vec, integrate(integrand1, lower = - Inf, upper = Inf,
                                   subdivisions = 2000)$value)
  
  integrand2 <- function(x) {x * x * (m_hat(x, b = 1, data_p = Y[, k],
                                            grid_p = grid_points_alt,
                                            bw = bw_abs/sqrt(t_len)) / norm[k])}
  tmp        <- integrate(integrand2, lower = - Inf, upper = Inf,
                          subdivisions = 2000)$value
  b_vec      <- c(b_vec, sqrt(tmp - a_vec[k]^2))
  c_vec      <- c(c_vec, norm[k] / b_vec[k])
  integrand3 <- function(x) {m_hat(a_vec[k] + b_vec[k] * x, b = 1, 
                                   data_p = Y[, k], grid_p = grid_points_alt,
                                   bw = bw_abs/sqrt(t_len)) / c_vec[k]}
  norm_p     <- c(norm_p, integrate(integrand3, lower = - Inf, upper = Inf,
                                    subdivisions = 2000)$value)
}

#Matrix with the distances: Step 3
Delta_hat_alt <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

for (i in 1:(n_ts - 1)){
  p_i_star <- function(x) {(m_hat(a_vec[i] + b_vec[i] * x, b = 1, data_p = Y[, i],
                                  grid_p = grid_points_alt,
                                  bw = bw_abs/sqrt(t_len)) / c_vec[i]) / norm_p[i]}
  for (j in (i + 1):n_ts){
    p_j_star <- function(x) {(m_hat(a_vec[j] + b_vec[j] * x, b = 1, data_p = Y[, j],
                                    grid_p = grid_points_alt,
                                    bw = bw_abs/sqrt(t_len)) / c_vec[j]) / norm_p[j]}
    integrand <- function(x) {(sqrt(p_i_star(x)) - sqrt(p_j_star(x)))^2}
    tmp <- integrate(integrand, lower = -Inf, upper = Inf,
                     subdivisions=2000)$value
    Delta_hat_alt[i, j] <- tmp
    Delta_hat_alt[j, i] <- tmp
    cat("i = ", i, ", j = ", j, " - success\n")
  }
}


colnames(Delta_hat_alt) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")
rownames(Delta_hat_alt) <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")

delta_dist_alt <- as.dist(Delta_hat_alt)
res_alt        <- hclust(delta_dist_alt)

pdf("plots/simulations/dendrogram_alternative.pdf", width = 15, height = 6, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot(res_alt, cex = 0.8, xlab = "", ylab = "")
rect.hclust(res_alt, k = 5, border = 2:6)
dev.off()
