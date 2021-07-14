###############################################
#Analysis of covid data - alternative approach#
###############################################
rm(list=ls())

library(tidyr)
library(aweek)
library(dendextend)
require(rworldmap)
source("functions.R")

#Half of the bandwidth window (in absolute terms)
bw_abs <- 7

#Loading the data
data  <- data_load(aligning = FALSE)
t_len <- data$t_len
n_ts  <- data$n_ts
covid_mat <- data$covid_mat
countries <- names(data$covid_list)

grid_points <- (1:t_len)/sqrt(t_len)

#Step 2
norm   <- c()
a_vec  <- c()
b_vec  <- c()
c_vec  <- c()
norm_p <- c()
for (k in 1:n_ts) {
  norm <- c(norm, integrate(m_hat_standard, lower = - Inf, upper = Inf,
                            data_p = covid_mat[, k], grid_p = grid_points,
                            bw = bw_abs/sqrt(t_len), subdivisions = 2000)$value)
  
  integrand1 <- function(x) {x * m_hat_standard(x, data_p = covid_mat[, k],
                                                grid_p = grid_points,
                                                bw = bw_abs/sqrt(t_len)) / norm[k]}
  a_vec      <- c(a_vec, integrate(integrand1, lower = - Inf, upper = Inf,
                                   subdivisions = 2000)$value)
  
  integrand2 <- function(x) {x * x * m_hat_standard(x, data_p = covid_mat[, k],
                                           grid_p = grid_points,
                                           bw = bw_abs/sqrt(t_len)) / norm[k]}
  tmp        <- integrate(integrand2, lower = - Inf, upper = Inf,
                          subdivisions = 2000)$value
  b_vec      <- c(b_vec, sqrt(tmp - a_vec[k]^2))
  c_vec      <- c(c_vec, norm[k] / b_vec[k])
  integrand3 <- function(x) {m_hat_standard(a_vec[k] + b_vec[k] * x,
                                            data_p = covid_mat[, k],
                                            grid_p = grid_points,
                                            bw = bw_abs/sqrt(t_len)) / c_vec[k]}
  norm_p <- c(norm_p, integrate(integrand3, lower = - Inf, upper = Inf,
                                subdivisions = 2000)$value)
}

# #Matrix with the distances: Step 3
# Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
# 
# for (i in 1:(n_ts - 1)){
#   p_i_star <- function(x) {(m_hat_standard(a_vec[i] + b_vec[i] * x, data_p = covid_mat[, i],
#                                   grid_p = grid_points,
#                                   bw = bw_abs/sqrt(t_len)) / c_vec[i]) / norm_p[i]}
#   for (j in (i + 1):n_ts){
#     p_j_star <- function(x) {(m_hat_standard(a_vec[j] + b_vec[j] * x, data_p = covid_mat[, j],
#                                     grid_p = grid_points,
#                                     bw = bw_abs/sqrt(t_len)) / c_vec[j]) / norm_p[j]}
#     integrand <- function(x) {(sqrt(p_i_star(x)) - sqrt(p_j_star(x)))^2}
#     tmp <- integrate(integrand, lower = -Inf, upper = Inf,
#                      subdivisions=3000)$value
#     Delta_hat[i, j] <- tmp
#     Delta_hat[j, i] <- tmp
#     cat("i = ", i, ", j = ", j, " - success\n")
#   }
# }
# 
# #And now the clustering itself
# colnames(Delta_hat) <- countries
# rownames(Delta_hat) <- countries
# 
# save(Delta_hat, file = "results_alt_approach_14days.RData")
load("results_alt_approach_14days.RData")

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)

g_hat <- function(x_vec, covid_mat, inds, a_vec, b_vec, c_vec, norm_p, bw_abs,
                  t_len, grid_points){
  g_hat_vec <- rep(0, length(x_vec))
  for (ind in inds){
    p_i_hat_star <- m_hat_standard(a_vec[ind] + b_vec[ind] * x_vec,
                                   covid_mat[, ind], grid_points,
                                   bw = bw_abs/sqrt(t_len)) / (c_vec[ind] * norm_p[ind])
    g_hat_vec    <- g_hat_vec + p_i_hat_star
  }
  return(g_hat_vec/length(inds))
}

max_n_cl <- 40
BIC_vec  <- c()

for (n_cl in 1:max_n_cl){
  subgroups         <- cutree(res, n_cl)
  sigma_hat_vec     <- rep(0, n_ts)
  for (cl in 1:n_cl){
    countries_cluster <- colnames(Delta_hat)[subgroups == cl]
    inds              <- match(countries_cluster, countries)
    for (ind in inds){
      x_vec <- (grid_points - a_vec[ind])/b_vec[ind]
      g_hat_normed <- c_vec[ind] * g_hat(x_vec = x_vec, covid_mat = covid_mat,
                                         inds = inds, a_vec = a_vec, b_vec = b_vec,
                                         c_vec = c_vec, norm_p = norm_p,
                                         bw_abs = bw_abs, t_len = t_len,
                                         grid_point = grid_points)
      tmp <- covid_mat[, ind] - g_hat_normed
      sigma_hat_vec[ind] <- sum(tmp^2)/t_len
    }
  }
  BIC <- t_len * sum(log(sigma_hat_vec))# - log(n_ts * t_len) * (n_cl * (n_ts + t_len) + n_ts)
  BIC_vec <- c(BIC_vec, BIC)
}


pdf("plots/BIC_alt_withoutpenalty.pdf", width=7, height=5, paper="special")
par(cex = 1, tck = -0.025)
par(mar = c(3.5, 2.5, 4, 0)) #Margins for each plot
plot(1:max_n_cl, BIC_vec, type = "l", xlab = "", ylab = "")
title(main = "COVID data", line = 2.5)
title(main = expression(BIC==T~Sigma[i==1]^n~log(sigma[i]^2~(K))), line = 1)
title(xlab = "Number of clusters, K", line = 2.5)
dev.off()

# pdf("plots/BIC_alt.pdf", width=7, height=5, paper="special")
# par(cex = 1, tck = -0.025)
# par(mar = c(3.5, 2.5, 4, 0)) #Margins for each plot
# plot(1:max_n_cl, BIC_vec, type = "l", xlab = "", ylab = "")
# title(main = "COVID data", line = 2.5)
# title(main = expression(BIC==T~Sigma[i==1]^n~log(sigma[i]^2~(K))-log(nT)(K(n+T)+n)), line = 1)
# title(xlab = "Number of clusters, K", line = 2.5)
# dev.off()

AIC_vec <- c()
BIC_vec_2  <- c()
for (n_cl in 1:max_n_cl){
  subgroups         <- cutree(res, n_cl)
  sigma_hat_vec     <- rep(0, n_ts)
  for (cl in 1:n_cl){
    countries_cluster <- colnames(Delta_hat)[subgroups == cl]
    inds              <- match(countries_cluster, countries)
    for (ind in inds){
      x_vec <- (grid_points - a_vec[ind])/b_vec[ind]
      g_hat_normed <- c_vec[ind] * g_hat(x_vec = x_vec, covid_mat = covid_mat,
                                         inds = inds, a_vec = a_vec, b_vec = b_vec,
                                         c_vec = c_vec, norm_p = norm_p,
                                         bw_abs = bw_abs, t_len = t_len,
                                         grid_point = grid_points)
      sigma_hat_vec[ind] <- sum((covid_mat[, ind] - g_hat_normed)^2)
    }
  }
  BIC <- (-1) * n_ts * t_len * log(sum(sigma_hat_vec)/(n_ts * t_len))# - n_cl * n_ts * t_len * (2 * bw_abs / t_len) * log(n_ts * t_len)
  BIC_vec_2 <- c(BIC_vec_2, BIC)
  AIC <- (-1) * n_ts * t_len * log(sum(sigma_hat_vec)/(n_ts * t_len))# - n_cl * n_ts * t_len * (2 * bw_abs / t_len)
  AIC_vec <- c(AIC_vec, AIC)
}

pdf("plots/BIC_2_alt_withoutpenalty.pdf", width=7, height=5, paper="special")
par(cex = 1, tck = -0.025)
par(mar = c(3.5, 2.5, 4, 0)) #Margins for each plot
plot(1:max_n_cl, BIC_vec_2, type = "l", xlab = "", ylab = "")
title(main = "COVID data", line = 2.5)
title(main = expression(BIC==-nT~log(tilde(sigma)^2)), line = 1)
title(xlab = "Number of clusters, K", line = 2.5)
dev.off()


pdf("plots/AIC_alt_wihtoutpenalty.pdf", width=7, height=5, paper="special")
par(cex = 1, tck = -0.025)
par(mar = c(3.5, 2.5, 4, 0)) #Margins for each plot
plot(1:max_n_cl, AIC_vec, type = "l", xlab = "", ylab = "")
title(main = "COVID data", line = 2.5)
title(main = expression(AIC==-nT~log(tilde(sigma)^2)), line = 1)
title(xlab = "Number of clusters, K", line = 2.5)
dev.off()


# pdf("plots/BIC_2_alt.pdf", width=7, height=5, paper="special")
# par(cex = 1, tck = -0.025)
# par(mar = c(3.5, 2.5, 4, 0)) #Margins for each plot
# plot(1:max_n_cl, BIC_vec_2, type = "l", xlab = "", ylab = "")
# title(main = "COVID data", line = 2.5)
# title(main = expression(BIC==-nT~log(tilde(sigma)^2)-Delta[K]~log(nT)), line = 1)
# title(xlab = "Number of clusters, K", line = 2.5)
# dev.off()
# 
# 
# pdf("plots/AIC_alt.pdf", width=7, height=5, paper="special")
# par(cex = 1, tck = -0.025)
# par(mar = c(3.5, 2.5, 4, 0)) #Margins for each plot
# plot(1:max_n_cl, AIC_vec, type = "l", xlab = "", ylab = "")
# title(main = "COVID data", line = 2.5)
# title(main = expression(AIC==-nT~log(tilde(sigma)^2)-Delta[K]), line = 1)
# title(xlab = "Number of clusters, K", line = 2.5)
# dev.off()


#results_output_alt(res, covid_mat, Delta_hat, n_cl, countries, path = "plots/",
#                   bw_abs, grid_points, t_len, a_vec, b_vec, c_vec, norm_p)