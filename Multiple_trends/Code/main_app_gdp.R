rm(list=ls())

library(tidyr)
library(MSinference)
library(dplyr)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


######################
#Necessary parameters#
######################
alpha     <- 0.05
sim_runs  <- 5000
q         <- 20 #Parameters for the estimation of sigma
r         <- 10
at        <- seq(5, 125, by = 20) #Ticks and labels on the x-axis for plotting

##############
#Data loading#
##############
load(file = "data/gdp_data.RData")
#Now matrix X_mat_filled contains all the necessary data

#Variables
dates <- unique(X_mat_filled$date)
n_ts  <- length(unique(X_mat_filled$WBcode))
t_len <- nrow(X_mat_filled) / n_ts


##############################
#Estimation of the parameters#
##############################

#Estimating alpha and beta
source("functions/parameters_estimation.R")
estimated      <- parameters_gdp(X_mat_filled, n_ts, t_len, countries)
gdp_mat_growth <- estimated$gdp_mat_growth #Original time series
gdp_mat_augm   <- estimated$gdp_mat_augm   #Augmented time series

#Estimating the long-run variance
source("functions/sigma_estimation.R")
sigma_vec <- sigma(gdp_mat_augm, n_ts = n_ts, q_ = q, r_ = r,
                   procedure = "BIC")


#########
#Testing#
#########
#Constructing the grid
u_grid <- seq(from = 9 / (2 * t_len), to = 1, by = 4 / t_len)
h_grid <- seq(from = 4 / t_len, to = 1 / 4, by = 4 / t_len)
h_grid <- h_grid[h_grid > log(t_len) / t_len]
grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

result <- multiscale_test(data = gdp_mat_augm, sigma_vec = sigma_vec,
                          alpha = alpha,  n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs)

#####################
#Plots for the paper#
#####################
source("functions/functions.R")

for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  if (result$stat_pairwise[i, j] > result$quant) {
    if ((countries[i] == "NOR") | ((countries[i] == "FRA") & (countries[j] != "NOR"))) {
      #For color consistency in the paper
      i <- result$ijset[l, 2]
      j <- result$ijset[l, 1]
    }
    produce_plots_gdp(results = result, data_i = gdp_mat_augm[, i],
                      data_j = gdp_mat_augm[, j], ticks_ = at,
                      labels_ = dates[at], name_i = countries[i],
                      name_j = countries[j])
  }
}

############
#Clustering#
############
#for the distance matrix we need a symmetrical one
Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
for (i in 1:(n_ts - 1)){
  for (j in (i + 1):n_ts){
    Delta_hat[i, j] <- result$stat_pairwise[i, j]
    Delta_hat[j, i] <- Delta_hat[i, j]
  }
}

colnames(Delta_hat) <- countries
rownames(Delta_hat) <- countries

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
subgroups  <- cutree(res, h = result$quant)

pdf("output/plots/gdp/dendrogram.pdf", width = 15, height = 6,
    paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
plot(res, cex = 1.5, xlab = "", ylab = "", yaxt = "n", hang = 0.55,
     main = "HAC dendrogram", cex.main = 2.2)
par(lwd = 2)
rect.hclust(res, h = result$quant, border = 2:(max(subgroups) + 1))
dev.off()

n_cl        <- max(subgroups)
grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
colors      <- c(3, 4, 2)

#Plotting all clusters on one plot
pdf("output/plots/gdp/all_clusters.pdf", width = 7, height = 4,
    paper="special")
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

plot(NA, ylab = "", xlab = "", xlim = c(0, 1),
     ylim = c(-0.027, 0.015), xaxt = 'n',
     mgp = c(2, 0.5, 0), cex = 0.8, tck = -0.025)  

for (cl in 1:max(subgroups)){
  locations_cluster <- colnames(Delta_hat)[subgroups == cl]
  
  for (column in locations_cluster){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                             MoreArgs = list(data_p = gdp_mat_augm[, column],
                                             grid_p = grid_points, bw = 0.1))
    lines(grid_points, smoothed_curve, col = colors[cl])
  }
}
axis(1, at = grid_points[at], labels = dates[at], cex = 0.7)
title(main = "Estimated time trends", line = 1, cex.main = 1.1)
dev.off()

#Plotting the trend functions of each cluster on a separate graph
for (cl in 1:max(subgroups)){
  locations_cluster <- colnames(Delta_hat)[subgroups == cl]
  pdf(paste0("output/plots/gdp/cluster_", cl, ".pdf"),
      width = 7, height = 4, paper = "special")
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  plot(NA, ylab = "", xlab = "", xlim = c(0, 1),
       ylim = c(-0.027, 0.015), xaxt = 'n',
       mgp = c(2, 0.5, 0), cex = 1.2, tck = -0.025)
  
  for (column in locations_cluster){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                             MoreArgs = list(data_p = gdp_mat_augm[, column],
                                             grid_p = grid_points, bw = 0.1))
    lines(grid_points, smoothed_curve, col = colors[cl])
  }
  axis(1, at = grid_points[at], labels = dates[at])
  title(main = paste0("Cluster ", cl), line = 1)
  dev.off()
}
