rm(list=ls())

library(haven)
library(tidyr)
library(MSinference)
library(dplyr)
library(zoo)

######################
#Necessary parameters#
######################
alpha     <- 0.05
sim_runs  <- 5000
q         <- 15 #Parameters for the estimation of sigma
r         <- 10


##############
#Data loading#
##############
source("functions/data_loading_hp.R") #Now matrix hp_data contains all data

#Variables
at        <- c(1, 31, 61, 91, 121)
dates     <- unique(hp_data$year)
n_ts      <- length(unique(hp_data$iso))
t_len     <- nrow(hp_data) / n_ts
countries <- unique(hp_data$iso)


##############################
#Estimation of the parameters#
##############################

#Estimating alpha and beta
source("functions/parameters_estimation.R")
estimated   <- parameters_hp(hp_data, n_ts, t_len, countries)
hp_log      <- estimated$hp_log #Original time series
hp_log_augm <- estimated$hp_log_augm   #Augmented time series

#Estimating the long-run variance
source("functions/sigma_estimation.R")
sigma_vec <- sigma(hp_log_augm, n_ts = n_ts, q_ = q, r_ = r,
                   procedure = "all_one")


#########
#Testing#
#########

#Constructing the grid
u_grid <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
h_grid <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
h_grid <- h_grid[h_grid > log(t_len) / t_len]
grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

result <- multiscale_test(data = hp_log_augm, sigma_vec = sigma_vec,
                          alpha = alpha,  n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)


#####################
#Plots for the paper#
#####################
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


source("functions/functions.R")

# #Producing plots with the final results
# for (l in seq_len(nrow(result$ijset))){
#   i <- result$ijset[l, 1]
#   j <- result$ijset[l, 2]
#   if (result$stat_pairwise[i, j] > result$quant){
#     produce_plots_hp(results = result, data_i = hp_log_augm[, i],
#                      data_j = hp_log_augm[, j], at_ = at,
#                      labels_ = dates[at], name_i = countries[i],
#                      name_j = countries[j], l = l)
#   }
# }

#Producing plots with the final results
for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  if (result$stat_pairwise[i, j] > result$quant){
    produce_plots_talk(results = result, l = l, data_i = hp_log[, i],
                     data_j = hp_log[, j], at_ = at,
                     labels_ = dates[at], dates_ = dates, name_i = countries[i],
                     name_j = countries[j])
  }
}


############
#CLUSTERING#
############
library(dendextend)

colors <- c("#EB811B", "#604c38", "#14B03D")
#colors <- c(2, 3, 4)

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
res_dend   <- as.dendrogram(res)

pdf("output/plots/talk/VOC/dendrogram.pdf", width = 15, height = 6,
    paper = "special")
par(mar = c(1.5, 1, 3, 1)) #Margins for each plot
par(oma = c(0.2, 0.5, 0.5, 0.5)) #Outer margins
plot(res, cex = 2, xlab = "", ylab = "", yaxt = "n",
     hang = 0.4, main = "")
     #main = "HAC dendrogram", xlim = c(-1, 10), cex.main = 3)
title("HAC dendrogram", cex.main = 3, line = 1)
par(lwd = 2)
rect.dendrogram(res_dend, k = 3, border = colors, lower_rect = -3.5)
dev.off()

n_cl        <- max(subgroups)
grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)

#Plotting all clusters on one plot
pdf("output/plots/talk/VOC/all_clusters.pdf", width = 7, height = 4,
    paper="special")

#Setting the layout of the graphs
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

plot(NA, ylab = "", xlab = "", xlim = c(0, 1),
     ylim = c(-1.2, 1.5), xaxt = 'n',
     mgp = c(2, 0.5, 0), cex = 0.8, tck = -0.025)  

for (cl in 1:max(subgroups)){
  locations_cluster <- colnames(Delta_hat)[subgroups == cl]
  for (column in locations_cluster){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                             MoreArgs = list(data_p = hp_log_augm[, column],
                                             grid_p = grid_points, bw = 7/t_len))
    lines(grid_points, smoothed_curve, col = colors[cl])
  }
}
axis(1, at = grid_points[at], labels = dates[at], cex = 0.7)
title(main = "Estimated time trends", line = 1, cex.main = 1.1)
dev.off()

#Plotting all clusters on one plot
pdf("output/plots/talk/VOC/all_time_series.pdf", width = 7, height = 4,
    paper="special")

#Setting the layout of the graphs
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

plot(NA, ylab = "", xlab = "", xlim = c(0, 1),
     ylim = c(-1.2, 1.5), xaxt = 'n',
     mgp = c(2, 0.5, 0), cex = 0.8, tck = -0.025)  

for (cl in 1:max(subgroups)){
  locations_cluster <- colnames(Delta_hat)[subgroups == cl]
  for (column in locations_cluster){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                             MoreArgs = list(data_p = hp_log_augm[, column],
                                             grid_p = grid_points, bw = 7/t_len))
    lines(grid_points, smoothed_curve, col = "#604c38")
  }
}
axis(1, at = grid_points[at], labels = dates[at], cex = 0.7)
title(main = "Estimated time trends", line = 1, cex.main = 1.1)
dev.off()


cluster_labels <- c(3, 1, 2)

#Plotting the trend functions of each cluster on a separate graph
for (cl in 1:max(subgroups)){
  locations_cluster <- colnames(Delta_hat)[subgroups == cl]
  pdf(paste0("output/plots/talk/VOC/cluster_", cluster_labels[cl], ".pdf"),
      width = 7, height = 4, paper = "special")
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  plot(NA, ylab = "", xlab = "", xlim = c(0, 1),
       ylim = c(-1.2, 1.5), xaxt = 'n',
       mgp = c(2, 0.5, 0), cex = 1.2, tck = -0.025)
  
  for (column in locations_cluster){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                             MoreArgs = list(data_p = hp_log_augm[, column],
                                             grid_p = grid_points, bw = 7/t_len))
    lines(grid_points, smoothed_curve, col = colors[cl])
  }
  axis(1, at = grid_points[at], labels = dates[at])
  title(main = paste0("Cluster ", cluster_labels[cl]), line = 1)
  dev.off()
}