rm(list=ls())

library(fda.usc)
library(MSinference)

######################
#Necessary parameters#
######################
alpha     <- 0.05
sim_runs  <- 1000
q         <- 20 #Parameters for the estimation of sigma
r         <- 15


##############
#Data loading#
##############
data(aemet)
wind_data <- t(aemet$wind.speed[['data']])
temp_data <- t(aemet$temp[['data']])
prec_data <- t(aemet$logprec[['data']])

rownames(temp_data) <- 1:365
rownames(wind_data) <- 1:365
rownames(prec_data) <- 1:365

#Variables
dates     <- rownames(temp_data)
n_ts      <- ncol(temp_data)
t_len     <- nrow(temp_data)
stations  <- gsub("1980-2009", "", colnames(temp_data))

stations_df <- data.frame("lat"= aemet$df$latitude,
                          "lon" = aemet$df$longitude,
                          row.names = stations)

colnames(temp_data) <- stations
colnames(wind_data) <- stations
colnames(prec_data) <- stations


##############################
#Estimation of the parameters#
##############################
temp_data_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(temp_data_augm) <- stations

wind_data_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(wind_data_augm) <- stations

prec_data_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(prec_data_augm) <- stations

alpha_temp_vec <- c()
alpha_wind_vec <- c()
alpha_prec_vec <- c()

for (i in 1:n_ts){
  alpha_temp_vec[i] <- mean(temp_data[, i])
  temp_data_augm[, i] <- temp_data[, i] - alpha_temp_vec[i]
  
  alpha_wind_vec[i] <- mean(wind_data[, i])
  wind_data_augm[, i] <- wind_data[, i] - alpha_wind_vec[i]

  alpha_prec_vec[i] <- mean(prec_data[, i])
  prec_data_augm[, i] <- prec_data[, i] - alpha_prec_vec[i]
}

#Estimating the long-run variance
source("functions/sigma_estimation.R")
sigma_temp_vec <- sigma(temp_data_augm, n_ts = n_ts, q_ = q, r_ = r,
                        procedure = "BIC")

sigma_wind_vec <- sigma(wind_data_augm, n_ts = n_ts, q_ = q, r_ = r,
                        procedure = "BIC")

sigma_prec_vec <- sigma(prec_data_augm, n_ts = n_ts, q_ = q, r_ = r,
                        procedure = "BIC")


#########
#Testing#
#########

u_grid <- seq(from = 11 / (2 * t_len), to = 1, by = 5 / t_len)
h_grid <- seq(from = 5 / t_len, to = 1 / 12, by = 5 / t_len)
h_grid <- h_grid[h_grid > log(t_len) / t_len]
grid   <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

result_temp <- multiscale_test(data = temp_data_augm, sigma_vec = sigma_temp_vec,
                               alpha = alpha,  n_ts = n_ts, grid = grid,
                               sim_runs = sim_runs, epidem = FALSE)

result_wind <- multiscale_test(data = wind_data_augm, sigma_vec = sigma_wind_vec,
                               alpha = alpha,  n_ts = n_ts, grid = grid,
                               sim_runs = sim_runs, epidem = FALSE)

result_prec <- multiscale_test(data = prec_data_augm, sigma_vec = sigma_prec_vec,
                               alpha = alpha,  n_ts = n_ts, grid = grid,
                               sim_runs = sim_runs, epidem = FALSE)


############
#CLUSTERING#
############
library(dendextend)
source("functions/functions.R")

#aggregating distance matrices plus making a symmetrical one
Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
for (i in 1:(n_ts - 1)){
  for (j in (i + 1):n_ts){
    Delta_hat[i, j] <- max(result_temp$stat_pairwise[i, j],
                           result_wind$stat_pairwise[i, j],
                           result_prec$stat_pairwise[i, j])
    Delta_hat[j, i] <- Delta_hat[i, j]
  }
}

colnames(Delta_hat) <- stations
rownames(Delta_hat) <- stations


delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
subgroups  <- cutree(res, h = result_temp$quant)
res_dend   <- as.dendrogram(res)

par(mar = c(1.5, 1, 3, 1)) #Margins for each plot
par(oma = c(0.2, 0.5, 0.5, 0.5)) #Outer margins
plot(res, cex = 0.5, xlab = "", ylab = "", yaxt = "n",
     hang = 0.4, main = "")
title("HAC dendrogram", cex.main = 1, line = 1)
par(lwd = 2)
rect.dendrogram(res_dend, h = result_temp$quant,
                border = 2:(max(subgroups) + 1), lower_rect = -5)

stations_df['cluster'] <- subgroups 

library(maps)
library(mapdata)
map("worldHires","Spain", xlim = c(-13.5, 8.5), ylim = c(34, 45), col="gray90", fill=TRUE)
points(x = stations_df$lon, y = stations_df$lat, pch = 20, col = subgroups)

library(ggplot2)
world_map <- map_data("world")
world_map <- subset(world_map,
                    (world_map$long > -20) & (world_map$long < 5) & (world_map$lat> 25) & (world_map$lat < 45))

#Add map to base plot
base_spain <- ggplot() +
  coord_fixed() +
  xlab("") +
  ylab("") +
  geom_polygon(data=world_map, aes(x=long, y=lat, group=group),
               colour="light green", fill="light green")

map_data_coloured <- 
  base_spain +
  geom_point(data=stations_df, 
             aes(x = lon, y= lat, colour=cluster),
             size=5, alpha=I(0.7))

map_data_coloured

n_cl        <- max(subgroups)
grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)

#Plotting all clusters on one plot
pdf("output/plots/temp/all_clusters.pdf", width = 7, height = 4,
    paper="special")

#Setting the layout of the graphs
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

plot(NA, ylab = "", xlab = "", xlim = c(0, 1),
     ylim = c(-10, 12), xaxt = 'n',
     mgp = c(2, 0.5, 0), cex = 0.8, tck = -0.025)  

for (cl in 1:max(subgroups)){
  locations_cluster <- colnames(Delta_hat)[subgroups == cl]
  for (column in locations_cluster){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                             MoreArgs = list(data_p = temp_data_augm[, column],
                                             grid_p = grid_points, bw = 10/t_len))
    lines(grid_points, smoothed_curve, col = cl + 1)
  }
}
axis(1, cex = 0.7)
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
                                             grid_p = grid_points, bw = 10/t_len))
    lines(grid_points, smoothed_curve, col = cl + 1)
  }
  axis(1, at = grid_points[at], labels = dates[at])
  title(main = paste0("Cluster ", cluster_labels[cl]), line = 1)
  dev.off()
}