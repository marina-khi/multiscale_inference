rm(list=ls())

library(ncdf4)
library(raster)
library(rgdal)
library(ggplot2)
library(reshape2)
library(multiscale)
library(dendextend)

source("functions/functions.R")

##############
#Coefficients#
##############

alpha     <- 0.05
sim_runs  <- 1000
q         <- 25 #Parameters for the estimation of sigma
r         <- 10

##############
#Data loading#
##############

temp_data <- nc_open('data/Complete_TAVG_LatLong1.nc')

lon <- ncvar_get(temp_data, "longitude")
lat <- ncvar_get(temp_data, "latitude", verbose = F)
t   <- ncvar_get(temp_data, "time")

tmp  <- ncvar_get(temp_data, "temperature") # store the data in a 3-dimensional array

#Taking 5x5 grid instead of 1x1 and yearly observation instead of monthly
tmp2 <- array(NA, dim = c(dim(tmp)[1] / 5, dim(tmp)[2] / 5, dim(tmp)[3] / 12))
for (i in 1:(dim(tmp)[1] / 5)){
  for (j in 1:(dim(tmp)[2] / 5)){
    for (k in 1:(dim(tmp)[3] / 12)){
      tmp2[i, j, k] <- mean(tmp[(5 *(i - 1) + 1):(5 * i), (5 *(j - 1) + 1):(5 * j),
                                (12 * (k - 1) + 1):(12 * k)])
    }
  }
}

dimnames(tmp2) <- list(Col1 = lon[seq(3, length(lon), by = 5)],
                       Col2 = lat[seq(3, length(lat), by = 5)],
                       Col3 = floor(t[seq(1, length(t), by = 12)]))
rm(tmp)

#Checking how many time series satisfy our requirements of no NaNs
n_ts <- 0
for (i in 1:(dim(tmp2)[1])){
  for (j in 1:(dim(tmp2)[2])){
    if (sum(is.na(tmp2[i, j, 101:(dim(tmp2)[3])])) == 0){
      n_ts <- n_ts + 1
    }
  }
}

#Our main temperature matrix
dates       <- dimnames(tmp2)$Col3[101:(dim(tmp2)[3])]
t_len       <- length(dates) #Length of time series
temp_matrix <- matrix(NA, ncol = n_ts, nrow = t_len)
col_names   <- c()

k <- 1
for (i in 1:(dim(tmp2)[1])){
  for (j in 1:(dim(tmp2)[2])){
    if (sum(is.na(tmp2[i, j, 101:(dim(tmp2)[3])])) == 0){
      temp_matrix[, k] <- tmp2[i, j, 101:(dim(tmp2)[3])] - mean(tmp2[i, j, 101:(dim(tmp2)[3])])
      col_names <- c(col_names, paste0(dimnames(tmp2)$Col1[i], ", ",
                                       dimnames(tmp2)$Col2[j]))
      k <- k + 1
    }
  }
}


rownames(temp_matrix) <- dates
colnames(temp_matrix) <- col_names

#####################
#Estimating variance#
#####################

# #Order selection
# order_results <- order_selection(matrix = temp_matrix, q_vec = 10:20,
#                                  r_vec = 10:15)

#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 1:n_ts){
  AR.struc        <- estimate_lrv(data = temp_matrix[, i], q = q, r_bar = r,
                                  #p = order_results[i])  
                                  p = 1)
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}

#Constructing the grid
grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len) #For plotting
at_         <- seq(from = 1, to = t_len, by = 20)

u_grid      <- (3:t_len)[c(TRUE, FALSE, TRUE, rep(FALSE, 2))]/t_len
#u_grid      <- seq(from = 5 / t_len, to = 1, by = 1 / t_len)
h_grid      <- seq(from = 2 / t_len, to = 1 / 4, by = 5 / t_len)
h_grid      <- h_grid[h_grid > log(t_len) / t_len]
grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

format(Sys.time(), "%a %b %d %X %Y")
result <- multiscale_test(data = temp_matrix, sigma_vec = sigmahat_vector,
                          alpha = alpha,  n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)
format(Sys.time(), "%a %b %d %X %Y")

save(result, file = "result_long_ts.RData")
#load(file = "result_long_ts.RData")

###########################
#CLustering of the results#
###########################

#for the distance matrix we need a symmetrical one
Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)
for (i in 1:(n_ts - 1)){
  for (j in (i + 1):n_ts){
    Delta_hat[i, j] <- result$stat_pairwise[i, j]
    Delta_hat[j, i] <- Delta_hat[i, j]
  }
}

colnames(Delta_hat) <- col_names
rownames(Delta_hat) <- col_names

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
subgroups  <- cutree(res, h = result$quant)

#Dendrogram
pdf("plots/dendrogram_long_ts.pdf", width = 15, height = 6, paper = "special")
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins
plot(res, cex = 0.8, xlab = "", ylab = "")
rect.hclust(res, h = result$quant, border = 2:(max(subgroups) + 1))
dev.off()


#Plotting all clusters on one plot
pdf("plots/results_all_clusters_long_ts.pdf", width=7, height=6, paper="special")

#Setting the layout of the graphs
par(cex = 1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

plot(NA, ylab = "", xlab = "", xlim = c(0, 1),
     ylim = c(-1.5, 3.1), xaxt = 'n',
     mgp = c(2, 0.5, 0), cex = 1.2, tck = -0.025)  

for (cl in 1:max(subgroups)){
  locations_cluster <- colnames(Delta_hat)[subgroups == cl]
  
  for (column in locations_cluster){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                             MoreArgs = list(data_p = temp_matrix[, column],
                                             grid_p = grid_points, bw = 0.1))
    lines(grid_points, smoothed_curve, col = c("red", "black", "blue")[cl])
  }
  axis(1, at = grid_points[at_], labels = dates[at_])
}
title(main = "All clusters", line = 1)
dev.off()


#Plotting the trend functions of each cluster on a separate graph
for (cl in 1:max(subgroups)){
  locations_cluster <- colnames(Delta_hat)[subgroups == cl]
  pdf(paste0("plots/results_cluster_", cl, "_long_ts.pdf"), width = 7,
      height = 6, paper="special")
  
  #Setting the layout of the graphs
  par(cex = 1, tck = -0.025)
  par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
  par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
  plot(NA, ylab = "", xlab = "", xlim = c(0, 1),
       ylim = c(-1.5, 3.1), xaxt = 'n',
       mgp = c(2, 0.5, 0), cex = 1.2, tck = -0.025)
  
  for (column in locations_cluster){
    smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                             MoreArgs = list(data_p = temp_matrix[, column],
                                             grid_p = grid_points, bw = 0.1))
    lines(grid_points, smoothed_curve, col = c("red", "black", "blue")[cl])
  }
  axis(1, at = grid_points[at_], labels = dates[at_])
  title(main = paste0("Cluster ", cl), line = 1)
  dev.off()
}


#Plotting world map
temp_tmp           <- data.frame(col_names)
temp_tmp$cluster   <- cutree(res, n_cl)
temp_tmp$lon       <- as.double(sub(", .*", "", temp_tmp$col_names))
temp_tmp$lat       <- as.double(sub(".*, ", "", temp_tmp$col_names))
temp_tmp$col_names <- NULL

temp_map <- reshape(temp_tmp, idvar = "lon", timevar = "lat", v.names = "cluster",
                    direction = "wide", sep = "_")
rownames(temp_map) <- temp_map$lon
colnames(temp_map) <- sub('^cluster_', '', colnames(temp_map))
temp_map$lon <- NULL
temp_map     <- as.matrix(temp_map)

rm(temp_tmp)

temp_slice <- tmp2[, , 1]
r <- raster(t(temp_map), xmn = min(lon), xmx = max(lon), ymn = min(lat),
            ymx = max(lat),
            crs = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
#r <- flip(r, direction='x')
plot(r)

