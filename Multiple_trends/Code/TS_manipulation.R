rm(list=ls())

library(ncdf4)
library(raster)
library(rgdal)
library(ggplot2)
library(reshape2)
library(multiscale)
library(dendextend)
library(corrplot)

source("functions/functions.R")

##############
#Coefficients#
##############

alpha     <- c(0.05, 0.1)
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
t   <- t[t < 2022] #deleting the last observations for year 2022

tmp  <- ncvar_get(temp_data, "temperature") # store the data in a 3-dimensional array

#Taking 5x5 grid instead of 1x1 and yearly observation instead of monthly
tmp2 <- array(NA, dim = c(dim(tmp)[1] / 5, dim(tmp)[2] / 5, dim(tmp)[3] / 12))
for (i in 1:(dim(tmp)[1] / 5)){
  for (j in 1:(dim(tmp)[2] / 5)){
    for (k in 1:(dim(tmp)[3] / 12)){
      tmp2[i, j, k] <- mean(tmp[5 * (i - 1) + 4, 5 * (j - 1) + 1, (12 * (k - 1) + 1):(12 * k)])
    }
  }
}
dimnames(tmp2) <- list(Col1 = lon[seq(4, length(lon), by = 5)],
                       Col2 = lat[seq(1, length(lat), by = 5)],
                       Col3 = floor(t[seq(1, length(t), by = 12)]))

tmp3 <- tmp2[, , 101:(dim(tmp2)[3])]
rm(tmp, tmp2)

#Checking how many time series satisfy our requirements of no NaNs
n_ts <- 0
for (i in 1:(dim(tmp3)[1])){
  for (j in 1:(dim(tmp3)[2])){
    if (sum(is.na(tmp3[i, j, ])) == 0){
      n_ts <- n_ts + 1
    }
  }
}

=
#Our main temperature matrix
dates       <- dimnames(tmp3)$Col3
t_len       <- length(dates) #Length of time series
temp_matrix <- matrix(NA, ncol = n_ts, nrow = t_len)
col_names   <- c()

k <- 1
for (i in 1:(dim(tmp3)[1])){
  for (j in 1:(dim(tmp3)[2])){
    if (sum(is.na(tmp3[i, j, ])) == 0){
      temp_matrix[, k] <- tmp3[i, j, ] - mean(tmp3[i, j, ])
      col_names <- c(col_names, paste0(dimnames(tmp3)$Col1[i], ", ",
                                       dimnames(tmp3)$Col2[j]))
      k <- k + 1
    }
  }
}

rownames(temp_matrix) <- dates
colnames(temp_matrix) <- col_names

residual_matrix <- matrix(NA, ncol = n_ts, nrow = t_len)
grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len) #For plotting

for (i in 1:n_ts){
  smoothed_curve <- mapply(local_linear_smoothing, grid_points,
                           MoreArgs = list(data_p = temp_matrix[, i],
                                           grid_p = grid_points,
                                           bw = 7 / t_len))
  residual_matrix[, i] <- temp_matrix[, i] - smoothed_curve
}

pearson_matrix <- matrix(NA, ncol = n_ts, nrow = n_ts)
for (i in 1:(n_ts - 1)){
  pearson_matrix[i, i] <- 1
  for (j in (i + 1):n_ts){
    mu_i <- mean(residual_matrix[, i])    
    mu_j <- mean(residual_matrix[, j])
    demeaned_i <- residual_matrix[, i] - mu_i
    demeaned_j <- residual_matrix[, j] - mu_j
    pearson_matrix[i, j] <- (sum(demeaned_i * demeaned_j)) / (sqrt(sum(demeaned_i ^ 2) * sum(demeaned_j ^ 2)))
    pearson_matrix[j, i] <- pearson_matrix[i, j]
  }
}

pdf("plots/clustering/pearson_correlation.pdf", width = 30, height = 30,
    paper = "special")
#par(cex = 1, tck = -0.025)
par(mar = c(0, 0, 0, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
corrplot(pearson_matrix, is.corr = FALSE, method = "color", tl.pos = "n",  cl.cex = 4.0)
dev.off()

#####################
#Estimating variance#
#####################

sigmahat_vector <- c()
for (i in 1:n_ts){
  AR.struc        <- estimate_lrv(data = temp_matrix[, i], q = q, r_bar = r,
                                  #p = order_results[i])  
                                  p = 1)
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}
