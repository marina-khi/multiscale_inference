# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

#library(microbenchmark)
library(Rcpp)
library(tictoc)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

# The following file contains functions that are used to estimate the AR parameters, in particular, to compute
# the estimators of the coefficients, a_1, ... ,a_p, the estimator of the variance of the innovation term, sigma_eta^2 
# and the estimator of the long-run variance sigma^2.
source("functions/long_run_variance.r")

#Load necessary functions  
source("functions/ConstructGrid.r")
source("functions/multiscale_statistics.r")
source("functions/multiscale_statistics_old.r")
source("functions/multiscale_quantiles.r")
source("functions/multiscale_testing.r")
source("functions/minimal_intervals.r")
source("functions/functions.r")
sourceCpp("functions/kernel_averages.cpp")
sourceCpp("functions/kernel_weights.cpp")



##############################
#Defining necessary constants#
##############################

N_ts          <- 34 #number of different time series for application to analyze
alpha         <- 0.05 #confidence level for application


###########################################
#Loading the real station data for England#
###########################################

for (i in 1:N_ts){
  filename = paste("data/txt", i, ".txt", sep = "")
  temperature_tmp  <- read.table(filename, header = FALSE, skip = 7,
                                 col.names = c("year", "month", "tmax", "tmin", "af", "rain", "sun", "aux"), fill = TRUE,  na.strings = c("---"))
  monthly_temp_tmp <- data.frame('1' = as.numeric(temperature_tmp[['year']]), '2' = as.numeric(temperature_tmp[['month']]),
                                 '3' = (temperature_tmp[["tmax"]] + temperature_tmp[["tmin"]]) / 2)
  colnames(monthly_temp_tmp) <- c('year', 'month', paste0("tmean", i))
  if (i == 1){
    monthly_temp <- monthly_temp_tmp
  } else {
    monthly_temp <- merge(monthly_temp, monthly_temp_tmp, by = c("year", "month"), all.x = TRUE, all.y = TRUE)
  }
}

rm(monthly_temp_tmp, temperature_tmp)

monthly_temp <- subset(monthly_temp, year >= 1986) #Subsetting years 1986 - 2018 because of closed and new stations
monthly_temp <- monthly_temp[,colSums(is.na(monthly_temp)) <= 2] #Ommitting the time series with too sparse data
monthly_temp <- na.omit(monthly_temp)#Deleting the rows with ommitted variables

date               <- paste(sprintf("%02d", monthly_temp$month), monthly_temp$year,  sep='-')
monthly_temp       <- cbind(date, monthly_temp)
TemperatureColumns <- setdiff(names(monthly_temp), c("year", "month", "date"))
T_tempr            <- nrow(monthly_temp)
N_ts               <- ncol(monthly_temp) - 3 #Updating the number of time series because of dropped stations


######################
#Deseasonalizing data#
######################

#monthly_temp[4:(N_ts + 3)] <- lapply(monthly_temp[4:(N_ts + 3)], function(x) x - ave(x, monthly_temp[['month']], FUN=mean))
#monthly_temp[TemperatureColumns] <- lapply(monthly_temp[TemperatureColumns], function(x) x - ave(x, monthly_temp[['month']], FUN=mean))


#####################
#Estimating variance#
#####################

#Setting tuning parameters for testing
order <- 1
q     <- 25
r.bar <- 10

#Construct grid
grid    <- grid_construction(T_tempr)
gset    <- grid$gset
u.grid  <- sort(unique(gset[,1]))
h.grid  <- sort(unique(gset[,2]))
correct <- sqrt(2*log(1/(2*gset[,2])))

# Compute kernel weights and critical value for multiscale test
T_tempr                <- as.integer(T_tempr) 
N                      <- as.integer(dim(grid$gset)[1])
gset_cpp               <- as.matrix(grid$gset)
gset_cpp               <- as.vector(gset_cpp) 
storage.mode(gset_cpp) <- "double"

#Calculating each sigma_i separately
sigmahat_vector_2 <- c()
for (i in TemperatureColumns){
  AR.struc          <- AR_lrv(data = monthly_temp[[i]], q = q, r.bar = r.bar, p=order)
  sigma_hat_i       <- sqrt(AR.struc$lrv)
  sigmahat_vector_2 <- c(sigmahat_vector_2, AR.struc$lrv)
}

SimRuns <- 10
Psi1    <- rep(NA,SimRuns)


cat("","\n")
cat("Computing critical value of the multiscale test","\n")
progbar <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")
tic("New method")
T = T_tempr
if(N_ts == 1){ 
  for(pos in 1:SimRuns){
    zeta     <- rnorm(T, mean=0, sd=1)
    res      <- multiscale_statistics(data=zeta, sigmahat=1, grid=grid) 
    Psi[pos] <- res$stat_ms
    setTxtProgressBar(progbar, pos)
  }
} else {
  for (pos in 1:SimRuns){
    Psi.ij <- matrix(NA, nrow = N_ts, ncol = N_ts)      
    z_matrix <- matrix(rnorm(T * N_ts, mean=0, sd=1), nrow = T, ncol = N_ts) 
    for(i in 1:(N_ts-1)){
      for (j in (i+1):N_ts){
        z_diff <- sqrt(sigmahat_vector_2[i]) * (z_matrix[, i] - mean (z_matrix[, i])) - sqrt(sigmahat_vector_2[j]) * (z_matrix[, j] - mean (z_matrix[, j]))
        Psi.ij[i, j] <- multiscale_statistics(data=z_diff, sigmahat=sqrt(sigmahat_vector_2[i] + sigmahat_vector_2[j]), grid=grid)$stat_ms
        if ((i == 1) & (j == 2)){
          Psi <- Psi.ij[i, j]
        } else if (Psi.ij[i, j] > Psi){
          Psi <- Psi.ij[i, j]          
        }
      }
    }
    Psi1[pos] <- Psi
    setTxtProgressBar(progbar, pos)
  }
}
close(progbar)
toc()



Psi2    <- rep(NA,SimRuns)

cat("","\n")
cat("Computing critical value of the multiscale test","\n")
progbar <- txtProgressBar(min = 1, max = SimRuns, style = 3, char = ".")
tic("Old method")
T = T_tempr
if(N_ts == 1){ 
  for(pos in 1:SimRuns){
    zeta     <- rnorm(T, mean=0, sd=1)
    res      <- multiscale_statistics(data=zeta, sigmahat=1, grid=grid) 
    Psi[pos] <- res$stat_ms
    setTxtProgressBar(progbar, pos)
  }
} else {
  wghts <- matrix(kernel_weights(T_tempr, gset_cpp, N), ncol = T_tempr, byrow = TRUE)
  for (pos in 1:SimRuns){
    Psi.ij <- matrix(NA, nrow = N_ts, ncol = N_ts)      
    z_matrix <- matrix(rnorm(T * N_ts, mean=0, sd=1), nrow = T, ncol = N_ts) 
    for(i in 1:(N_ts-1)){
      for (j in (i+1):N_ts){
        z_diff <- sqrt(sigmahat_vector_2[i]) * (z_matrix[, i] - mean (z_matrix[, i])) - sqrt(sigmahat_vector_2[j]) * (z_matrix[, j] - mean (z_matrix[, j]))
        Psi.ij[i, j] <- multiscale_statistics_old(data=z_diff, weights=wghts, sigmahat=sqrt(sigmahat_vector_2[i] + sigmahat_vector_2[j]), grid=grid)$stat_ms
      }
    }
    Psi2[pos] <- max(Psi.ij, na.rm = TRUE)
    setTxtProgressBar(progbar, pos)
  }
}
close(progbar)
toc()




################
#Testing itself#
################

y_data = monthly_temp[-c(1, 2, 3)]
month_column = monthly_temp['month']


data <- as.vector(y_data$tmean2)
b <- kernel_averages(T_tempr, gset_cpp, data, N)


  #Calculating smoothed curve for the data using local linear estimator#
  grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for estimating

  #Calculating gaussian quantile for T_tempr#
  
  g_t_set <- creating_g_set(T_tempr, kernel_method)
  filename = paste("quantiles/distr_multiple_T_", T_tempr, "_and_N_", N_ts, ".RData", sep = "")
  if(!file.exists(filename)) {
    z_matrix <- matrix(0, nrow = T_tempr, ncol = N_ts + 1)
    z_matrix <- cbind(z_matrix, month_column) #Adding auxiliary month column to "deseasonalize" simulated data
    z_matrix[1:N_ts] <- lapply(z_matrix[1:N_ts], function(x) x - ave(x, z_matrix[N_ts + 1], FUN=mean))
    z_matrix <- z_matrix[-(N_ts + 1)] #Dropping auxiliary month column
    psistar_statistic(z_matrix, N_ts, g_t_set, sigmahat_vector_2, kernel_method)
    save(gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }

  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)  
  cat("Gaussian quantile is", gaussian_quantile, "\n")
  
  
  #Calculating the statistic for real data#
  statistic <- psihat_statistic(y_data, N_ts, g_t_set, sigmahat_vector_2, kernel_method)
  
  #And now the testing itself
  if (statistic[[2]] > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", statistic[[2]],
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", statistic[[2]],
        "Gaussian quantile value = ", gaussian_quantile, "\n")
  }
  
#return(list(statistic[[2]], gaussian_quantile))

