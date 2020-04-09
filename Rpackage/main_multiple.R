# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

#library(microbenchmark)
library(Rcpp)
library(tictoc)
#library(xtable)
#options(xtable.floating = FALSE)
#options(xtable.timestamp = "")

# The following file contains functions that are used to estimate the AR parameters, in particular, to compute
# the estimators of the coefficients, a_1, ... ,a_p, the estimator of the variance of the innovation term, sigma_eta^2 
# and the estimator of the long-run variance sigma^2.
source("functions/long_run_variance.r")

#Load necessary functions  
source("functions/ConstructGrid.r")
source("functions/multiscale_statistics.r")
source("functions/multiscale_quantiles.r")
source("functions/multiscale_testing.r")
source("functions/minimal_intervals.r")
source("functions/functions.r")
sourceCpp("functions/multiscale_statistics.cpp")


##############################
#Defining necessary constants#
##############################

N_ts    <- 15 #number of different time series for application to analyze
alpha   <- 0.05 #confidence level for application
SimRuns <- 100


###########################################
#Loading the real station data for England#
###########################################

data_frame <- read.csv("data/returns.csv", stringsAsFactors = FALSE)
a <- split(data_frame[, c(1, 2, 9)], data_frame$PERMNO)


for (i in 1:N_ts){
  returns_tmp <- a[[i]][, c(2, 3)]
  returns_tmp[,2] <- log(as.numeric(returns_tmp[, 2])^2)
  colnames(returns_tmp) <- c('DATE', paste0("returns", i))
  if (i == 1){
    returns <- returns_tmp
  } else {
    returns <- merge(returns, returns_tmp, by = 'DATE', all.x = TRUE, all.y = TRUE)
  }
}

rm(a, returns_tmp)

returns <- as.matrix(returns)

returns <- returns[,colSums(is.na(returns)) <= 2] #Ommitting the time series with too sparse data
returns[!is.finite(returns)] <- NA
returns <- na.omit(returns)#Deleting the rows with ommitted variables

for (i in 2:(N_ts+1)){
  returns[, i] <- returns[, i] - mean(returns[, i])
}

Tlen          <- nrow(returns)
N_ts          <- ncol(returns) - 1 #Updating the number of time series because of dropped stations

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


#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 2:(N_ts+1)){
  AR.struc        <- AR_lrv(data = returns[, i], q = q, r.bar = r.bar, p=order)
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}



#Calculating the statistic for real data
grid <- grid_construction(Tlen)
gset <- grid$gset

result <- multiscale_testing(alpha, returns[, -1], grid, sigma_vec = sigmahat_vector, N_ts = N_ts)





################
#Testing itself#
################

y_data = monthly_temp[-c(1, 2, 3)]
month_column = monthly_temp['month']


  #Calculating smoothed curve for the data using local linear estimator#
  grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for estimating

  #Calculating gaussian quantile for T_tempr#
  
  g_t_set <- creating_g_set(T_tempr, kernel_method)
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

