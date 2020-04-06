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
sourceCpp("functions/multiscale_quantiles.cpp")
sourceCpp("functions/kernel_averages.cpp")


N_ts = 25
SimRuns = 1000
Tlen = 300
grid <- grid_construction(Tlen)
gset <- grid$gset

N                      <- as.integer(dim(grid$gset)[1])
gset_cpp               <- as.matrix(gset)
gset_cpp               <- as.vector(gset_cpp) 
storage.mode(gset_cpp) <- "double"

set.seed(1)
sigma_vec <- rep(1, N_ts)
tic("With C method")
a <- gaussian_stat_distr(T = Tlen, N_ts = N_ts, SimRuns = SimRuns, gset = gset_cpp, N = N,
                    sigma_vec = sigma_vec)
toc()
set.seed(1)
b <- multiscale_quantiles(T = Tlen, grid = grid, sigma_vector = sigma_vec, SimRuns= SimRuns, N_ts = N_ts)
  
all.equal(a, b)

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
Tlen               <- nrow(monthly_temp)
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
grid    <- grid_construction(Tlen)
gset    <- grid$gset
u.grid  <- sort(unique(gset[,1]))
h.grid  <- sort(unique(gset[,2]))
correct <- sqrt(2*log(1/(2*gset[,2])))

# Compute kernel weights and critical value for multiscale test
Tlen                   <- as.integer(Tlen) 
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


#Calculating the statistic for real data
statistic <- multiscale_testing(alpha, y_data, N_ts, g_t_set, sigmahat_vector_2, kernel_method)






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

