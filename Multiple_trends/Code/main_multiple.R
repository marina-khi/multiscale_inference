# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

library(multiscale)
library(tictoc)
#library(xtable)
#options(xtable.floating = FALSE)
#options(xtable.timestamp = "")

##############################
#Defining necessary constants#
##############################

alpha    <- 0.05 #confidence level for application
sim_runs <- 5000


###########################################
#Loading the real station data for England#
###########################################

data_frame <- read.csv("data/returns.csv", stringsAsFactors = FALSE)
data_frame$RETX <- as.numeric(data_frame$RETX)
a <- split(data_frame[, c("PERMNO", "DATE", "RETX")], data_frame$PERMNO)

n_ts <- length(a)

for (i in 1:n_ts){
  returns_tmp <- a[[i]][, c(2, 3)]
  colnames(returns_tmp) <- c('DATE', paste0("returns", i))
  returns_tmp[, 2] <- log((returns_tmp[, 2] - mean(returns_tmp[, 2], na.rm = TRUE))^2)
  if (i == 1){
    returns <- returns_tmp
  } else {
    returns <- merge(returns, returns_tmp, by = 'DATE', all.x = TRUE, all.y = TRUE)
  }
}

rm(a, returns_tmp)

returns <- as.matrix(returns)

returns <- returns[,colSums(is.na(returns)) <= 600] #Ommitting the time series with too sparse data
returns[!is.finite(returns)] <- NA
returns <- na.omit(returns)#Deleting the rows with ommitted variables

t_len         <- nrow(returns)
n_ts          <- ncol(returns) - 1 #Updating the number of time series because of dropped stations

for (i in 2:(n_ts+1)){
  returns[, i] <- returns[, i] - mean(returns[, i])
}


#####################
#Estimating variance#
#####################

#Order selection
q <- 20:55
r <- 10:15

j = 3
criterion_matrix <- expand.grid(q = q, r = r)

criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))

for (i in 1:nrow(criterion_matrix)){
  FPE <- c()
  AIC <- c()
  AICC <- c()
  SIC <- c()
  HQ <- c()
  
  different_orders <- (1:9)
  
  for (order in different_orders){
    AR.struc      <- AR_lrv(data=returns[, j], q=criterion_matrix$q[[i]], r.bar=criterion_matrix$r[[i]], p=order)
    sigma_eta_hat <- sqrt(AR.struc$vareta)
    FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
    AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
    AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order +2)/t_len))
    SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
    HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
  }
  criterion_matrix$FPE[[i]]  <- which.min(FPE)
  criterion_matrix$AIC[[i]]  <- which.min(AIC)
  criterion_matrix$AICC[[i]] <- which.min(AICC)
  criterion_matrix$SIC[[i]]  <- which.min(SIC)
  criterion_matrix$HQ[[i]]   <- which.min(HQ)
}
cat("For stock ", colnames(returns)[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")





#Setting tuning parameters for testing
order <- c(1, 3, 2, 6, 4, 1, 4, 1, 1, 1)
q     <- 50
r.bar <- 10


#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 2:(n_ts+1)){
  AR.struc        <- AR_lrv(data = returns[, i], q = q, r.bar = r.bar, p=order[i-1])
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}



#Calculating the statistic for real data
result <- multiscale_testing(alpha = alpha, data = returns[, -1], sigma_vec = sigmahat_vector, SimRuns = SimRuns, n_ts = n_ts)

#And now the testing itself
if (max(result$Psi_ij) > result$quant) {
  cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", max(result$Psi_ij),
      "Gaussian quantile value = ", result$quant, "\n")
} else {
  cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", max(result$Psi_ij),
      "Gaussian quantile value = ", result$quant, "\n")
}



#############
#Temperature#
#############
# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

library(multiscale)
library(tictoc)
#library(xtable)
#options(xtable.floating = FALSE)
#options(xtable.timestamp = "")

n_ts          <- 34 #number of different time series for application to first analyze
alpha         <- 0.05 #confidence level for application
sim_runs      <- 1000

###########################################
#Loading the real station data for England#
###########################################

for (i in 1:n_ts){
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

rm(monthly_temp_tmp)

monthly_temp <- subset(monthly_temp, year >= 1986) #Subsetting years 1986 - 2018 because of closed and new stations
monthly_temp <- monthly_temp[,colSums(is.na(monthly_temp)) <= 2] #Ommitting the time series with too sparse data
monthly_temp <- na.omit(monthly_temp)#Deleting the rows with ommitted variables

date               <- paste(sprintf("%02d", monthly_temp$month), monthly_temp$year,  sep='-')
monthly_temp       <- cbind(date, monthly_temp)
TemperatureColumns <- setdiff(names(monthly_temp), c("year", "month", "date"))
T_tempr            <- nrow(monthly_temp)
n_ts               <- ncol(monthly_temp) - 3 #Updating the number of time series because of dropped stations


######################
#Deseasonalizing data#
######################

monthly_temp[4:(n_ts + 3)] <- lapply(monthly_temp[4:(n_ts + 3)], function(x) x - ave(x, monthly_temp[['month']], FUN=mean))
monthly_temp[TemperatureColumns] <- lapply(monthly_temp[TemperatureColumns], function(x) x - ave(x, monthly_temp[['month']], FUN=mean))


#####################
#Estimating variance#
#####################

#Tuning parameters
order <- 1
q     <- 25
r_bar <- 10

#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in TemperatureColumns){
  AR.struc        <- estimate_lrv(data = monthly_temp[[i]], q = q, r_bar = r_bar, p=order)
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}

#Constructing the grid
grid <- construct_grid(t = T_tempr)

#Calculating the statistic for real data
monthly_temp <- do.call(cbind, monthly_temp)
result <- multiscale_test(data = monthly_temp[, -c(1, 2, 3)], sigma_vec = sigmahat_vector,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

#And now the testing itself
if (max(result$Psi_ij) > result$quant) {
  cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", max(result$Psi_ij),
      "Gaussian quantile value = ", result$quant, "\n")
} else {
  cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", max(result$Psi_ij),
      "Gaussian quantile value = ", result$quant, "\n")
}