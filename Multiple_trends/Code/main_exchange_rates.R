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


#################################
#Loading the exchange rates data#
#################################

exchange_rates <- read.csv("data/exchange_rates.csv", sep = ",", dec = ",", quote = '"', stringsAsFactors = FALSE)
exchange_rates <- within(exchange_rates, rm("exusal", "exusec", "exuseu", "exusir", "exusnz", "exusuk", "twexb", "twexm", "twexo", "indexgx", "exvzus"))

#exchange_rates <- as.matrix(exchange_rates)
colSums(is.na(exchange_rates))

exchange_rates <- exchange_rates[,colSums(is.na(exchange_rates)) <= 100] #Ommitting the time series with too sparse data
exchange_rates <- na.omit(exchange_rates)#Deleting the rows with ommitted variables

t_len         <- nrow(exchange_rates)
n_ts          <- ncol(exchange_rates) - 1 #Updating the number of time series because of dropped stations

#column_names <- names(exchange_rates)

#exchange_rates <-matrix(unlist(exchange_rates), nrow = t_len)

#colnames(exchange_rates) <- column_names
exchange_rates[, 2:(n_ts + 1)] <- scale(exchange_rates[, 2:(n_ts + 1)], scale = FALSE)


#####################
#Estimating variance#
#####################

#Order selection
q <- 30:60
r <- 10:15
order_results <- c()

for (j in 2:(n_ts + 1)){
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
      AR.struc      <- estimate_lrv(data=exchange_rates[[j]], q=criterion_matrix$q[[i]], r_bar=criterion_matrix$r[[i]], p=order)
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
  maxim <- max(criterion_matrix[, 3:7])
  order_results <- c(order_results, maxim)
  cat("For stock ", names(exchange_rates)[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
}


#Setting tuning parameters for testing
q     <- 55
r_bar <- 10


#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 2:(n_ts+1)){
  AR.struc        <- estimate_lrv(data = exchange_rates[[i]], q = q, r_bar = r_bar, p=order_results[i-1])
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}

#Constructing the grid
grid <- construct_grid(t = t_len)

#Calculating the statistic for real data
result <- multiscale_test(data = matrix(unlist(exchange_rates[, 2:(n_ts + 1)]), ncol = n_ts, byrow = FALSE),
                          sigma_vec = sigmahat_vector,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

