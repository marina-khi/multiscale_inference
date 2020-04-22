# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

library(Rcpp)
library(tidyr)
library(tictoc)

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

alpha   <- 0.05 #confidence level for application
SimRuns <- 5000


###########################################
#Loading the real station data for England#
###########################################

covid <- read.csv("data/COVID-19.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
covid <- within(covid, rm("day", "month", "year", "countryterritoryCode", "countriesAndTerritories", "continentExp"))

covid$dateRep <- as.Date(covid$dateRep, format = "%d/%m/%Y")

covid <- complete(covid, dateRep = seq.Date(min(dateRep), max(dateRep), by='day'), geoId, fill = list(cases = 0, deaths = 0))

covid$cumcases  <- 0
covid$cumdeaths <- 0

countries <- unique(covid$geoId)
dates <- unique(covid$dateRep)

covid_list <- list()
#covid_df$dates <- as.Date(covid_df$dates, format = "%d.%m.%Y")
#covid_df <- covid_df[order(covid_df$dates), ]

for (country in countries){
  covid[covid$geoId == country, "cumcases"] <- cumsum(covid[covid$geoId == country, "cases"])
  covid[covid$geoId == country, "cumdeaths"] <- cumsum(covid[covid$geoId == country, "deaths"])
  tmp <- covid[covid$geoId == country, c("dateRep", "cases", "deaths", "cumcases", "cumdeaths")]
  if (max(tmp$cumdeaths) >= 1000){
    covid_list[[country]] <- covid[(covid$geoId == country & covid$cumcases >= 100), c("dateRep", "cases", "deaths", "cumcases", "cumdeaths")]
  }
}

rm(tmp)

min_dates <- min(sapply(covid_list[names(covid_list) != "CN"], NROW))
max_dates <- max(sapply(covid_list[names(covid_list) != "CN"], NROW))

countries_num <- length(covid_list)

covid_mat_min <- matrix(NA, ncol = countries_num, nrow = min_dates)
covid_mat_max <- matrix(NA, ncol = countries_num, nrow = max_dates)
colnames(covid_mat_min) <- names(covid_list)
colnames(covid_mat_max) <- names(covid_list)


i = 1
for (country in names(covid_list)) {
  covid_mat_min[, i] <- log(covid_list[[country]]$cumcases[1:min_dates])
  len <- length(covid_list[[country]]$deaths)
  if (len < max_dates){
    covid_mat_max[, i] <- c(log(covid_list[[country]]$cumcases), rep(0, max_dates - len))
  } else {
    covid_mat_max[, i] <- log(covid_list[[country]]$cumcases[1:max_dates])
  }
  i = i + 1
}

dev.new()
matplot(1:min_dates, covid_mat_min, type = 'l', lty = 1, col = 1:min_dates, xlab = 'Number of deaths since 100th case', ylab = 'Deaths')
legend(1, 1200, legend = names(covid_list), lty = 1, col = 1:min_dates, cex = 0.8)

#dev.off()

Tlen          <- min_dates
N_ts          <- countries_num #Updating the number of time series because of dropped stations

covid_mat_min <- scale(covid_mat_min, scale = FALSE)


#####################
#Estimating variance#
#####################

#Order selection
q <- 5:15
r <- 5:10
order_results <- c()

for (j in 1:N_ts){
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
      AR.struc      <- AR_lrv(data=covid_mat_min[, j], q=criterion_matrix$q[[i]], r.bar=criterion_matrix$r[[i]], p=order)
      sigma_eta_hat <- sqrt(AR.struc$vareta)
      FPE <- c(FPE, (sigma_eta_hat^2 * (Tlen + order)) / (Tlen - order))
      AIC <- c(AIC, Tlen * log(sigma_eta_hat^2) + 2 * order)
      AICC <- c(AICC, Tlen * log(sigma_eta_hat^2) + Tlen * (1 + order / Tlen)/(1 - (order +2)/Tlen))
      SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(Tlen) / Tlen)
      HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(Tlen)) / Tlen)
    }
    criterion_matrix$FPE[[i]]  <- which.min(FPE)
    criterion_matrix$AIC[[i]]  <- which.min(AIC)
    criterion_matrix$AICC[[i]] <- which.min(AICC)
    criterion_matrix$SIC[[i]]  <- which.min(SIC)
    criterion_matrix$HQ[[i]]   <- which.min(HQ)
  }
  maxim <- max(criterion_matrix[, 3:7])
  order_results <- c(order_results, maxim)
  cat("For country ", colnames(covid_mat_min)[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
}



#Setting tuning parameters for testing
q     <- 10
r.bar <- 10


#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 1:N_ts){
  AR.struc        <- AR_lrv(data = covid_mat_min[, i], q = q, r.bar = r.bar, p=1)
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}



#Calculating the statistic for real data

result <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector, SimRuns = SimRuns, N_ts = N_ts)

#And now the testing itself
if (max(result$Psi_ij) > result$quant) {
  cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", max(result$Psi_ij),
      "Gaussian quantile value = ", result$quant, "\n")
} else {
  cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", max(result$Psi_ij),
      "Gaussian quantile value = ", result$quant, "\n")
}