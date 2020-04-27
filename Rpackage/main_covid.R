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
SimRuns <- 10


####################################
#Loading the world coronavirus data#
####################################

covid         <- read.csv("data/covid.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
covid$dateRep <- as.Date(covid$dateRep, format = "%d/%m/%Y")
covid         <- complete(covid, dateRep = seq.Date(min(dateRep), max(dateRep), by='day'), geoId, fill = list(cases = 0, deaths = 0))

covid$cumcases  <- 0
covid$cumdeaths <- 0

covid_list <- list()
for (country in unique(covid$geoId)){
  covid[covid$geoId == country, "cumcases"] <- cumsum(covid[covid$geoId == country, "cases"])
  covid[covid$geoId == country, "cumdeaths"] <- cumsum(covid[covid$geoId == country, "deaths"])
  tmp <- max(covid[covid$geoId == country, "cumdeaths"])
  if (tmp >= 1000){
    #We restrict our attention only to the contries with more than 1000 deaths and only starting from 100th case
    covid_list[[country]] <- covid[(covid$geoId == country & covid$cumcases >= 100), c("dateRep", "cases", "deaths", "cumcases", "cumdeaths")]
  }
}


#We are taking BR, CA, TR out
covid_list[['BR']]<-NULL
covid_list[['CA']]<-NULL
covid_list[['TR']]<-NULL

#If we are taking logs, then we need to remove those countries that have 0 cases some of the days
#covid_list[['IR']]<-NULL
#covid_list[['CH']]<-NULL


min_dates <- min(sapply(covid_list[names(covid_list) != "CN"], NROW)) #We are not considering China as it has too long dataset
max_dates <- max(sapply(covid_list[names(covid_list) != "CN"], NROW))

countries     <- names(covid_list)
dates         <- unique(covid$dateRep)
countries_num <- length(covid_list)

covid_mat_min <- matrix(NA, ncol = countries_num, nrow = min_dates)
covid_mat_max <- matrix(NA, ncol = countries_num, nrow = max_dates)

colnames(covid_mat_min) <- countries
colnames(covid_mat_max) <- countries


i = 1
for (country in countries) {
  covid_mat_min[, i] <- covid_list[[country]]$cases[1:min_dates]
  len <- length(covid_list[[country]]$deaths)
  if (len < max_dates){
    covid_mat_max[, i] <- c(covid_list[[country]]$cases, rep(NA, max_dates - len))
  } else {
    covid_mat_max[, i] <- covid_list[[country]]$cases[1:max_dates]
  }
  i = i + 1
}

#Plotting different countries on one plot
#dev.new()
#matplot(1:min_dates, covid_mat_min, type = 'l', lty = 1, col = 1:min_dates, xlab = 'Number of deaths since 100th case', ylab = 'Deaths')
#legend(1, 13, legend = names(covid_list), lty = 1, col = 1:min_dates, cex = 0.8)
#dev.off()

#For now, we focus our attention on the "short" time series

Tlen          <- min_dates
N_ts          <- countries_num #Updating the number of time series because of dropped stations

covid_mat_min <- scale(covid_mat_min, scale = FALSE)


#################
#Order selection#
#################

source("functions/order_selection.r")
orders <- order_selection(covid_mat_min, q = 5:20, r = 5:15)


#Setting tuning parameters
q     <- 10
r.bar <- 10


#########################################################
#Different ways to calculate long-run variance estimator#
#########################################################

sigmahat_vector_order_1 <- c()
for (i in 1:N_ts){
  AR.struc                <- AR_lrv(data = covid_mat_min[, i], q = q, r.bar = r.bar, p=1)
  sigma_hat_i             <- sqrt(AR.struc$lrv)
  sigmahat_vector_order_1 <- c(sigmahat_vector_order_1, sigma_hat_i)
}

sigmahat_vector_order_2 <- c()
for (i in 1:N_ts){
  AR.struc                <- AR_lrv(data = covid_mat_min[, i], q = q, r.bar = r.bar, p=2)
  sigma_hat_i             <- sqrt(AR.struc$lrv)
  sigmahat_vector_order_2 <- c(sigmahat_vector_order_2, sigma_hat_i)
}

sigmahat_vector_iid  <- sigmahat_vec_iid(covid_mat_min)
sigmahat_vector_iid2 <- sigmahat_vec_iid2(covid_mat_min)


#Calculating the statistic for real data

result_1    <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector_order_1, SimRuns = SimRuns, N_ts = N_ts)
result_2    <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector_order_2, SimRuns = SimRuns, N_ts = N_ts)
result_iid  <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector_iid, SimRuns = SimRuns, N_ts = N_ts)
result_iid2 <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector_iid2, SimRuns = SimRuns, N_ts = N_ts)