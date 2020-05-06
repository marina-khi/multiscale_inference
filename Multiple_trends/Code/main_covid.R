########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(multiscale)

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

covid_list <- covid_list[c('BE', 'CN', 'DE', 'ES', 'FR', 'IT', 'NL', 'RU', 'UK', 'US')]
#CN = CHINA

#Calculate the number of days that we have data for all countries.
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
  len <- length(covid_list[[country]]$cases)
  if (len < max_dates){
    covid_mat_max[, i] <- c(covid_list[[country]]$cases, rep(NA, max_dates - len))
  } else {
    covid_mat_max[, i] <- covid_list[[country]]$cases[1:max_dates]
  }
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
covid_mat_min[covid_mat_min < 0] <- abs(covid_mat_min[covid_mat_min < 0])

#Plotting different countries on one plot
dev.new()
matplot(1:min_dates, covid_mat_min, type = 'l', lty = 1, col = 1:min_dates, xlab = 'Number of deaths since 100th case', ylab = 'Deaths')
legend(1, 13, legend = names(covid_list), lty = 1, col = 1:min_dates, cex = 0.8)
#dev.off()

#For now, we focus our attention on the "short" time series

Tlen          <- min_dates
N_ts          <- countries_num #Updating the number of time series because of dropped stations

#covid_mat_min <- scale(covid_mat_min, scale = FALSE)


#################
#Order selection#
#################

#orders <- select_order(covid_mat_min, q = 5:20, r = 5:15)


#########################################################
#Different ways to calculate long-run variance estimator#
#########################################################

##Setting tuning parameters
#q     <- 10
#r.bar <- 10
#
#sigmahat_vector_order_1 <- c()
#for (i in 1:N_ts){
#  AR.struc                <- estimate_lrv(data = covid_mat_min[, i], q = q, r.bar = r.bar, p=1)
#  sigma_hat_i             <- sqrt(AR.struc$lrv)
#  sigmahat_vector_order_1 <- c(sigmahat_vector_order_1, sigma_hat_i)
#}
#
#sigmahat_vector_order_2 <- c()
#for (i in 1:N_ts){
#  AR.struc                <- estimate_lrv(data = covid_mat_min[, i], q = q, r.bar = r.bar, p=2)
#  sigma_hat_i             <- sqrt(AR.struc$lrv)
#  sigmahat_vector_order_2 <- c(sigmahat_vector_order_2, sigma_hat_i)
#}
#
#sigmahat_vector_iid  <- estimate_iid_sd(covid_mat_min, method = 'first')
#sigmahat_vector_iid2 <- estimate_iid_sd(covid_mat_min, method = 'second')


#Calculating the statistic for real data


#result_1    <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector_order_1, SimRuns = SimRuns, N_ts = N_ts)
#result_2    <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector_order_2, SimRuns = SimRuns, N_ts = N_ts)
#result_iid  <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector_iid, SimRuns = SimRuns, N_ts = N_ts)
#result_iid2 <- multiscale_testing(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vector_iid2, SimRuns = SimRuns, N_ts = N_ts)

#sigmahat_vec <- rep(354.9, N_ts)
#result <- multiscale_test(alpha = alpha, data = covid_mat_min, sigma_vec = sigmahat_vec, SimRuns = SimRuns, N_ts = N_ts, deriv_order = 0)
Tlen <- as.integer(Tlen)

grid <- construct_grid(Tlen, u.grid = seq(from = 1/Tlen, to = 1, by = 1/Tlen),
                         h.grid = seq(from = 5/Tlen, to = 1/4, by = 1/Tlen))
  
# Select (1-alpha) quantile of the multiscale statistic under the null
quantiles <- compute_quantiles(Tlen, grid, N_ts, SimRuns = SimRuns, epidem = TRUE)
  
probs    <- as.vector(quantiles[1,])
quant    <- as.vector(quantiles[2,])
  
if(sum(probs == (1-alpha)) == 0)
  pos <- which.min(abs(probs-(1-alpha)))
if(sum(probs == (1-alpha)) != 0)
  pos <- which.max(probs == (1-alpha)) 

quant <- quant[pos]
  
# Compute test results
gset                   <- grid$gset
N                      <- as.integer(dim(gset)[1])
gset_cpp               <- as.matrix(gset)
gset_cpp               <- as.vector(gset_cpp) 
storage.mode(gset_cpp) <- "double"

Psi_ij <- multiscale_stats_multiple(T = Tlen, N_ts = N_ts, data = covid_mat_min,
                                    gset = gset_cpp, N = N, sigma_vec = rep(1, N_ts), deriv_order = 0, epidem = TRUE)
    if (max(Psi_ij$stat) > quant) {
      cat("We reject H_0 with probability", alpha, ". Psihat_statistic = ", max(Psi_ij$stat), ". Number of pairwise rejections = ", 
          sum(Psi_ij$stat > quant), ". Gaussian quantile value = ", quant, "\n")
    } else {
      cat("We fail to reject H_0 with probability", alpha, ". Psihat_statistic = ", max(Psi_ij$stat),
          ". Gaussian quantile value = ", quant, "\n")
    }
    return(list(quant = quant, statistics = Psi_ij$stat))
  }