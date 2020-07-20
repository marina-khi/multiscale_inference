########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(multiscale)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


#Defining necessary constants
alpha    <- 0.05 #confidence level for application
sim_runs <- 5000
bw       <- 0.05 #Bandwidth for calculating the residuals for \hat{\sigma}

#We load government response index as well
gov_responces  <- read.csv("data/OxCGRT_latest.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
gov_responces$Date <- as.Date(as.character(gov_responces$Date), format = "%Y%m%d")
names(gov_responces)[names(gov_responces) == 'CountryCode'] <- 'countryterritoryCode'
names(gov_responces)[names(gov_responces) == 'Date'] <- 'dateRep'

#Loading the world coronavirus data

covid_tmp         <- read.csv("data/covid.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "")
covid_tmp         <- covid_tmp[!is.na(covid_tmp$countryterritoryCode), ]
covid_tmp$dateRep <- as.Date(covid_tmp$dateRep, format = "%d/%m/%Y")
covid_tmp         <- complete(covid_tmp, dateRep = seq.Date(min(dateRep), max(dateRep), by='day'),
                              countryterritoryCode, fill = list(cases = 0, deaths = 0))

covid <- merge(covid_tmp, gov_responces, by  = c('countryterritoryCode', 'dateRep'), all.x = TRUE)

rm(covid_tmp)
rm(gov_responces)

covid$cumcases        <- 0
covid$cumdeaths       <- 0
covid$lagged_gov_resp <- 0

covid_list <- list()
for (country in unique(covid$countryterritoryCode)){
  covid[covid$countryterritoryCode == country, "cumcases"] <- cumsum(covid[covid$countryterritoryCode == country, "cases"])
  covid[covid$countryterritoryCode == country, "cumdeaths"] <- cumsum(covid[covid$countryterritoryCode == country, "deaths"])
  gov_resp_column <- covid[covid$countryterritoryCode == country, "GovernmentResponseIndex"]
  #gov_resp_column <- covid[covid$countryterritoryCode == country, "StringencyIndex"]
  time_range      <- length(gov_resp_column)
  covid[covid$countryterritoryCode == country, "lagged_gov_resp"] <- c(rep(NA, 14), gov_resp_column[1:(time_range - 14)])
  
  tmp <- max(covid[covid$countryterritoryCode == country, "cumdeaths"])
  if (tmp >= 1000){
    #We restrict our attention only to the contries with more than 1000 deaths and only starting from 100th case
    covid_list[[country]] <- covid[(covid$countryterritoryCode == country & covid$cumcases >= 100),
                                   c("dateRep", "cases", "deaths", "cumcases", "cumdeaths", "GovernmentResponseIndex", "lagged_gov_resp")]
  }
}

#We are mainly interested in the european countries
#We can eliminate the countries we are not interesed in
#covid_list <- covid_list[names(covid_list) %in% c("BRA","CAN", "CHL", "CHN", "COL", "ECU", "EGY", "IDN", "IND", "IRN", "MEX", "PAK", "PER", "POL", "ROU", "TUR", "USA") == FALSE]
#Or we can just leave those that are interesting
covid_list <- covid_list[names(covid_list) %in% c("DEU", "FRA", "GBR", "ESP", "ITA") == TRUE]


#Calculate the number of days that we have data for all countries.
#We are not considering CHN = China as it has too long dataset
min_dates <- min(sapply(covid_list[names(covid_list) != "CHN"], NROW))

countries     <- names(covid_list)
dates         <- unique(covid$dateRep)
countries_num <- length(covid_list)

#In order not to work with lists, we construct three matrices,
#one with number of cases for all countries, one with Government Responce Index for all countries
#and one with lagged Government responce Index for all countries.
#It is not necessary, but it is more convenient to work with.
covid_mat           <- matrix(NA, ncol = countries_num, nrow = min_dates)
colnames(covid_mat) <- countries

gov_resp            <- matrix(NA, ncol = countries_num, nrow = min_dates)
colnames(gov_resp) <- countries

lagged_gov_resp           <- matrix(NA, ncol = countries_num, nrow = min_dates)
colnames(lagged_gov_resp) <- countries


i = 1
for (country in countries) {
  covid_mat[, i] <- covid_list[[country]]$cases[1:min_dates]
  gov_resp[, i]  <- covid_list[[country]]$GovernmentResponseIndex[1:min_dates]
  lagged_gov_resp[, i]  <- covid_list[[country]]$lagged_gov_resp[1:min_dates]
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
covid_mat[covid_mat < 0] <- 0

Tlen          <- min_dates
n_ts          <- countries_num #Updating the number of time series because of dropped stations
grid          <- construct_weekly_grid(Tlen) 

#We plot the family of the considered intervals for the paper
all_intervals <- data.frame('startpoint' = (grid$gset$u - grid$gset$h) * Tlen,
                            'endpoint' = (grid$gset$u + grid$gset$h) * Tlen, 'values' = 0)
all_intervals$values <- (1:nrow(all_intervals)) / nrow(all_intervals)

pdf("plots/all_intervals.pdf", width=5, height=5, paper="special")
par(mar = c(3, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
plot(NA, xlim=c(0,Tlen),  ylim = c(0, 1 + 1/nrow(all_intervals)), xlab="days", ylab = "", yaxt= "n", mgp=c(2,0.5,0))
title(main = expression(The ~ family ~ of ~ intervals ~ italic(F)), line = 1)
segments(all_intervals$startpoint, all_intervals$values, all_intervals$endpoint, all_intervals$values, lwd = 2)
dev.off()


source("functions/functions.R")
grid_points <- seq(from = 1 / Tlen, to = 1, length.out = Tlen) #grid points for estimating

smoothed_curve           <- matrix(NA, ncol = n_ts, nrow = min_dates)
colnames(smoothed_curve) <- countries

# #This is the first method of estimating sigma from the variance
for (i in 1:n_ts){
  smoothed_curve[, i] <- mapply(nadaraya_watson_smoothing, grid_points,
                                MoreArgs = list(covid_mat[, i], grid_points, bw = 3.5 / Tlen))
  #r_it <- (covid_mat[, i] - smoothed_curve[, i]) / sqrt(smoothed_curve[, i])
  #for deaths we have zero values for the smoothed curve, so we don't take them into account
  #r_it <- (covid_mat[smoothed_curve[, i] > 0, i] - smoothed_curve[smoothed_curve[, i] > 0, i]) / sqrt(smoothed_curve[smoothed_curve[, i] > 0, i])
  #sigma_vec[i] <- sqrt(mean(r_it^2))
}

sigma_vec <- rep(0, n_ts)
#We estimate the overdispersion parameter:
for (i in 1:n_ts){
  sigma_squared <- sum((covid_mat[2:Tlen, i] - covid_mat[1:(Tlen - 1), i])^2) / (2 * sum(covid_mat[, i]))
  sigma_vec[i] <- sqrt(sigma_squared)
}

sigmahat <- mean(sigma_vec * sigma_vec)

result_bwfree <- multiscale_test(data = covid_mat, sigma = sigmahat,
                                n_ts = n_ts, grid = grid,
                                sim_runs = sim_runs, epidem = TRUE)

#THIS IS VERY SPECIFIC, WE RENAME THE COUNTRIES FOR THE PLOT
countries_names <- c("Germany", "Spain", "France", "United Kingdom", "Italy")

#Plotting pairwise comparison for difference based sigma estimator
for (l in seq_len(nrow(result_bwfree$ijset))){
  i <- result_bwfree$ijset[l, 1]
  j <- result_bwfree$ijset[l, 2]
  filename = paste0("plots/", countries[i], "_vs_", countries[j], ".pdf")
  produce_plots(results = result_bwfree, l = l, data_i = covid_mat[, i], data_j = covid_mat[, j],
                smoothed_i = smoothed_curve[, i], smoothed_j = smoothed_curve[, j],
                gov_resp_i = gov_resp[, i], gov_resp_j = gov_resp[, j],
                lagged_gov_resp_i = lagged_gov_resp[, i], lagged_gov_resp_j = lagged_gov_resp[, j],
                country_i = countries_names[i], country_j = countries_names[j],
                filename = filename)
}
