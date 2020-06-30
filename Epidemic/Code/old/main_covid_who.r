########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(multiscale)

#Defining necessary constants
alpha    <- 0.05 #confidence level for application
sim_runs <- 1000
bw       <- 0.05 #Bandwidth for calculating the residuals for \hat{\sigma}

#We load government response index as well
gov_responces  <- read.csv("data/OxCGRT_latest.csv", sep = ";", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
gov_responces$Date <- as.Date(as.character(gov_responces$Date), format = "%Y%m%d")
names(gov_responces)[names(gov_responces) == 'CountryName'] <- 'Country'
names(gov_responces)[names(gov_responces) == 'Date'] <- 'Date_reported'

#Additional cleaning so that the names of the countries are coded in the same way
gov_responces[gov_responces$Country == "United Kingdom", "Country"] <- "The United Kingdom"


#Loading the world coronavirus data

covid_tmp          <- read.csv("data/WHO-COVID-19-global-data.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "")
dtparts            <- t(as.data.frame(strsplit(covid_tmp$Date_reported,'T')))
row.names(dtparts) <- NULL

covid_tmp$Date_reported <- as.Date(dtparts[, 1], format = "%Y-%m-%d")
covid_tmp$WHO_region    <- NULL
covid_tmp               <- complete(covid_tmp, Date_reported = seq.Date(min(Date_reported), max(Date_reported), by='day'),
                                    Country, fill = list(New_cases = 0, New_deaths = 0))

covid <- merge(covid_tmp, gov_responces, by  = c('Country', "Date_reported"), all.x = TRUE)

rm(covid_tmp)
rm(gov_responces)

covid$cumcases        <- 0
covid$cumdeaths       <- 0
covid$lagged_gov_resp <- 0

covid_list <- list()
for (country in unique(covid$Country)){
  covid[covid$Country == country, "cumcases"] <- cumsum(covid[covid$Country == country, "New_cases"])
  covid[covid$Country == country, "cumdeaths"] <- cumsum(covid[covid$Country == country, "New_deaths"])
  gov_resp_column <- covid[covid$Country == country, "GovernmentResponseIndex"]
  time_range      <- length(gov_resp_column)
  covid[covid$Country == country, "lagged_gov_resp"] <- c(rep(NA, 14), gov_resp_column[1:(time_range - 14)])
  tmp <- max(covid[covid$Country == country, "cumdeaths"])
  if (tmp >= 1000){
    #We restrict our attention only to the contries with more than 1000 deaths and only starting from 100th case
    covid_list[[country]] <- covid[(covid$Country == country & covid$cumcases >= 100),
                                   c("Date_reported", "New_cases", "New_deaths", "cumcases", "cumdeaths", "GovernmentResponseIndex", "lagged_gov_resp")]
  }
}

#We are mainly interested in the european countries
#We can eliminate the countries we are not interesed in
#covid_list <- covid_list[names(covid_list) %in% c("BRA","CAN", "CHL", "CHN", "COL", "ECU", "EGY", "IDN", "IND", "IRN", "MEX", "PAK", "PER", "POL", "ROU", "TUR", "USA") == FALSE]
#Or we can just leave those that are interesting
covid_list <- covid_list[names(covid_list) %in% c("Germany", "France", "The United Kingdom", "Spain", "Italy") == TRUE]


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
  covid_mat[, i] <- covid_list[[country]]$New_cases[1:min_dates]
  gov_resp[, i]  <- covid_list[[country]]$GovernmentResponseIndex[1:min_dates]
  lagged_gov_resp[, i]  <- covid_list[[country]]$lagged_gov_resp[1:min_dates]
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
covid_mat[covid_mat < 0] <- abs(covid_mat[covid_mat < 0])

#Plotting different countries on one plot
#dev.new()
#matplot(1:min_dates, covid_mat, type = 'l', lty = 1, col = 1:min_dates, xlab = 'Number of deaths since 100th case', ylab = 'Deaths')
#legend(1, 40000, legend = names(covid_list), lty = 1, col = 1:min_dates, cex = 0.8)
#dev.off()

Tlen          <- min_dates
n_ts          <- countries_num #Updating the number of time series because of dropped stations
grid          <- construct_weekly_grid(Tlen) 

source("functions/functions.R")
grid_points <- seq(from = 1 / Tlen, to = 1, length.out = Tlen) #grid points for estimating

smoothed_curve           <- matrix(NA, ncol = n_ts, nrow = min_dates)
colnames(smoothed_curve) <- countries

# #This is the first method of estimating sigma from the variance
sigma_vec <- rep(0, n_ts)
for (i in 1:n_ts){
  smoothed_curve[, i] <- mapply(nadaraya_watson_smoothing, grid_points, MoreArgs = list(covid_mat[, i], grid_points, bw = 3.5 / Tlen))
  r_it <- (covid_mat[, i] - smoothed_curve[, i]) / sqrt(smoothed_curve[, i])
  #for deaths we have zero values for the smoothed curve, so we don't take them into account
  #r_it <- (covid_mat[smoothed_curve[, i] > 0, i] - smoothed_curve[smoothed_curve[, i] > 0, i]) / sqrt(smoothed_curve[smoothed_curve[, i] > 0, i])
  sigma_vec[i] <- sqrt(mean(r_it^2))
}

#This is the second method of estimating sigma from the variance
sigma_bwfree_vec <- rep(0, n_ts)

for (i in 1:n_ts){
  sigma_bwfree_squared <- sum((covid_mat[2:Tlen, i] - covid_mat[1:(Tlen - 1), i])^2) / (2 * sum(covid_mat[, i]))
  sigma_bwfree_vec[i] <- sqrt(sigma_bwfree_squared)
}

result_bwfree <- multiscale_test(data = covid_mat, sigma_vec = sigma_bwfree_vec,
                                 n_ts = n_ts, grid = grid,
                                 sim_runs = sim_runs, epidem = TRUE)

#Pairwise coparison DEU vs ESP, DEU vs FRA, DEU vs GBR and DEU vs ITA
#here we code all the pairs that we want to compare in a matrix
#ijset <- matrix(data = c(1, 2, 1, 3, 1, 4, 1, 5), ncol = 2, byrow = TRUE)
#
#result_some_countries <- multiscale_test(data = covid_mat, sigma_vec = sigma_vec,
#                          n_ts = n_ts, grid = grid, ijset = ijset,
#                          sim_runs = sim_runs, epidem = TRUE)
#
#result_bwfree_some_countries <- multiscale_test(data = covid_mat, sigma_vec = sigma_bwfree_vec,
#                                 n_ts = n_ts, grid = grid, ijset = ijset,
#                                 sim_runs = sim_runs, epidem = TRUE)
# 
# 
# 

#Plotting pairwise comparison for difference based sigma estimator
pdf(paste0("plots/bandwidth_free_sigma_all_countries_who.pdf"), width=7, height=9, paper="special")
for (l in seq_len(nrow(result_bwfree$ijset))){
  i <- result_bwfree$ijset[l, 1]
  j <- result_bwfree$ijset[l, 2]
  produce_plots(results = result_bwfree, l = l, data_i = covid_mat[, i], data_j = covid_mat[, j],
                smoothed_i = smoothed_curve[, i], smoothed_j = smoothed_curve[, j],
                gov_resp_i = gov_resp[, i], gov_resp_j = gov_resp[, j],
                lagged_gov_resp_i = lagged_gov_resp[, i], lagged_gov_resp_j = lagged_gov_resp[, j],
                country_i = countries[i], country_j = countries[j],
                text_ = "Same bandwidth-free sigmas for all the time trends")
}
dev.off()