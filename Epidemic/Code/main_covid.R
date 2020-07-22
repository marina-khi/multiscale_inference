########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(multiscale)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/functions.R")

#Defining necessary constants
alpha    <- 0.05 #confidence level for application
sim_runs <- 5000 #Number of simulation runs to produce the Gaussian qauntiles

#Loading the world coronavirus data
covid_tmp         <- read.csv("data/covid.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "")
covid_tmp         <- covid_tmp[!is.na(covid_tmp$countryterritoryCode), ]
covid_tmp$dateRep <- as.Date(covid_tmp$dateRep, format = "%d/%m/%Y")
covid_tmp         <- complete(covid_tmp, dateRep = seq.Date(min(dateRep), max(dateRep), by='day'),
                              countryterritoryCode, fill = list(cases = 0, deaths = 0))

#We load government response index as well
gov_responces      <- read.csv("data/OxCGRT_latest.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
gov_responces$Date <- as.Date(as.character(gov_responces$Date), format = "%Y%m%d")
names(gov_responces)[names(gov_responces) == 'CountryCode'] <- 'countryterritoryCode'
names(gov_responces)[names(gov_responces) == 'Date']        <- 'dateRep'

#Merging the two datasets
covid <- merge(covid_tmp, gov_responces, by  = c('countryterritoryCode', 'dateRep'), all.x = TRUE)
rm(covid_tmp)
rm(gov_responces)

#Now we "normalize" the data counting only the countries with more than 1000 deaths overall
#and taking the day of 100th case as the starting point
covid$cumcases        <- 0
covid$cumdeaths       <- 0
covid$lagged_gov_resp <- 0

covid_list <- list()
for (country in unique(covid$countryterritoryCode)){
  covid[covid$countryterritoryCode == country, "cumcases"] <- cumsum(covid[covid$countryterritoryCode == country, "cases"])
  covid[covid$countryterritoryCode == country, "cumdeaths"] <- cumsum(covid[covid$countryterritoryCode == country, "deaths"])
  tmp <- max(covid[covid$countryterritoryCode == country, "cumdeaths"])
  if (tmp >= 1000){
    covid_list[[country]] <- covid[(covid$countryterritoryCode == country & covid$cumcases >= 100),
                                   c("dateRep", "cases", "deaths", "cumcases", "cumdeaths",
                                     "GovernmentResponseIndex")]
  }
}

#We are mainly interested in the "main" european countries,
#so we leave all the others out of the analysis
covid_list <- covid_list[names(covid_list) %in% c("DEU", "FRA", "GBR", "ESP", "ITA") == TRUE]


#Calculate the number of days that we have data for all fivecountries.
#We are not considering CHN = China as it has too long dataset.
t_len     <- min(sapply(covid_list[names(covid_list) != "CHN"], NROW))
countries <- names(covid_list)
dates     <- unique(covid$dateRep)
n_ts      <- length(covid_list) #Number of time series

#In order not to work with lists, we construct two matrices,
#one with number of cases for all countries and 
#one with Government Responce Index for all countries.
#It is not necessary, but it is more convenient to work with.
covid_mat           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(covid_mat) <- countries

gov_resp            <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(gov_resp)  <- countries

i = 1
for (country in countries) {
  covid_mat[, i]        <- covid_list[[country]]$cases[1:t_len]
  gov_resp[, i]         <- covid_list[[country]]$GovernmentResponseIndex[1:t_len]
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
sum(covid_mat < 0)
covid_mat[covid_mat < 0] <- 0

#We construct the grid (family of intervals) for our test statistics
grid <- construct_weekly_grid(t_len) 

#And we plot it
all_intervals <- data.frame('startpoint' = (grid$gset$u - grid$gset$h) * t_len,
                            'endpoint' = (grid$gset$u + grid$gset$h) * t_len,
                            'values' = 0)
all_intervals$values <- (1:nrow(all_intervals)) / nrow(all_intervals)

pdf("plots/all_intervals.pdf", width=5, height=5, paper="special")
par(mar = c(3, 0.5, 2, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
plot(NA, xlim=c(0,t_len),  ylim = c(0, 1 + 1/nrow(all_intervals)), xlab="days", ylab = "", yaxt= "n", mgp=c(2,0.5,0))
title(main = expression(The ~ family ~ of ~ intervals ~ italic(F)), line = 1)
segments(all_intervals$startpoint, all_intervals$values, all_intervals$endpoint, all_intervals$values, lwd = 2)
dev.off()

#We also need to estimate the overdispersion parameter
sigma_vec <- rep(0, n_ts)
for (i in 1:n_ts){
  sigma_squared <- sum((covid_mat[2:t_len, i] - covid_mat[1:(t_len - 1), i])^2) / (2 * sum(covid_mat[, i]))
  sigma_vec[i]  <- sqrt(sigma_squared)
}
sigmahat <- sqrt(mean(sigma_vec * sigma_vec))

#Now we are ready to perform the test.
result <- multiscale_test(data = covid_mat, sigma = sigmahat,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs)

#Rename the countries for the plot
countries_names <- c("Germany", "Spain", "France", "United Kingdom", "Italy")

#Plotting pairwise comparison 
for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  filename = paste0("plots/", countries[i], "_vs_", countries[j], ".pdf")
  produce_plots(results = result, l = l, data_i = covid_mat[, i], data_j = covid_mat[, j],
                gov_resp_i = gov_resp[, i], gov_resp_j = gov_resp[, j],
                country_i = countries_names[i], country_j = countries_names[j],
                filename = filename)
}
