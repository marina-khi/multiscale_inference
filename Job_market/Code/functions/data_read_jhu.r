#################################
#Analysis of covid data from JHU#
#################################
rm(list=ls())

library(tidyr)
library(multiscale)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/functions.R")

#Defining necessary constants
alpha    <- 0.05 #confidence level for application
sim_runs <- 5000

#We load government response index as well
gov_responces      <- read.csv("data/OxCGRT_latest.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
gov_responces$Date <- as.Date(as.character(gov_responces$Date), format = "%Y%m%d")
names(gov_responces)[names(gov_responces) == 'CountryCode'] <- 'countryterritoryCode'
names(gov_responces)[names(gov_responces) == 'Date']        <- 'dateRep'
gov_responces      <- gov_responces[gov_responces$RegionCode == "", ]

#Loading the world coronavirus data
covid_tmp          <- read.csv("data/time_series_covid19_confirmed_global.csv", sep = ",", 
                               stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)
#covid_tmp          <- covid_tmp[is.na(covid_tmp$`Province/State`), ]
names(covid_tmp)[names(covid_tmp) == "Country/Region"] <- 'CountryName'
covid_tmp          <- covid_tmp[, -c(1, 3, 4)]

new_covid          <- aggregate(. ~ CountryName, covid_tmp, sum)
t_len              <- ncol(new_covid) - 1
new_covid2         <- gather(new_covid, key = "dateRep", value = "cumcases", 2:(t_len + 1))
new_covid2$dateRep <- as.Date(new_covid2$dateRep, format = "%m/%d/%y")


covid <- merge(new_covid2, gov_responces, by  = c('CountryName', 'dateRep'), all.x = TRUE)

rm(covid_tmp, new_covid2, new_covid)
rm(gov_responces)

covid$cases           <- 0
covid$weekday         <- weekdays(covid$dateRep)

covid_list <- list()
for (country in unique(covid$CountryName)){
  cumcases_column <- covid[covid$CountryName == country, "cumcases"]
  gov_resp_column <- covid[covid$CountryName == country, "GovernmentResponseIndex"]
  time_range      <- length(gov_resp_column)
  covid[covid$CountryName == country, "cases"] <- c(0, cumcases_column[2:time_range] - cumcases_column[1:(time_range - 1)])
  tmp <- max(covid[covid$CountryName == country, "cumcases"])
  if (tmp >= 1000){
    tmp_df <- covid[(covid$CountryName == country & covid$cumcases >= 100),
                    c("dateRep", "cases", "cumcases", "weekday",
                      "GovernmentResponseIndex")]
    tmp_index <- match("Monday", tmp_df[, "weekday"])
    covid_list[[country]] <- tmp_df[tmp_index:nrow(tmp_df), ]
  }
}


#We are mainly interested in the european countries
#We can eliminate the countries we are not interesed in
#covid_list <- covid_list[names(covid_list) %in% c("BRA","CAN", "CHL", "CHN", "COL", "ECU", "EGY", "IDN", "IND", "IRN", "MEX", "PAK", "PER", "POL", "ROU", "TUR", "USA") == FALSE]
#Or we can just leave those that are interesting
covid_list <- covid_list[names(covid_list) %in% c("Germany", "France", "United Kingdom", "Spain", "Italy") == TRUE]


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

i = 1
for (country in countries) {
  covid_mat[, i] <- covid_list[[country]]$cases[1:min_dates]
  gov_resp[, i]  <- covid_list[[country]]$GovernmentResponseIndex[1:min_dates]
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
sum(covid_mat < 0)
covid_mat[covid_mat < 0] <- 0

#Plotting different countries on one plot
dev.new()
matplot(1:min_dates, covid_mat, type = 'l', lty = 1, col = 1:min_dates, xlab = 'Number of days since 100th case', ylab = 'Deaths')
legend(1, 30000, legend = names(covid_list), lty = 1, col = 1:min_dates, cex = 0.8)
dev.off()
