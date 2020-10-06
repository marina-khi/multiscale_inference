#################################
#Analysis of covid data from WHO#
#################################
rm(list=ls())

library(tidyr)
library(multiscale)

#Defining necessary constants
alpha    <- 0.05 #confidence level for application
sim_runs <- 1000
source("functions/functions.R")

#We load government response index as well
gov_responces      <- read.csv("data/OxCGRT_latest.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
gov_responces$Date <- as.Date(as.character(gov_responces$Date), format = "%Y%m%d")
names(gov_responces)[names(gov_responces) == 'CountryName'] <- 'Country'
names(gov_responces)[names(gov_responces) == 'Date'] <- 'Date_reported'
gov_responces      <- gov_responces[gov_responces$RegionCode == "", ]

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
covid$weekday         <- weekdays(covid$Date_reported)

covid_list <- list()
for (country in unique(covid$Country)){
  covid[covid$Country == country, "cumcases"] <- cumsum(covid[covid$Country == country, "New_cases"])
  covid[covid$Country == country, "cumdeaths"] <- cumsum(covid[covid$Country == country, "New_deaths"])
  gov_resp_column <- covid[covid$Country == country, "GovernmentResponseIndex"]
  tmp <- max(covid[covid$Country == country, "cumdeaths"])
  if (tmp >= 1000){
    tmp_df <- covid[(covid$Country == country & covid$cumcases >= 100),
                    c("Date_reported", "New_cases", "New_deaths", "cumcases", "cumdeaths", "weekday",
                      "GovernmentResponseIndex")]
    tmp_index <- match("Monday", tmp_df[, "weekday"])
    covid_list[[country]] <- tmp_df[tmp_index:nrow(tmp_df), ]
    
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
dates         <- unique(covid$Date_reported)
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
  covid_mat[, i] <- covid_list[[country]]$New_cases[1:min_dates]
  gov_resp[, i]  <- covid_list[[country]]$GovernmentResponseIndex[1:min_dates]
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
sum(covid_mat < 0)
covid_mat[covid_mat < 0] <- 0

#Plotting different countries on one plot
matplot(1:min_dates, covid_mat, type = 'l', lty = 1, col = 1:min_dates, xlab = 'Number of days since 100th case', ylab = 'Deaths')
legend(1, 20000, legend = names(covid_list), lty = 1, col = 1:min_dates, cex = 0.8)