library(tidyr)

covid         <- read.csv("covid.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
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

#covid_list <- covid_list[c('BE', 'CN', 'DE', 'ES', 'FR', 'IT', 'NL', 'RU', 'UK', 'US')]
#CN = CHINA
covid_list <- covid_list[c('DE', 'IT')]

#Calculate the number of days that we have data for all countries.
min_dates <- min(sapply(covid_list[names(covid_list) != "CN"], NROW)) #We are not considering China as it has too long dataset

countries     <- names(covid_list)
dates         <- unique(covid$dateRep)
countries_num <- length(covid_list)

covid_mat <- matrix(NA, ncol = countries_num, nrow = min_dates)
colnames(covid_mat) <- countries

i = 1
for (country in countries) {
  covid_mat[, i] <- covid_list[[country]]$cases[1:min_dates] #You can put 'cases', 'cumcases', or 'cumdeaths' here
  i = i + 1
}

plot(covid_mat[,1],ylim=c(min(covid_mat),max(covid_mat)),type="l")
for(i in 2:dim(covid_mat)[2])
   lines(covid_mat[,i],col=i)  




