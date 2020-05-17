########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(multiscale)

#Defining necessary constants
set.seed(123)
alpha   <- 0.05 #confidence level for application
SimRuns <- 1000
bw      <- 0.05 #Bandwidth for calculating the residuals for \hat{\sigma}


#Loading the world coronavirus data

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
    #covid_list[[country]] <- covid[(covid$geoId == country & covid$cumcases >= 100 & covid$cumdeaths > 0), c("dateRep", "cases", "deaths", "cumcases", "cumdeaths")]
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
#dev.new()
#matplot(1:min_dates, covid_mat_min, type = 'l', lty = 1, col = 1:min_dates, xlab = 'Number of deaths since 100th case', ylab = 'Deaths')
#legend(1, 40000, legend = names(covid_list), lty = 1, col = 1:min_dates, cex = 0.8)
#dev.off()

#For now, we focus our attention on the "short" time series

Tlen          <- min_dates
N_ts          <- countries_num #Updating the number of time series because of dropped stations

#covid_mat_min <- scale(covid_mat_min, scale = FALSE)

u.grid <-  seq(from = 3.5/Tlen, to = 1, by = 3.5/Tlen)
h.grid <- seq(from = 3.5/Tlen, to = 1/4, by = 3.5/Tlen)
grid_matrix <- expand.grid(u.grid, h.grid)
deletions <- ((grid_matrix$Var1 - grid_matrix$Var2 < 0) | grid_matrix$Var1 + grid_matrix$Var2 > 1 )

grid <- construct_grid(Tlen, u.grid = u.grid, h.grid = h.grid, deletions = !deletions) 

source("functions/functions.R")
grid_points <- seq(from = 1/Tlen, to = 1, length.out = Tlen) #grid points for estimating

smoothed_curve           <- matrix(NA, ncol = N_ts, nrow = min_dates)
colnames(smoothed_curve) <- countries

sigma_vec <- rep(0, N_ts)

for (i in 1:N_ts){
  #plot(NA, ylab="", xlab = "", xlim = c(0,1), ylim = c(0, max(covid_mat_min[, i]) + 100), xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
  smoothed_curve[, i] <- mapply(nadaraya_watson_smoothing, grid_points, MoreArgs = list(covid_mat_min[, i], grid_points, bw))
  r_it <- (covid_mat_min[, i] - smoothed_curve[, i]) / sqrt(smoothed_curve[, i])
  #for deaths we have zero values for the smoothed curve, so we don't take them into account
  #r_it <- (covid_mat_min[smoothed_curve[, i] > 0, i] - smoothed_curve[smoothed_curve[, i] > 0, i]) / sqrt(smoothed_curve[smoothed_curve[, i] > 0, i])
  sigma_vec[i] <- sqrt(mean(r_it^2))
  #lines(grid_points, covid_mat_min[, i])
  #lines(grid_points, smoothed_curve[, i], col = 'red')
}

sigma_vec_same <- rep(sqrt(mean(sigma_vec^2)), N_ts)

sigma_mat <- matrix(0.0, nrow = nrow(grid$gset), ncol = N_ts)
sigma_mat_same <- matrix(0.0, nrow = nrow(grid$gset), ncol = N_ts)

for (i in 1:N_ts){
  for (k in 1:nrow(grid$gset)){
    smoothed_curve[, i] <- mapply(nadaraya_watson_smoothing, grid_points, MoreArgs = list(covid_mat_min[, i], grid_points, bw = grid$gset[k, 2]))
    r_it <- (covid_mat_min[, i] - smoothed_curve[, i]) / sqrt(smoothed_curve[, i])
    sigma_mat[k, i] <- sqrt(mean(r_it^2))
  }
}

for (k in 1:nrow(grid$gset)){
  sigma_mat_same[k,] <- rep(sqrt(mean(sigma_mat[k,]^2)), N_ts)  
}

result_same <- multiscale_test(data = covid_mat_min, sigma_vec = sigma_vec_same, N_ts = N_ts, grid = grid,
                          SimRuns = SimRuns, epidem = TRUE)
result_diff <- multiscale_test(data = covid_mat_min, sigma_vec = sigma_vec, N_ts = N_ts, grid = grid,
                               SimRuns = SimRuns, epidem = TRUE)

result_same_bw <- multiscale_test(data = covid_mat_min, sigma_vec = rep(1, N_ts), N_ts = N_ts, grid = grid,
                                  SimRuns = SimRuns, epidem = TRUE, bw_dependent = TRUE, sigma_mat = sigma_mat_same)

result_diff_bw <- multiscale_test(data = covid_mat_min, sigma_vec = rep(1, N_ts), N_ts = N_ts, grid = grid,
                                  SimRuns = SimRuns, epidem = TRUE, bw_dependent = TRUE, sigma_mat = sigma_mat)



l = 1
for (i in 1:(N_ts - 1)){
  for (j in (i + 1):N_ts){
    if (l %in% c(6, 10, 13, 18, 19, 20, 34, 38, 43)){
      pdf(paste0("plots/cases_bw_00", bw*1000, "_", countries[i], "_vs_", countries[j], ".pdf"), width=7, height=9, paper="special")
      
      produce_plots(results = result_same, data_i = covid_mat_min[, i], data_j = covid_mat_min[, j],
                    smoothed_i = smoothed_curve[, i], smoothed_j = smoothed_curve[, j],
                    country_i = countries[i], country_j = countries[j],
                    text_ = "Same sigmas for all the time trends")
        
      produce_plots(results = result_diff, data_i = covid_mat_min[, i], data_j = covid_mat_min[, j],
                    smoothed_i = smoothed_curve[, i], smoothed_j = smoothed_curve[, j],
                    country_i = countries[i], country_j = countries[j],
                    text_ = "Different sigmas for each time trend")

      produce_plots(results = result_same_bw, data_i = covid_mat_min[, i], data_j = covid_mat_min[, j],
                    smoothed_i = smoothed_curve[, i], smoothed_j = smoothed_curve[, j],
                    country_i = countries[i], country_j = countries[j],
                    text_ = "Same bandwidth-specific sigmas for all time trends")
      produce_plots(results = result_diff_bw, data_i = covid_mat_min[, i], data_j = covid_mat_min[, j],
                    smoothed_i = smoothed_curve[, i], smoothed_j = smoothed_curve[, j],
                    country_i = countries[i], country_j = countries[j],
                    text_ = "Different bandwidth-specific sigmas for each time trend")
      dev.off()
    }
    l = l + 1
  }
}

#dev.off()
#axis(1, at = grid_points[seq(5, 385, by = 20)], labels = monthly_temp$date[seq(5, 385, by = 20)])
