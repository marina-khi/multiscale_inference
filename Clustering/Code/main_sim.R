########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(multiscale)
library(xtable)
library(aweek)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

library(dendextend)

source("functions.R")

#Defining necessary constants
b_bar  <- 2
bw_abs <- 7
t_len  <- 200
n_ts   <- 6
sigma  <- 15

# functions for data simulations
lambda_fct <- function(u, c = 1000, height = 5000, position = 10) {
  return (height * exp(-(position * u - 3) ^ 2 / 2) + c)
}

r_doublepois <- function(n, mu, theta) {
  rnbinom(n = n, mu = mu, size = mu/(theta - 1))
}

simulate_data <- function(n_ts, t_len, lambda_vec, sigma) {
  data <- matrix(0, ncol = n_ts, nrow = t_len)
  for(t in 1:t_len) {
    data[t, ] <- r_doublepois(n = n_ts, lambda_vec[t], sigma^2)
  }
  return(data)
}


lambda_vec_1 <- lambda_fct((1:t_len) / t_len, c = 1000, height = 6000, position = 10)
lambda_vec   <- lambda_fct((1:t_len) / t_len, c = 1000, height = 5000, position = 10)

par(mar = c(3, 2, 2, 0)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins
plot((1:t_len) / t_len, lambda_vec,  ylim = c(0, max(lambda_vec) + 100), xlab="u",
     ylab = "", mgp=c(2,0.5,0), type = "l")
title(main = expression(Plot ~ of ~ the ~ "function" ~ lambda), line = 1)

Y1 <- simulate_data(n_ts = 1, t_len = t_len, lambda_vec = lambda_vec_1, sigma = sigma)
Y2 <- simulate_data(n_ts = n_ts - 1, t_len = t_len, lambda_vec = lambda_vec_2, sigma = sigma)
Y  <- cbind(Y1, Y2)



b_grid <- seq(1, b_bar, by = 0.01)



m_hat <- function(vect_u, b, data_p, grid_p, bw){
  m_hat_vec <- c()
  for (u in vect_u){
    result = sum((abs((grid_p - u * b) / bw) <= 1) * data_p)
    norm = sum((abs((grid_p - u * b) / bw) <= 1))
    m_hat_vec <- c(m_hat_vec, result/norm)
  }
  return(m_hat_vec)
}

grid_points <- seq(1/t_len, 1, by = 1/t_len)
# integral_points <- seq(1/t_len, 1, by = 0.01/t_len)
# 
# m_hat_vec <- m_hat(integral_points, b = 1, covid_mat[, 1], grid_points, bw = bw_abs/t_len)
# plot(grid_points, covid_mat[, 1], type = "l")
# lines(integral_points, m_hat_vec, col = "red")
# 
# m_hat_vec <- m_hat(integral_points, b = 1, covid_mat[, 2], grid_points, bw = bw_abs/t_len)
# plot(grid_points, covid_mat[, 2], type = "l")
# lines(integral_points, m_hat_vec, col = "blue")
# 
# m_hat_vec <- m_hat(integral_points, b = 1, covid_mat[, 3], grid_points, bw = bw_abs/t_len)
# plot(grid_points, covid_mat[, 3], type = "l")
# lines(integral_points, m_hat_vec, col = "green")


integrand <- function(vect_u, b, data_points_i, data_points_j, 
                      norm_b, norm, grid_points, bw) {
  tmp <- m_hat(vect_u, b, data_points_i, grid_points, bw)/(norm_b * 1/b) - m_hat(vect_u, b = 1, data_points_j, grid_points, bw) / (norm * 1/b)
  return(tmp^2)
}

Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

for (b in b_grid){
  norm_b <- c()
  norm   <- c()
  for (k in 1:n_ts){
    norm_b <- c(norm_b, integrate(m_hat, lower = 0, upper = 1/b, b = b,
                                  data_p = covid_mat[, k], grid_p = grid_points,
                                  bw = bw_abs/t_len, subdivisions=2000)$value)
    norm <- c(norm, integrate(m_hat, lower = 0, upper = 1/b, b = 1,
                              data_p = covid_mat[, k], grid_p = grid_points, 
                              bw = bw_abs/t_len, subdivisions=2000)$value)
  }
  for (i in 1:(n_ts - 1)){
    for (j in (i + 1):n_ts){
      delta_ij <- 1/b * integrate(integrand, lower = 0, upper = 1/b, b = b,
                                  data_points_i = covid_mat[, i],
                                  data_points_j = covid_mat[, j],
                                  norm_b = norm_b[i], norm = norm[j],
                                  grid_points = grid_points, bw = bw_abs/t_len,
                                  subdivisions=2000)$value
      delta_ji <- 1/b * integrate(integrand, lower = 0, upper = 1/b, b = b,
                                  data_points_i = covid_mat[, j],
                                  data_points_j = covid_mat[, i],
                                  norm_b = norm_b[j], norm = norm[i],
                                  grid_points = grid_points, bw = bw_abs/t_len,
                                  subdivisions=2000)$value
      if (b == 1) {
        Delta_hat[i, j]   <- min(delta_ij, delta_ji)
      } else {
        if (min(delta_ij, delta_ji) < Delta_hat[i, j]) {
          Delta_hat[i, j] <- min(delta_ij, delta_ji)
        }
      }
      Delta_hat[j, i] <- Delta_hat[i, j]
      #cat("b = ", b, ", Delta_hat = ", Delta_hat[i, j], "\n")
    }
  }  
}

colnames(Delta_hat) <- countries
rownames(Delta_hat) <- countries

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
plot(res)########################
#Analysis of covid data#
########################
rm(list=ls())

library(tidyr)
library(multiscale)
library(xtable)
library(aweek)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

library(dendextend)

source("functions.R")

#Defining necessary constants
b_bar <- 2
bw_abs <- 7

#Loading the world coronavirus data
covid         <- read.csv("data/covid.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "")
covid         <- covid[!is.na(covid$countryterritoryCode), ]
covid$dateRep <- as.Date(covid$dateRep, format = "%d/%m/%Y")
covid         <- complete(covid, dateRep = seq.Date(min(dateRep), max(dateRep), by='day'),
                          countryterritoryCode, fill = list(cases = 0, deaths = 0))

#Now we "normalize" the data counting only the countries with more than 1000 deaths overall
#and taking the day of 100th case as the starting point
covid$weekday         <- weekdays(covid$dateRep)
covid$cumcases        <- 0
covid$cumdeaths       <- 0

covid_list <- list()
for (country in unique(covid$countryterritoryCode)){
  covid[covid$countryterritoryCode == country, "cumcases"]  <- cumsum(covid[covid$countryterritoryCode == country, "cases"])
  covid[covid$countryterritoryCode == country, "cumdeaths"] <- cumsum(covid[covid$countryterritoryCode == country, "deaths"])
  tmp <- max(covid[covid$countryterritoryCode == country, "cumdeaths"])
  if (tmp >= 1000){
    tmp_df <- covid[(covid$countryterritoryCode == country & covid$cumcases >= 100),
                    c("dateRep", "cases", "deaths", "cumcases", "cumdeaths", "weekday")]
    tmp_index <- match("Monday", tmp_df$weekday)
    #tmp_index = 1 #If we do not want to normalize by Mondays
    covid_list[[country]] <- tmp_df[tmp_index:nrow(tmp_df), ]
  }
}


#Calculate the number of days that we have data for all fivecountries.
#We are not considering CHN = China as it has too long dataset.
t_len     <- min(sapply(covid_list[names(covid_list) != "CHN"], NROW))
#t_len     <- 150 #We consider the first five months of the pandemic
countries <- names(covid_list)
dates     <- unique(covid$dateRep)
n_ts      <- length(covid_list) #Number of time series

#In order not to work with lists, we construct a matrix
#with number of cases for all countries and.
#It is not necessary, but it is more convenient to work with.
covid_mat           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(covid_mat) <- countries

i = 1
for (country in countries) {
  covid_mat[, i]        <- covid_list[[country]]$cases[1:t_len]
  i = i + 1
}

#Cleaning the data: there are weird cases in the dataset when the number of new cases is negative! 
sum(covid_mat < 0)
covid_mat[covid_mat < 0] <- 0

b_grid <- seq(1, b_bar, by = 0.01)



m_hat <- function(vect_u, b, data_p, grid_p, bw){
  m_hat_vec <- c()
  for (u in vect_u){
    result = sum((abs((grid_p - u * b) / bw) <= 1) * data_p)
    norm = sum((abs((grid_p - u * b) / bw) <= 1))
    m_hat_vec <- c(m_hat_vec, result/norm)
  }
  return(m_hat_vec)
}

grid_points <- seq(1/t_len, 1, by = 1/t_len)
# integral_points <- seq(1/t_len, 1, by = 0.01/t_len)
# 
# m_hat_vec <- m_hat(integral_points, b = 1, covid_mat[, 1], grid_points, bw = bw_abs/t_len)
# plot(grid_points, covid_mat[, 1], type = "l")
# lines(integral_points, m_hat_vec, col = "red")
# 
# m_hat_vec <- m_hat(integral_points, b = 1, covid_mat[, 2], grid_points, bw = bw_abs/t_len)
# plot(grid_points, covid_mat[, 2], type = "l")
# lines(integral_points, m_hat_vec, col = "blue")
# 
# m_hat_vec <- m_hat(integral_points, b = 1, covid_mat[, 3], grid_points, bw = bw_abs/t_len)
# plot(grid_points, covid_mat[, 3], type = "l")
# lines(integral_points, m_hat_vec, col = "green")


integrand <- function(vect_u, b, data_points_i, data_points_j, 
                      norm_b, norm, grid_points, bw) {
  tmp <- m_hat(vect_u, b, data_points_i, grid_points, bw)/(norm_b * 1/b) - m_hat(vect_u, b = 1, data_points_j, grid_points, bw) / (norm * 1/b)
  return(tmp^2)
}

Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

for (b in b_grid){
  norm_b <- c()
  norm   <- c()
  for (k in 1:n_ts){
    norm_b <- c(norm_b, integrate(m_hat, lower = 0, upper = 1/b, b = b,
                                  data_p = covid_mat[, k], grid_p = grid_points,
                                  bw = bw_abs/t_len, subdivisions=2000)$value)
    norm <- c(norm, integrate(m_hat, lower = 0, upper = 1/b, b = 1,
                              data_p = covid_mat[, k], grid_p = grid_points, 
                              bw = bw_abs/t_len, subdivisions=2000)$value)
  }
  for (i in 1:(n_ts - 1)){
    for (j in (i + 1):n_ts){
      delta_ij <- 1/b * integrate(integrand, lower = 0, upper = 1/b, b = b,
                                  data_points_i = covid_mat[, i],
                                  data_points_j = covid_mat[, j],
                                  norm_b = norm_b[i], norm = norm[j],
                                  grid_points = grid_points, bw = bw_abs/t_len,
                                  subdivisions=2000)$value
      delta_ji <- 1/b * integrate(integrand, lower = 0, upper = 1/b, b = b,
                                  data_points_i = covid_mat[, j],
                                  data_points_j = covid_mat[, i],
                                  norm_b = norm_b[j], norm = norm[i],
                                  grid_points = grid_points, bw = bw_abs/t_len,
                                  subdivisions=2000)$value
      if (b == 1) {
        Delta_hat[i, j]   <- min(delta_ij, delta_ji)
      } else {
        if (min(delta_ij, delta_ji) < Delta_hat[i, j]) {
          Delta_hat[i, j] <- min(delta_ij, delta_ji)
        }
      }
      Delta_hat[j, i] <- Delta_hat[i, j]
      #cat("b = ", b, ", Delta_hat = ", Delta_hat[i, j], "\n")
    }
  }  
}

colnames(Delta_hat) <- countries
rownames(Delta_hat) <- countries

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)
plot(res)