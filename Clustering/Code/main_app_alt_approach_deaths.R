###############################################
#Analysis of covid data - alternative approach#
###############################################
rm(list=ls())

library(tidyr)
library(aweek)
library(dendextend)
require(rworldmap)
source("functions.R")

#Defining necessary constants
bw_abs <- 7

#Loading the world coronavirus data
covid_deaths_tmp <- read.csv("data/time_series_covid19_deaths_global.csv", sep = ",", 
                             stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)
covid_cases_tmp  <- read.csv("data/time_series_covid19_confirmed_global.csv", sep = ",", 
                             stringsAsFactors = FALSE, na.strings = "", check.names = FALSE)

names(covid_deaths_tmp)[names(covid_deaths_tmp) == "Country/Region"] <- 'CountryName'
covid_deaths_tmp <- covid_deaths_tmp[, -c(1, 3, 4)]
names(covid_cases_tmp)[names(covid_cases_tmp) == "Country/Region"] <- 'CountryName'
covid_cases_tmp  <- covid_cases_tmp[, -c(1, 3, 4)]


new_covid_cases  <- aggregate(. ~ CountryName, covid_cases_tmp, sum)
new_covid_deaths <- aggregate(. ~ CountryName, covid_deaths_tmp, sum)

covid_cases  <- gather(new_covid_cases, key = "dateRep",
                       value = "cumcases", 2:480)
covid_deaths <- gather(new_covid_deaths, key = "dateRep",
                       value = "cumdeaths", 2:480)

covid <- merge(covid_cases, covid_deaths, by = c("CountryName", "dateRep"))
rm(covid_cases_tmp, covid_deaths_tmp, new_covid_deaths, new_covid_cases,
   covid_cases, covid_deaths)

covid$dateRep <- as.Date(covid$dateRep, format = "%m/%d/%y")
covid$cases   <- 0
covid$deaths  <- 0
covid$weekday <- weekdays(covid$dateRep)

covid <- covid[order(covid$CountryName, covid$dateRep), ]

covid_list <- list()
for (country in unique(covid$CountryName)){
  cumdeaths_column <- covid[covid$CountryName == country, "cumdeaths"]
  cumcases_column  <- covid[covid$CountryName == country, "cumcases"]
  time_range       <- length(cumdeaths_column)
  covid[covid$CountryName == country, "deaths"] <- c(0, cumdeaths_column[2:time_range] - cumdeaths_column[1:(time_range - 1)])
  covid[covid$CountryName == country, "cases"]  <- c(0, cumcases_column[2:time_range] - cumcases_column[1:(time_range - 1)])
  tmp <- max(covid[covid$CountryName == country, "cumcases"])
  tmp_deaths <- max(cumdeaths_column)
  if (tmp >= 1000 & tmp_deaths >= 100 & (country != "Cambodia") & (country != "Mongolia")){
    #We restrict our attention only to the contries with more than 1000 cases and only starting from 100th case
    tmp_df <- covid[(covid$CountryName == country & covid$cumcases >= 100),
                    c("dateRep", "deaths", "cumdeaths", "cases", "cumcases", "weekday")]
    tmp_index <- match("Monday", tmp_df$weekday)
    if (nrow(tmp_df) > 300) {
      covid_list[[country]] <- tmp_df[tmp_index:nrow(tmp_df), ]
    }
  }
}

#Calculate the number of days that we have data for all countries.
#We are not considering CHN = China as it has too long dataset.
t_len     <- min(sapply(covid_list, NROW))
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
  covid_mat[, i] <- covid_list[[country]]$deaths[1:t_len]
  i = i + 1
}

#Cleaning the data: there are weird cases when the number of new cases is negative! 
sum(covid_mat < 0)
covid_mat[covid_mat < 0] <- 0

m_hat <- function(vect_u, data_p, grid_p, bw){
  m_hat_vec <- c()
  for (u in vect_u){
    result <- sum((((u - grid_p) / bw < 1) & ((u - grid_p) / bw >= -1)) * data_p)
    norm   <- sum((((u - grid_p) / bw < 1) & ((u - grid_p) / bw >= -1)))
    if (norm == 0){
      m_hat_vec <- c(m_hat_vec, 0)
    } else {
      m_hat_vec <- c(m_hat_vec, result/norm)
    }
  }
  return(m_hat_vec)
}

grid_points <- (1:t_len)/sqrt(t_len)

#Step 2
norm   <- c()
a_vec  <- c()
b_vec  <- c()
c_vec  <- c()
norm_p <- c()
for (k in 1:n_ts) {
  norm <- c(norm, integrate(m_hat, lower = - Inf, upper = Inf,
                            data_p = covid_mat[, k], grid_p = grid_points,
                            bw = bw_abs/sqrt(t_len), subdivisions = 2000)$value)
  
  integrand1 <- function(x) {x * m_hat(x, data_p = covid_mat[, k],
                                       grid_p = grid_points,
                                       bw = bw_abs/sqrt(t_len)) / norm[k]}
  a_vec      <- c(a_vec, integrate(integrand1, lower = - Inf, upper = Inf,
                                   subdivisions = 2000)$value)
  
  integrand2 <- function(x) {x * x * m_hat(x, data_p = covid_mat[, k],
                                           grid_p = grid_points,
                                           bw = bw_abs/sqrt(t_len)) / norm[k]}
  tmp        <- integrate(integrand2, lower = - Inf, upper = Inf,
                          subdivisions = 2000)$value
  b_vec      <- c(b_vec, sqrt(tmp - a_vec[k]^2))
  c_vec      <- c(c_vec, norm[k] / b_vec[k])
  integrand3 <- function(x) {m_hat(a_vec[k] + b_vec[k] * x,
                                   data_p = covid_mat[, k], grid_p = grid_points,
                                   bw = bw_abs/sqrt(t_len)) / c_vec[k]}
  norm_p <- c(norm_p, integrate(integrand3, lower = - Inf, upper = Inf,
                                subdivisions = 2000)$value)
}


#Matrix with the distances: Step 3
Delta_hat <- matrix(data = rep(0, n_ts * n_ts), nrow = n_ts, ncol = n_ts)

for (i in 1:(n_ts - 1)){
  p_i_star <- function(x) {(m_hat(a_vec[i] + b_vec[i] * x, data_p = covid_mat[, i],
                                  grid_p = grid_points,
                                  bw = bw_abs/sqrt(t_len)) / c_vec[i]) / norm_p[i]}
  for (j in (i + 1):n_ts){
    p_j_star <- function(x) {(m_hat(a_vec[j] + b_vec[j] * x, data_p = covid_mat[, j],
                                    grid_p = grid_points,
                                    bw = bw_abs/sqrt(t_len)) / c_vec[j]) / norm_p[j]}
    integrand <- function(x) {(sqrt(p_i_star(x)) - sqrt(p_j_star(x)))^2}
    tmp <- integrate(integrand, lower = -Inf, upper = Inf,
                     subdivisions=3000)$value
    Delta_hat[i, j] <- tmp
    Delta_hat[j, i] <- tmp
  }
}

#And now the clustering itself
colnames(Delta_hat) <- countries
rownames(Delta_hat) <- countries

#save(Delta_hat, file = "results_alt_approach_14days_deaths.RData")
load("results_alt_approach_14days_deaths.RData")

delta_dist <- as.dist(Delta_hat)
res        <- hclust(delta_dist)

g_hat <- function(x_vec, covid_mat, inds, a_vec, b_vec, c_vec, norm_p, bw_abs,
                  t_len, grid_points){
  g_hat_vec <- rep(0, length(x_vec))
  for (ind in inds){
    p_i_hat_star <- m_hat_standard(a_vec[ind] + b_vec[ind] * x_vec,
                                   covid_mat[, ind], grid_points,
                                   bw = bw_abs/sqrt(t_len)) / (c_vec[ind] * norm_p[ind])
    g_hat_vec <- g_hat_vec + p_i_hat_star
  }
  return(g_hat_vec/length(inds))
}

max_n_cl <- 20
BIC_vec  <- c()

for (n_cl in 1:max_n_cl){
  subgroups         <- cutree(res, n_cl)
  sigma_hat_vec     <- rep(0, n_ts)
  for (cl in 1:n_cl){
    countries_cluster <- colnames(Delta_hat)[subgroups == cl]
    inds              <- match(countries_cluster, countries)
    for (ind in inds){
      x_vec <- (grid_points - a_vec[ind])/b_vec[ind]
      g_hat_normed <- c_vec[ind] * g_hat(x_vec = x_vec, covid_mat = covid_mat,
                                         inds = inds, a_vec = a_vec, b_vec = b_vec,
                                         c_vec = c_vec, norm_p = norm_p,
                                         bw_abs = bw_abs, t_len = t_len,
                                         grid_point = grid_points)
      sigma_hat_vec[ind] <- mean((covid_mat[, ind] - g_hat_normed)^2)
    }
  }
  BIC <- t_len * sum(log(sigma_hat_vec)) - log(n_ts * t_len) * (n_cl * (n_ts + t_len) + n_ts)
  BIC_vec <- c(BIC_vec, BIC)
}

#results_output_alt(res, covid_mat, Delta_hat, n_cl, countries, path = "plots/deaths/",
#                   bw_abs, grid_points, t_len, a_vec, b_vec, c_vec, norm_p)