rm(list=ls())

library(readxl)
library(haven)
library(tidyr)
library(multiscale)
library(dplyr)


######################
#Necessary parameters#
######################
alpha     <- 0.05
sim_runs  <- 5000
q         <- 15 #Parameters for the estimation of sigma
r         <- 10


##############
#Data loading#
##############

hp  <- read_excel("data/JSTdatasetR5.xlsx", sheet = "Data", na = "") #Macroeconomic dataset for 1870-2017
hp2 <- read_dta("data/NPLH.dta")                                     #Knoll et al. dataset

hp <- 
  hp %>%
  select('year', 'iso', 'country', 'cpi', 'rgdppc', 'pop', 'ltrate')

hp2 <-
  hp2 %>%
  select('year', 'iso', 'hpnom')

hp <- merge(hp, hp2, by = c('iso', "year"), all = TRUE)

hp_data <-
  hp %>%
  mutate(hpreal = hpnom / (cpi / 100)) %>%   #Calculating real housing prices
  mutate(infl = cpi - dplyr::lag(cpi, n = 1, default = NA)) %>%
  mutate(infl = coalesce(infl, 0)) %>%       #Inflation as changes in CPI
  mutate(delta_log_pop = log(pop) - log(dplyr::lag(pop, n = 1, default = NA))) %>%
  mutate(delta_log_gdp = log(rgdppc) - log(dplyr::lag(rgdppc, n = 1, default = NA))) %>%
  mutate(delta_ltrate = ltrate - dplyr::lag(ltrate, n = 1, default = NA)) %>%
  mutate(delta_infl = infl - dplyr::lag(infl, n = 1, default = NA)) %>%
  mutate(delta_log_pop = coalesce(delta_log_pop, 0)) %>%
  mutate(delta_log_gdp = coalesce(delta_log_gdp, 0)) %>%
  mutate(delta_ltrate = coalesce(delta_ltrate, 0)) %>%
  mutate(delta_infl = coalesce(delta_infl, 0)) %>%
  subset((year <= 2012) & (year >= 1890)) %>%
  subset(iso %in% c('AUS', 'BEL', 'DNK', 'FRA', 'NLD', 'NOR', 'SWE', 'USA')) %>%
  group_by(iso) %>%
  transform(hpreal = na.approx(hpreal)) %>%
  mutate(log_hp = log(hpreal)) %>%
  mutate(log_gdp = log(rgdppc)) %>%
  mutate(log_pop = log(pop)) %>%
  mutate(delta_log_hp = log(hpreal) - log(dplyr::lag(hpreal, n = 1, default = NA))) %>%
  mutate(delta_log_hp = coalesce(delta_log_hp, 0))


#Variables
dates     <- unique(hp_data$year)
n_ts      <- length(unique(hp_data$iso))
t_len     <- nrow(hp_data) / n_ts
countries <- unique(hp_data$iso)

for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  cat("Country: ", country, "; ", sum(is.na(tmp$hpreal)), "\n")
}


##########################
#Plotting the time series#
##########################

pdf("plots/real_housing_prices.pdf", width=10, height=10, paper="special")
par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
par(mar = c(0, 0.5, 1.5, 0)) #Margins for each plot
par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins

plot(NA, ylab="", xlab = "", xlim = c(0, t_len), ylim = c(0, max(hp_data$hpreal)),
     xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2,
     tck = -0.025, main = "Real house prices")
for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  tmp <- tmp[order(tmp$year),]
  lines(1:t_len, tmp$hpreal)
}

plot(NA, ylab="", xlab = "", xlim = c(0,t_len), ylim = c(1, max(hp_data$log_hp)),
     xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2,
     tck = -0.025, main = "Logarithm of the real house prices")
for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  tmp <- tmp[order(tmp$year),]
  lines(1:t_len, tmp$log_hp)
}

plot(NA, ylab="", xlab = "", xlim = c(0,t_len), ylim = c(-0.5, 0.75),
     xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2,
     tck = -0.025, main = "Growth rate of the logarithm of the real house prices")
for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  tmp <- tmp[order(tmp$year),]
  lines(1:t_len, tmp$delta_log_hp)
}
dev.off()


#################################
#Parameters and necessary values#
################################$

hp_mat_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(hp_mat_augm) <- countries

# hp_mat_growth           <- matrix(NA, ncol = n_ts, nrow = t_len)
# colnames(hp_mat_growth) <- countries

hp_mat_log           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(hp_mat_log) <- countries

beta      <- matrix(data = NA, nrow = 4, ncol = n_ts)
alpha_vec <- c()
i         <- 1
for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  tmp <- tmp[order(tmp$year),]
  
  #Calculating first differences
  tmp <- 
    tmp %>%
    mutate(delta_delta_hp = delta_log_hp - dplyr::lag(delta_log_hp, n = 1, default = NA))%>%
    mutate(delta_delta_gdp = delta_log_gdp - dplyr::lag(delta_log_gdp, n = 1, default = NA))%>%
    mutate(delta_delta_pop = delta_log_pop - dplyr::lag(delta_log_pop, n = 1, default = NA))%>%
    mutate(delta_delta_rate = delta_ltrate - dplyr::lag(delta_ltrate, n = 1, default = NA)) %>%
    mutate(delta_delta_infl = delta_infl - dplyr::lag(delta_infl, n = 1, default = NA))
  
  #Estimating beta_i
  # y_vec_tmp <- as.matrix(tmp[-1, 'delta_delta_hp'])
  # x_mat_tmp <- as.matrix(tmp[-1 , c('delta_delta_gdp', 'delta_delta_pop',
  #                                   'delta_delta_rate', 'delta_delta_infl')])
  y_vec_tmp <- as.matrix(tmp[, 'delta_log_hp'])
  x_mat_tmp <- as.matrix(tmp[ , c('delta_log_gdp', 'delta_log_pop',
                                    'delta_ltrate', 'delta_infl')])
  
  beta_tmp  <- solve(t(x_mat_tmp) %*% x_mat_tmp) %*% t(x_mat_tmp) %*% y_vec_tmp
  beta[, i] <- beta_tmp
  
  #Estimating alpha_i
  # alpha_tmp    <- mean(tmp$delta_log_hp - as.vector(as.matrix(tmp[, c('delta_log_gdp', 'delta_log_pop',
  #                                                                     'delta_ltrate', 'delta_infl')]) %*% beta_tmp))
  alpha_tmp    <- mean(tmp$log_hp - as.vector(as.matrix(tmp[, c('log_gdp', 'log_pop',
                                                                'ltrate', 'infl')]) %*% beta_tmp))
  
  alpha_vec[i] <- alpha_tmp
  
  #Working with adjusted time series and storing the original one
  # y_vec_adj             <- tmp$delta_log_hp - as.vector(as.matrix(tmp[, c('delta_log_gdp', 'delta_log_pop',
  #                                                                         'delta_ltrate', 'delta_infl')]) %*% beta_tmp) - alpha_tmp
  # hp_mat_growth[, i]   <- tmp$delta_log_hp
  y_vec_adj         <- tmp$log_hp - as.vector(as.matrix(tmp[, c('log_gdp', 'log_pop',
                                                                'ltrate', 'infl')]) %*% beta_tmp) - alpha_tmp
  hp_mat_log[, i]   <- tmp$log_hp
  hp_mat_augm[, i]  <- y_vec_adj
  i = i + 1
}

#Estimating the variance
source("sigma_estimation.R")
sigma_vec <- sigma(gdp_mat_augm, n_ts = n_ts, q = q, r = r)


#########
#Testing#
#########

#Constructing the grid
u_grid      <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
h_grid      <- seq(from = 5 / t_len, to = 1 / 4, by = 5 / t_len)
h_grid      <- h_grid[h_grid > log(t_len) / t_len]
grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)


result <- multiscale_test(data = gdp_mat_augm,
                          sigma_vec = sigmahat_vector,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

#####################
#Plots for the paper#
#####################

#Producing the smoothed curves using local linear estimator

source("functions.R")
pdfname <- "plots/smoothed_gdp_data.pdf"
produce_smoothed_plots(gdp_mat_growth, pdfname, dates)

if (capital_variable == "gfcf") {
  pdfname_augm <- "plots/smoothed_gdp_data_augmented_gfcf.pdf"
} else {
  pdfname_augm <- "plots/smoothed_gdp_data_augmented.pdf"
}
produce_smoothed_plots(gdp_mat_augm, pdfname_augm, dates)

#Producing plots with the final results

grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
ticks <- c(0, 40, 80, 120)
for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  produce_plots(results = result, l = l, data_i, data_j, dates = dates,
                name_i = countries[i], name_j = countries[j], filename)
}

