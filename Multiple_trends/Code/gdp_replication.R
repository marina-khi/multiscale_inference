rm(list=ls())

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
countries <- c("AUT", "AUS", "CAN", "CHE", "DEU", "FIN", "FRA", "GBR",
               "JPN", "NOR", "USA") #Necessary for loading the data

capital_variable <- "stock"
#capital_variable <- "gfcf" #Either we are using capital stock ("stock") or
                           #gross fixed capital formation ("gfcf")
                           #as the variable k


##############
#Data loading#
##############

source("data_loading.R") #Now we have the matrix X_mat_filled with all the data

#Variables
dates     <- unique(X_mat_filled$date)
dates     <- dates[dates >= as.Date("01-01-1976", format = "%d-%m-%Y")]
n_ts      <- length(unique(X_mat_filled$WBcode))
t_len     <- nrow(X_mat_filled) / n_ts

if (capital_variable == "stock") {
  X_mat_filled$k_it <- X_mat_filled$stock_it
} else {
  X_mat_filled$k_it <- X_mat_filled$gfcf_it
}


#################################
#Parameters and necessary values#
################################$

#Estimating alpha and beta
source("parameters_estimation.R")
estimated      <- parameters(X_mat_filled, n_ts, t_len, countries)
gdp_mat_growth <- estimated$gdp_mat_growth
gdp_mat_augm   <- estimated$gdp_mat_augm

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

