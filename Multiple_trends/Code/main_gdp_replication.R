rm(list=ls())

library(tidyr)
library(multiscale)
library(dplyr)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")


######################
#Necessary parameters#
######################
alpha     <- 0.05
sim_runs  <- 5000
q         <- 15 #Parameters for the estimation of sigma
r         <- 10
countries <- c("AUT", "AUS", "CAN", "CHE", "DEU", "FIN", "FRA", "GBR",
               "JPN", "NOR", "USA") #Necessary for loading the data

#capital_variable <- "stock"
capital_variable <- "gfcf" #Either we are using capital stock ("stock") or
                           #gross fixed capital formation ("gfcf")
                           #as the variable k


##############
#Data loading#
##############

source("functions/data_loading.R") #Now the matrix X_mat_filled contains all the data

#Variables
dates     <- unique(X_mat_filled$date)
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
source("functions/parameters_estimation.R")
estimated      <- parameters(X_mat_filled, n_ts, t_len, countries)
gdp_mat_growth <- estimated$gdp_mat_growth
gdp_mat_augm   <- estimated$gdp_mat_augm

#Estimating the variance
source("functions/sigma_estimation.R")
sigma_vec <- sigma(gdp_mat_augm, n_ts = n_ts, q = q, r = r)


#########
#Testing#
#########

#Constructing the grid
u_grid      <- seq(from = 5 / (2 * t_len), to = 1, by = 2 / t_len)
h_grid      <- seq(from = 4 / t_len, to = 1 / 4, by = 4 / t_len)
h_grid      <- h_grid[h_grid > log(t_len) / t_len]
grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

result <- multiscale_test(data = gdp_mat_augm,
                          sigma_vec = sigma_vec,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

#####################
#Plots for the paper#
#####################

#Producing the smoothed curves using local linear estimator

source("functions/functions.R")
at <- seq(5, 125, by = 20)

pdfname <- "plots/smoothed_gdp_data.pdf"
produce_smoothed_plots(gdp_mat_growth, pdfname,
                       y_min = min(gdp_mat_growth) + 0.035,
                       y_max = max(gdp_mat_growth) - 0.02,
                       ticks_at =  at, ticks_labels = dates[at],
                       yaxp_ = c(-0.02, 0.02, 4))

pdfname_augm <- "plots/smoothed_gdp_data_augmented.pdf"
produce_smoothed_plots(gdp_mat_augm, pdfname_augm,
                       y_min = min(gdp_mat_augm) + 0.03,
                       y_max = max(gdp_mat_augm) - 0.02,
                       ticks_at =  at, ticks_labels = dates[at],
                       yaxp_ = c(-0.03, 0.01, 4))

#Producing plots with the final results
for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]

  if (result$stat_pairwise[i, j] > result$quant) {
    if ((countries[i] %in% c("NOR", "USA")) & (countries[j] %in% c("USA", "NOR"))){
      #For color consistency in the paper
      i <- result$ijset[l, 2]
      j <- result$ijset[l, 1]
    }
    produce_plots(results = result, data_i = gdp_mat_augm[, i],
                  data_j = gdp_mat_augm[, j], ticks_ = at,
                  labels_ = dates[at], name_i = countries[i],
                  name_j = countries[j])
  }
}