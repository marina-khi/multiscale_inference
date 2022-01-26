rm(list=ls())

library(readxl)
library(haven)
library(tidyr)
library(multiscale)
library(dplyr)
library(zoo)

source("functions/functions.R")


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

hp1 <- read_excel("data/JSTdatasetR5.xlsx", sheet = "Data", na = "") #Macroeconomic dataset for 1870-2017
hp2 <- read_dta("data/NPLH.dta")                                     #Knoll et al. dataset

hp  <- merge(select(hp1, 'year', 'iso', 'country', 'cpi', 'rgdppc', 'pop', 'ltrate'),
             select(hp2, 'year', 'iso', 'hpnom'), by = c('iso', "year"), all = TRUE)

hp_data <-
  hp %>%
  subset((year <= 2012) & (year >= 1890)) %>%
  mutate(hpreal = hpnom / (cpi / 100)) %>%   #Calculating real housing prices
  group_by(iso) %>%
  mutate(infl = (cpi - dplyr::lag(cpi, n = 1, default = NA)) / 100) %>%
  mutate(infl = coalesce(infl, 0)) %>%       #Inflation as changes in CPI
  transform(hpreal = na.approx(hpreal)) %>%  #Imputing the missing values
  transform(ltrate = na.approx(ltrate)) %>%
  mutate(log_hp = log(hpreal)) %>%
  mutate(log_gdp = log(rgdppc)) %>%
  mutate(log_pop = log(pop)) %>%
  group_by(iso) %>%
  mutate(delta_log_hp = log(hpreal) - log(dplyr::lag(hpreal, n = 1, default = NA))) %>%
  mutate(delta_log_pop = log(pop) - log(dplyr::lag(pop, n = 1, default = NA))) %>%
  mutate(delta_log_gdp = log(rgdppc) - log(dplyr::lag(rgdppc, n = 1, default = NA))) %>%
  mutate(delta_ltrate = ltrate - dplyr::lag(ltrate, n = 1, default = NA)) %>%
  mutate(delta_infl = infl - dplyr::lag(infl, n = 1, default = NA)) %>%
  subset(iso %in% c('AUS', 'BEL', 'DNK', 'FRA', 'NLD', 'NOR', 'SWE', 'USA'))


####################################
#Necessary data-dependent variables#
####################################

ticks     <- c(1, 31, 61, 91, 121)
dates     <- unique(hp_data$year)
n_ts      <- length(unique(hp_data$iso))
t_len     <- nrow(hp_data) / n_ts
countries <- unique(hp_data$iso)

# colSums(hp_data[, 4:ncol(hp_data)])
# 
# for (country in countries){
#   tmp <- hp_data[hp_data$iso == country, ]
#   cat("Country: ", country, "; ", sum(is.na(tmp$hpreal)), "\n")
# }


##########################
#Plotting the time series#
##########################

pdf("plots/real_housing_prices.pdf", width=10, height=10, paper="special")
par(mfrow = c(3,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
par(mar = c(0, 0.5, 1.5, 0)) #Margins for each plot
par(oma = c(2.5, 1.5, 0.2, 0.2)) #Outer margins

plot(NA, ylab="", xlab = "", xlim = c(0, t_len),
     ylim = c(0, max(hp_data$hpreal, na.rm = TRUE)), xaxt = 'n',
     mgp=c(2,0.5,0), cex = 1.2, tck = -0.025, main = "Real house prices")
for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  tmp <- tmp[order(tmp$year),]
  lines(1:t_len, tmp$hpreal)
}

plot(NA, ylab="", xlab = "", xlim = c(0,t_len),
     ylim = c(1, max(hp_data$log_hp, na.rm = TRUE)), xaxt = 'n',
     mgp=c(2,0.5,0), cex = 1.2, tck = -0.025,
     main = "Logarithm of the real house prices")
for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  tmp <- tmp[order(tmp$year),]
  lines(1:t_len, tmp$log_hp)
}

plot(NA, ylab="", xlab = "", xlim = c(0,t_len), ylim = c(-0.5, 0.75),
     xaxt = 'n', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025,
     main = "Growth rate of the logarithm of the real house prices")
for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  tmp <- tmp[order(tmp$year),]
  lines(1:t_len, tmp$delta_log_hp)
}
axis(1, at = ticks, labels = dates[ticks])

dev.off()


#################################
#Parameters and necessary values#
################################$

hp_log           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(hp_log) <- countries

hp_log_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(hp_log_augm) <- countries

hp_growth_rate           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(hp_growth_rate) <- countries

hp_growth_rate_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(hp_growth_rate_augm) <- countries

beta_log  <- matrix(data = NA, nrow = 4, ncol = n_ts)
alpha_log <- c()

beta_growth_rate  <- matrix(data = NA, nrow = 4, ncol = n_ts)
alpha_growth_rate <- c()

i         <- 1

for (country in countries){
  tmp <- hp_data[hp_data$iso == country, ]
  tmp <- tmp[order(tmp$year),]
  
  #Calculating first difference of the growth rate
  tmp <-
    tmp %>%
    mutate(delta_log_hp = coalesce(delta_log_hp, 0)) %>%
    mutate(delta_log_gdp = coalesce(delta_log_gdp, 0)) %>%
    mutate(delta_ltrate = coalesce(delta_ltrate, 0)) %>%
    mutate(delta_log_pop = coalesce(delta_log_pop, 0)) %>%
    mutate(delta_infl = coalesce(delta_infl, 0)) %>%
    mutate(delta_delta_hp = delta_log_hp - dplyr::lag(delta_log_hp, n = 1, default = NA))%>%
    mutate(delta_delta_gdp = delta_log_gdp - dplyr::lag(delta_log_gdp, n = 1, default = NA))%>%
    mutate(delta_delta_pop = delta_log_pop - dplyr::lag(delta_log_pop, n = 1, default = NA))%>%
    mutate(delta_delta_rate = delta_ltrate - dplyr::lag(delta_ltrate, n = 1, default = NA)) %>%
    mutate(delta_delta_infl = delta_infl - dplyr::lag(delta_infl, n = 1, default = NA))

  #Estimating beta_i
  y_vec_tmp <- as.matrix(tmp[-1, 'delta_log_hp'])
  x_mat_tmp <- as.matrix(tmp[-1, c('delta_log_gdp', 'delta_log_pop',
                                   'delta_ltrate', 'delta_infl')])

  y_vec_growth_rate_tmp <- as.matrix(tmp[-1, 'delta_delta_hp'])
  x_mat_growth_rate_tmp <- as.matrix(tmp[-1 , c('delta_delta_gdp', 'delta_delta_pop',
                                                'delta_delta_rate', 'delta_delta_infl')])
  
  
  beta_tmp  <- solve(t(x_mat_tmp) %*% x_mat_tmp) %*% t(x_mat_tmp) %*% y_vec_tmp
  beta_log[, i] <- beta_tmp
  
  beta_growth_rate_tmp  <- solve(t(x_mat_growth_rate_tmp) %*% x_mat_growth_rate_tmp) %*% t(x_mat_growth_rate_tmp) %*% y_vec_growth_rate_tmp
  beta_growth_rate[, i] <- beta_growth_rate_tmp
  
  #Estimating alpha_i
  alpha_growth_rate_tmp <- mean(tmp$delta_log_hp - as.vector(as.matrix(tmp[, c('delta_log_gdp',
                                                                               'delta_log_pop',
                                                                               'delta_ltrate',
                                                                               'delta_infl')]) %*% beta_growth_rate_tmp))
  alpha_growth_rate[i]  <- alpha_growth_rate_tmp
  
  alpha_tmp     <- mean(tmp$log_hp - as.vector(as.matrix(tmp[, c('log_gdp',
                                                                 'log_pop',
                                                                 'ltrate',
                                                                 'infl')]) %*% beta_tmp))
  alpha_log[i]  <- alpha_tmp
  
  
  #Working with adjusted time series and storing the original one
  y_vec_growth_rate_adj    <- tmp$delta_log_hp - as.vector(as.matrix(tmp[, c('delta_log_gdp',
                                                                             'delta_log_pop',
                                                                             'delta_ltrate',
                                                                             'delta_infl')]) %*% beta_growth_rate_tmp) - alpha_growth_rate_tmp
  hp_growth_rate[, i]      <- tmp$delta_log_hp
  hp_growth_rate_augm[, i] <- y_vec_growth_rate_adj
  
  y_vec_adj        <- tmp$log_hp - as.vector(as.matrix(tmp[, c('log_gdp',
                                                               'log_pop',
                                                               'ltrate',
                                                               'infl')]) %*% beta_tmp) - alpha_tmp
  hp_log[, i]      <- tmp$log_hp
  hp_log_augm[, i] <- y_vec_adj
  i = i + 1
}

################################################################################
########################  LOG of house prices  #################################
################################################################################

#####################
#Estimating variance#
#####################

#Order selection
q_vec <- 10:20
r_vec <- 10:15
order_results <- c()

for (j in 1:n_ts){
  criterion_matrix <- expand.grid(q = q_vec, r = r_vec)
  
  criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
  criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
  
  for (i in 1:nrow(criterion_matrix)){
    FPE <- c()
    AIC <- c()
    AICC <- c()
    SIC <- c()
    HQ <- c()
    
    different_orders <- (1:9)
    
    for (order in different_orders){
      AR.struc      <- estimate_lrv(data = hp_log_augm[, j], q = criterion_matrix$q[[i]],
                                    r_bar = criterion_matrix$r[[i]], p = order)
      sigma_eta_hat <- sqrt(AR.struc$vareta)
      FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
      AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
      AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order +2)/t_len))
      SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
      HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
    }
    criterion_matrix$FPE[[i]]  <- which.min(FPE)
    criterion_matrix$AIC[[i]]  <- which.min(AIC)
    criterion_matrix$AICC[[i]] <- which.min(AICC)
    criterion_matrix$SIC[[i]]  <- which.min(SIC)
    criterion_matrix$HQ[[i]]   <- which.min(HQ)
  }
  maxim <- max(criterion_matrix[, 3:7])
  order_results <- c(order_results, maxim)
  cat("For the country ", colnames(hp_log_augm)[j],
      "for the log of house prices the results are as follows: ", max(criterion_matrix$FPE), " ",
      max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ",
      max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
}

#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 1:n_ts){
  AR.struc        <- estimate_lrv(data = hp_log_augm[, i], q = q, r_bar = r,
                                  #p = order_results[i])  
                                  p = 1)
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}


#########
#Testing#
#########

#Constructing the grid
u_grid      <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
h_grid      <- seq(from = 5 / t_len, to = 1 / 4, by = 5 / t_len)
h_grid      <- h_grid[h_grid > log(t_len) / t_len]
grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)

result <- multiscale_test(data = hp_log_augm, sigma_vec = sigmahat_vector,
                          alpha = alpha,  n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)


#####################
#Plots for the paper#
#####################

#Producing the smoothed curves using local linear estimator

produce_smoothed_plots(matrix = hp_log,
                       pdfname = "plots/smoothed_hp_data.pdf",
                       y_min = min(hp_log), y_max = max(hp_log),
                       ticks_at =  ticks, ticks_labels = dates[ticks])

produce_smoothed_plots(matrix = hp_log_augm,
                       pdfname = "plots/smoothed_hp_data_augmented.pdf",
                       y_min = min(hp_log_augm), y_max = max(hp_log_augm),
                       ticks_at =  ticks, ticks_labels = dates[ticks])


#Producing plots with the final results
for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  if (result$stat_pairwise[i, j] > result$quant){
    produce_plots_hp(results = result, data_i = hp_log_augm[, i],
                     data_j = hp_log_augm[, j], at_ = ticks,
                     labels_ = dates[ticks], name_i = countries[i],
                     name_j = countries[j])
  }
}

# ################################################################################
# #################  GROWTH RATE of the log of house prices  #####################
# ################################################################################
# 
# #####################
# #Estimating variance#
# #####################
# 
# #Order selection
# q_vec <- 10:20
# r_vec <- 10:15
# order_results_2 <- c()
# 
# for (j in 1:n_ts){
#   criterion_matrix <- expand.grid(q = q_vec, r = r_vec)
#   
#   criterion_matrix$FPE  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$AIC  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$SIC  <- numeric(length = nrow(criterion_matrix))
#   criterion_matrix$HQ   <- numeric(length = nrow(criterion_matrix))
#   
#   for (i in 1:nrow(criterion_matrix)){
#     FPE <- c()
#     AIC <- c()
#     AICC <- c()
#     SIC <- c()
#     HQ <- c()
#     
#     different_orders <- (1:9)
#     
#     for (order in different_orders){
#       AR.struc      <- estimate_lrv(data = hp_growth_rate_augm[, j],
#                                     q = criterion_matrix$q[[i]],
#                                     r_bar = criterion_matrix$r[[i]], p = order)
#       sigma_eta_hat <- sqrt(AR.struc$vareta)
#       FPE <- c(FPE, (sigma_eta_hat^2 * (t_len + order)) / (t_len - order))
#       AIC <- c(AIC, t_len * log(sigma_eta_hat^2) + 2 * order)
#       AICC <- c(AICC, t_len * log(sigma_eta_hat^2) + t_len * (1 + order / t_len)/(1 - (order + 2)/t_len))
#       SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(t_len) / t_len)
#       HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(t_len)) / t_len)
#     }
#     criterion_matrix$FPE[[i]]  <- which.min(FPE)
#     criterion_matrix$AIC[[i]]  <- which.min(AIC)
#     criterion_matrix$AICC[[i]] <- which.min(AICC)
#     criterion_matrix$SIC[[i]]  <- which.min(SIC)
#     criterion_matrix$HQ[[i]]   <- which.min(HQ)
#   }
#   maxim <- max(criterion_matrix[, 3:7])
#   order_results_2 <- c(order_results_2, maxim)
#   cat("For the country ", colnames(hp_growth_rate_augm)[j],
#       "for the growth rate of house prices the results are as follows: ",
#       max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ",
#       max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ",
#       max(criterion_matrix$HQ), " \n")
# }
# 
# #Calculating each sigma_i separately
# sigmahat_vector_2 <- c()
# for (i in 1:n_ts){
#   AR.struc        <- estimate_lrv(data = hp_log_augm[, i], q = q, r_bar = r,
#                                   #p = order_results_2[i])  
#                                   p = 1)
#   sigma_hat_i     <- sqrt(AR.struc$lrv)
#   sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_hat_i)
# }
# 
# 
# #########
# #Testing#
# #########
# 
# #Constructing the grid
# u_grid      <- seq(from = 1 / t_len, to = 1, by = 1 / t_len)
# h_grid      <- seq(from = 5 / t_len, to = 1 / 4, by = 1 / t_len)
# h_grid      <- h_grid[h_grid > log(t_len) / t_len]
# grid        <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
# 
# result <- multiscale_test(data = hp_growth_rate_augm,
#                           sigma_vec = sigmahat_vector_2, alpha = alpha,
#                           n_ts = n_ts, grid = grid, sim_runs = sim_runs,
#                           epidem = FALSE)
# 
# 
# #####################
# #Plots for the paper#
# #####################
# 
# #Producing the smoothed curves using local linear estimator
# 
# 
# produce_smoothed_plots(matrix = hp_growth_rate,
#                        pdfname = "plots/growth_rate/smoothed_hp_data.pdf",
#                        y_min = min(hp_growth_rate), y_max = max(hp_growth_rate),
#                        ticks_at =  ticks, ticks_labels = dates[ticks])
# 
# produce_smoothed_plots(matrix = hp_growth_rate_augm,
#                        pdfname = "plots/growth_rate/smoothed_hp_data_augmented.pdf",
#                        y_min = min(hp_growth_rate_augm),
#                        y_max = max(hp_growth_rate_augm),
#                        ticks_at =  ticks, ticks_labels = dates[ticks])
# 
# 
# #Producing plots with the final results
# for (l in seq_len(nrow(result$ijset))){
#   i <- result$ijset[l, 1]
#   j <- result$ijset[l, 2]
#   if (result$stat_pairwise[i, j] > result$quant){
#     produce_plots(results = result, data_i = hp_growth_rate_augm[, i],
#                   data_j = hp_growth_rate_augm[, j], at_ = ticks,
#                   labels_ = dates[ticks], name_i = countries[i],
#                   name_j = countries[j])
#   }
# }
# 
