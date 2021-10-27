rm(list=ls())

library(tidyr)
library(multiscale)
library(zoo)
library(dplyr)

#Defining necessary constants
alpha    <- 0.05 #confidence level for application
sim_runs <- 5000 #Number of simulation runs to produce the Gaussian qauntiles

################################
#Loading the human capital data#
################################

h      <- read.csv("data/human_capital_long.csv", sep = ",", dec = ".",
                   stringsAsFactors = FALSE, na.strings = "")
h$date <- as.Date(paste0("01-01-", h$year), format = "%d-%m-%Y")
h      <- subset(h, select = - c(BLcode, sex, agefrom, ageto, year, region_code))

h_data <- complete(h, date = seq.Date(min(date), max(date), by='quarter'),
                   country)
h_data <- h_data[order(h_data$country, h_data$date),]

h_data        <- fill(h_data, WBcode, .direction = "up")
h_data$yr_sch <- na.approx(h_data$yr_sch)

h_data <- 
  h_data %>%
  group_by(WBcode) %>%
  mutate(delta_h_it = log(yr_sch) - log(dplyr::lag(yr_sch, n = 1, default = NA)))

h_data <- subset(h_data, select = c('date', 'WBcode', 'delta_h_it'))

rm(h)

#########################
#Loading the labour data#
#########################

l           <- read.csv("data/labour.csv", sep = ",", dec = ".",
                           stringsAsFactors = FALSE, na.strings = "")
l_data      <- subset(l, Subject == 'Employed population, Aged 15 and over, All persons')
l_data$date <- as.Date(as.yearqtr(l_data$TIME, format = '%Y-Q%q'), format = "%d-%m-%Y")
l_data <- complete(l_data, date = seq.Date(min(date), max(date), by='quarter'),
                   LOCATION)
l_data      <- l_data[order(l_data$LOCATION, l_data$date),]
l_data$Value <- na.approx(l_data$Value)


l_data <- 
  l_data %>%
  group_by(LOCATION) %>%
  mutate(delta_l_it = log(Value) - log(dplyr::lag(Value, n = 1, default = NA)))

l_data <- subset(l_data, select = c('date', 'LOCATION', 'delta_l_it'))
colnames(l_data) <- c('date', 'WBcode', 'delta_l_it')

X_mat <- merge(h_data, l_data, all = TRUE)
rm(l)

######################
#Loading the gdp data#
######################

gdp      <- read.delim("data/GDP_Quarterly_Quarterly.txt",
                       stringsAsFactors = FALSE, na.strings = "")
colnames(gdp) <- c("date", "AUT", "AUS", "BEL", "BRA", "CAN", "CHE", "CHL", "CZE", "DEU",
                   "DEU_log", "DNK", "EST", "ESP", "FIN", "FRA", "GBR", "GRC",
                   "HUN","IND", "IRL", "ISR", "IND", "ITA", "JPN", "KOR", "LUX",
                   "MEX", "NLD", "NOR", "NZL", "POL", "PRT", "RUS", "SWE", "SVN",
                   "SVK", "TUR", "USA", "ZAF")
gdp      <- subset(gdp, select = -c(DEU_log))
gdp$date <- as.Date(gdp$date, format = "%Y-%m-%d")

gdp_data <- gather(gdp, WBcode, Value, AUT:ZAF, factor_key = FALSE)
gdp_data <- gdp_data[order(gdp_data$WBcode, gdp_data$date),]

gdp_data <- 
  gdp_data %>%
  group_by(WBcode) %>%
  mutate(delta_gdp_it = log(Value) - log(dplyr::lag(Value, n = 1, default = NA)))

gdp_data <- subset(gdp_data, select = c('date', 'WBcode', 'delta_gdp_it'))

X_mat <- merge(X_mat, gdp_data, all = TRUE)
rm(gdp)

################################
#Loading the capital stock data#
################################

k      <- read.delim("data/Capital_stock_Annual.txt",
                     stringsAsFactors = FALSE, na.strings = "")
colnames(k) <- c("date", "AUT", "AUS", "BEL", "CAN", "CHE", "DEU",
                 "DNK", "ESP", "FIN", "FRA", "GBR", "GRC",
                 "IRL", "ISL", "ITA", "JPN", "LUX",
                 "NLD", "NOR", "NZL", "PRT", "SWE", "TUR", "USA")
k$date <- as.Date(k$date, format = "%Y-%m-%d")

k_data       <- gather(k, WBcode, Value, AUT:USA, factor_key = FALSE)
k_data       <- complete(k_data, date = seq.Date(min(date), max(date), by='quarter'),
                         WBcode)
k_data       <- k_data[order(k_data$WBcode, k_data$date),]
k_data$Value <- na.approx(k_data$Value)

k_data <- 
  k_data %>%
  group_by(WBcode) %>%
  mutate(delta_k_it = log(Value) - log(dplyr::lag(Value, n = 1, default = NA)))

k_data <- subset(k_data, select = c('date', 'WBcode', 'delta_k_it'))

X_mat  <- merge(X_mat, k_data, all = TRUE)
rm(k)

##########################
#All of the data together#
##########################


X_mat <- subset(X_mat, date > as.Date("30-03-1980", format = "%d-%m-%Y") & date < as.Date("02-01-2015", format = "%d-%m-%Y"))
#X_mat <- subset(X_mat, WBcode %in% c("AUT", "AUS", "BEL", "CAN", "CHE", "DEU",
#                                   "DNK", "ESP", "FIN", "FRA", "GBR", "GRC",
#                                   "IRL", "ITA", "JPN", "LUX",
#                                   "NLD", "NOR", "NZL", "PRT", "SWE", "TUR", "USA"))

#X_mat <- subset(X_mat, WBcode %in% c("AUT", "AUS", "CAN", "DEU",
#                                     "DNK", "ESP", "FRA", "GBR", "GRC",
#                                     "IRL", "ITA", "JPN",
#                                     "NLD", "NZL", "PRT", "USA"))

X_mat     <- subset(X_mat,WBcode %in% c("AUT", "AUS", "CAN", "DEU",
                                        "GBR", "JPN", "USA"))
countries <- unique(X_mat$WBcode)
country_names <- c("Australia", "Austria", "Canada", "Germany", "UK", "Japan", "USA")
dates     <- unique(X_mat$date)
n_ts      <- length(unique(X_mat$WBcode))
t_len     <- nrow(X_mat) / n_ts


#############################
#Estimating the coefficients#
#############################

gdp_mat           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(gdp_mat) <- countries

gdp_mat_original           <- matrix(NA, ncol = n_ts, nrow = t_len)
colnames(gdp_mat_original) <- countries


beta      <- matrix(data = NA, nrow = 3, ncol = n_ts)
alpha_vec <- c()
i         <- 1

for (country in countries){
  tmp <- X_mat[X_mat$WBcode == country, ]
  tmp <- tmp[order(tmp$date),]
  
  #Calculating first differences
  tmp <- 
    tmp %>%
    mutate(delta_delta_gdp = delta_gdp_it- dplyr::lag(delta_gdp_it, n = 1, default = NA))%>%
    mutate(delta_delta_h = delta_h_it- dplyr::lag(delta_h_it, n = 1, default = NA))%>%
    mutate(delta_delta_l = delta_l_it- dplyr::lag(delta_l_it, n = 1, default = NA))%>%
    mutate(delta_delta_k = delta_k_it- dplyr::lag(delta_k_it, n = 1, default = NA))

  #Estimating beta_i
  y_vec_tmp <- as.matrix(tmp[-1, 'delta_delta_gdp'])
  x_mat_tmp <- as.matrix(tmp[-1 , c('delta_delta_h', 'delta_delta_l', 'delta_delta_k')])
  beta_tmp  <- solve(t(x_mat_tmp) %*% x_mat_tmp) %*% t(x_mat_tmp) %*% y_vec_tmp
  beta[, i] <- beta_tmp
  
  #Estimating alpha_i
  alpha_tmp    <- mean(tmp$delta_gdp_it - as.vector(as.matrix(tmp[, c('delta_h_it', 'delta_l_it', 'delta_k_it')]) %*% beta_tmp))
  alpha_vec[i] <- alpha_tmp

  #Working with adjusted time series and storing the original one
  y_vec_adj             <- tmp$delta_gdp_it - as.vector(as.matrix(tmp[, c('delta_h_it', 'delta_l_it', 'delta_k_it')]) %*% beta_tmp) - alpha_tmp
  gdp_mat_original[, i] <- tmp$delta_gdp_it
  gdp_mat[, i]          <- y_vec_adj
  i = i + 1
}

rm(tmp, y_vec_tmp, x_mat_tmp, beta_tmp, alpha_tmp)

############################
#Plots for the presentation#
############################

#Original time series

country1 <- which(colnames(gdp_mat_original) == "AUT")
country2 <- which(colnames(gdp_mat_original) == "DEU")

pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], ".pdf"),
    width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = dates, y = gdp_mat_original[, country1],
     ylim=c(min(gdp_mat_original[, country1], gdp_mat_original[, country2]),
            max(gdp_mat_original[, country1], gdp_mat_original[, country2] + 0.03)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
lines(x = dates, y = gdp_mat_original[, country2], col="#604c38")
legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()


pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], "_1.pdf"),
    width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = dates, y = gdp_mat_original[, country1],
     ylim=c(min(gdp_mat_original[, country1], gdp_mat_original[, country2]),
            max(gdp_mat_original[, country1], gdp_mat_original[, country2] + 0.03)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
rect(xleft=as.Date("1988-04-01", format = "%Y-%m-%d"),
     xright = as.Date("1994-10-01", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
lines(x = dates, y = gdp_mat_original[, country1], col="#EB811B")
lines(x = dates, y = gdp_mat_original[, country2], col="#604c38")
segments(as.Date("1988-04-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("1994-10-01", format = "%Y-%m-%d"), par("usr")[3],
         col="#604c38", lwd = 3)
legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()


pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], "_2.pdf"),
    width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = dates, y = gdp_mat_original[, country1],
     ylim=c(min(gdp_mat_original[, country1], gdp_mat_original[, country2]),
            max(gdp_mat_original[, country1], gdp_mat_original[, country2] + 0.03)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
rect(xleft=as.Date("1988-04-01", format = "%Y-%m-%d"),
     xright = as.Date("1994-10-01", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
rect(xleft=as.Date("2006-01-01", format = "%Y-%m-%d"),
     xright = as.Date("2010-04-01", format = "%Y-%m-%d"),
     ybottom=par("usr")[3], ytop=par("usr")[4],
     density=40, col = "#d3d3d3", border = NA)
lines(x = dates, y = gdp_mat_original[, country1], col="#EB811B")
lines(x = dates, y = gdp_mat_original[, country2], col="#604c38")
segments(as.Date("1988-04-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("1994-10-01", format = "%Y-%m-%d"), par("usr")[3],
         col="#604c38", lwd = 3)
segments(as.Date("2006-01-01", format = "%Y-%m-%d"), par("usr")[3],
         as.Date("2010-04-01", format = "%Y-%m-%d"), par("usr")[3], col="#604c38", lwd = 3)
legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()

#Adjusted time series

country1 <- which(colnames(gdp_mat_original) == "AUT")
country2 <- which(colnames(gdp_mat_original) == "DEU")

pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], "_adj.pdf"),
    width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = dates, y = gdp_mat[, country1],
     ylim=c(min(gdp_mat[, country1], gdp_mat[, country2]),
            max(gdp_mat[, country1], gdp_mat[, country2] + 0.04)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
lines(x = dates, y = gdp_mat[, country2], col="#604c38")
legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()

#Original time series for Canada/USA

country1 <- which(colnames(gdp_mat_original) == "CAN")
country2 <- which(colnames(gdp_mat_original) == "USA")

pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], ".pdf"),
    width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = dates, y = gdp_mat_original[, country1],
     ylim=c(min(gdp_mat_original[, country1], gdp_mat_original[, country2]),
            max(gdp_mat_original[, country1], gdp_mat_original[, country2] + 0.03)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
lines(x = dates, y = gdp_mat_original[, country2], col="#604c38")
legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()

#Adjusted time series for Canada/USA

country1 <- which(colnames(gdp_mat) == "CAN")
country2 <- which(colnames(gdp_mat) == "USA")

pdf(paste0("plots/gdp_", countries[country1], "_", countries[country2], "_adj.pdf"),
    width=5.5, height=3, paper="special")

par(cex = 1, tck = -0.025)
par(mar = c(3, 3, 0.5, 0.5)) #Margins for each plot
par(oma = c(0.2, 0.2, 0.2, 0.2)) #Outer margins

plot(x = dates, y = gdp_mat[, country1],
     ylim=c(min(gdp_mat[, country1], gdp_mat[, country2]),
            max(gdp_mat[, country1], gdp_mat[, country2] + 0.03)), type="l",
     col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
lines(x = dates, y = gdp_mat[, country2], col="#604c38")
legend("topright", inset = 0.02, legend=c(country_names[country1], country_names[country2]),
       col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
dev.off()

#####################
#Estimating variance#
#####################

#Order selection
q <- 7:14
r <- 10:15
order_results <- c()

for (j in 1:n_ts){
  criterion_matrix <- expand.grid(q = q, r = r)
  
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
      AR.struc      <- estimate_lrv(data = gdp_mat[, j], q = criterion_matrix$q[[i]],
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
  cat("For the country ", colnames(gdp_mat)[j], " the results are as follows: ", max(criterion_matrix$FPE), " ", max(criterion_matrix$AIC), " ", max(criterion_matrix$AICC), " ", max(criterion_matrix$SIC), " ", max(criterion_matrix$HQ), " \n")
}


#Setting tuning parameters for testing
q     <- 12
r_bar <- 10


#Calculating each sigma_i separately
sigmahat_vector <- c()
for (i in 1:n_ts){
  AR.struc        <- estimate_lrv(data = gdp_mat[, i], q = q, r_bar = r_bar, p = order_results[i])
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector <- c(sigmahat_vector, sigma_hat_i)
}

#Calculating each sigma_i separately
sigmahat_vector2 <- c()
for (i in 1:n_ts){
  AR.struc        <- estimate_lrv(data = gdp_mat[, i], q = q, r_bar = r_bar, p = 1)
  sigma_hat_i     <- sqrt(AR.struc$lrv)
  sigmahat_vector2 <- c(sigmahat_vector2, sigma_hat_i)
}

#Constructing the grid
u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
h_grid <- seq(from = 5 / t_len, to = 1 / 4, by = 5 / t_len)
h_grid <- h_grid[h_grid > log(t_len) / t_len]
grid <- construct_grid(t = t_len, u_grid = u_grid, h_grid = h_grid)
grid_points <- seq(from = 1 / t_len, to = 1, by = 1 / t_len) #for plotting

#Calculating the statistic for real data
result <- multiscale_test(data = gdp_mat,
                          sigma_vec = sigmahat_vector2,
                          alpha = alpha,
                          n_ts = n_ts, grid = grid,
                          sim_runs = sim_runs, epidem = FALSE)

#Rename the countries for the plots
countries_names <- c("Australia", "Austria", "Canada", "Germany", "UK", "Japan", "USA")

gset  <- results$gset_with_values[[l]]
ticks <- c(0, 40, 80, 120)

source("functions.R")


for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  if (result$stat_pairwise[i, j] > result$quant){
    filename = paste0("plots/", countries[i], "_vs_", countries[j], ".pdf")
    pdf(filename, width=5.5, height=10.5, paper="special")
    layout(matrix(c(1, 2, 3),ncol=1), widths=c(2.2, 2.2, 2.2),
           heights=c(1.5, 1.5, 1.8), TRUE)
    
    #Setting the layout of the graphs
    
    par(cex = 1, tck = -0.025)
    par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
    par(oma = c(0.2, 1.5, 2, 0.2)) #Outer margins
    
    plot(gdp_mat[, i], ylim=c(min(gdp_mat[, i], gdp_mat[, j]),
                              max(gdp_mat[, i], gdp_mat[, j])),
         type="l", col="blue", ylab="", xlab="", xaxt = "n", mgp=c(1, 0.5, 0))
    lines(gdp_mat[, j], col="red")
    axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
         cex.axis = 0.95, mgp=c(1, 0.5, 0))
    
    title(main = "(a) adjusted GDP", font.main = 1, line = 0.5)
    legend("topright", inset = 0.02, legend=c(countries[i], countries[j]),
           col = c("blue", "red"), lty = 1, cex = 0.95, ncol = 1)
    
    par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot
    
    #Plotting the smoothed version of the time series that we have
    smoothed_i  <- mapply(nadaraya_watson_smoothing, grid_points,
                          MoreArgs = list(gdp_mat[, i], grid_points, bw = 5 / t_len))
    smoothed_j  <- mapply(nadaraya_watson_smoothing, grid_points,
                          MoreArgs = list(gdp_mat[, j], grid_points, bw = 5 / t_len))
    
    plot(smoothed_i, ylim=c(min(gdp_mat[, i], gdp_mat[, j]),
                            max(gdp_mat[, i], gdp_mat[, j])), type="l",
         col="black", ylab="", xlab = "", xaxt = "n", mgp=c(1,0.5,0))
    axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
         cex.axis = 0.95, mgp=c(1, 0.5, 0))
    title(main = "(b) smoothed curves from (a)", font.main = 1, line = 0.5)
    lines(smoothed_j, col="red")
    
    par(mar = c(2.7, 0.5, 3, 0)) #Margins for each plot
    gset    <- result$gset_with_values[[l]]
    a_t_set <- subset(gset, test == TRUE, select = c(u, h))
    if (nrow(a_t_set) > 0){
      p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * t_len + 0.5,
                            'endpoint' = (a_t_set$u + a_t_set$h) * t_len - 0.5, 'values' = 0)
      p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
      
      #Produce minimal intervals
      p_t_set2  <- compute_minimal_intervals(p_t_set)
      
      plot(NA, xlim=c(0, t_len),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab="", xaxt = "n",
           mgp=c(2, 0.5, 0), yaxt = "n")
      axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
           cex.axis = 0.95, mgp=c(1, 0.5, 0))
      title(main = "(c) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
      #title(xlab = "quarter", line = 1.7, cex.lab = 0.9)
      segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint, p_t_set2$values, lwd = 2)
      segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint, p_t_set$values, col = "gray")
    } else {
      #If there are no intervals where the test rejects, we produce empty plots
      plot(NA, xlim=c(0, t_len),  ylim = c(0, 1), xlab="", ylab = "", xaxt = "n",
           mgp=c(2,0.5,0), yaxt = "n")
      axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
           cex.axis = 0.95, mgp=c(1, 0.5, 0))
      title(main = "(c) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
      #title(xlab = "quarter", line = 1.7, cex.lab = 0.9)
    }
    mtext(paste0("Comparison of ", countries_names[i], " and ", countries_names[j]), side = 3,
          line = 0, outer = TRUE, font = 1, cex = 1.2)
    dev.off()
  }
}


#Producing plots for the talk

for (l in seq_len(nrow(result$ijset))){
  i <- result$ijset[l, 1]
  j <- result$ijset[l, 2]
  if (result$stat_pairwise[i, j] > result$quant){
    filename = paste0("plots/gdp_", countries[i], "_vs_", countries[j], "_talk.pdf")
    pdf(filename, width = 5, height = 6.5, paper="special")
    layout(matrix(c(1, 2), ncol=1), widths=c(2.4, 2.4),
           heights=c(1.5, 1.8), TRUE)
    
    #Setting the layout of the graphs
    par(cex = 1, tck = -0.025)
    par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
    par(oma = c(0.2, 1.5, 0.2, 0.2)) #Outer margins

    plot(x = dates, y = gdp_mat[, i], ylim=c(min(gdp_mat[, i], gdp_mat[, j]),
                                             max(gdp_mat[, i], gdp_mat[, j]) + 0.04),
         type="l", col="#EB811B", ylab="", xlab="", mgp=c(1, 0.5, 0))
    lines(x = dates, y = gdp_mat[, j], col="#604c38")

    title(main = "(a) augmented GDP series", font.main = 1, line = 0.5)
    legend("topright", inset = 0.02, legend = c(countries_names[i], countries_names[j]),
           col = c("#EB811B", "#604c38"), lty = 1, cex = 0.95, ncol = 1)
    
    par(mar = c(0.5, 0.5, 3, 0)) #Margins for each plot
    gset    <- result$gset_with_values[[l]]
    a_t_set <- subset(gset, test == TRUE, select = c(u, h))
    if (nrow(a_t_set) > 0){
      p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h) * t_len + 0.5,
                            'endpoint' = (a_t_set$u + a_t_set$h) * t_len - 0.5, 'values' = 0)
      p_t_set$values <- (1:nrow(p_t_set))/nrow(p_t_set)
      
      #Produce minimal intervals
      p_t_set2  <- compute_minimal_intervals(p_t_set)
      
      plot(NA, xlim=c(0, t_len),  ylim = c(0, 1 + 1 / nrow(p_t_set)), xlab="", xaxt = "n", yaxt = "n",
           mgp=c(2, 0.5, 0), yaxt = "n")
      title(main = "(b) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
      #title(xlab = "quarter", line = 1.7, cex.lab = 0.9)
      segments(p_t_set2$startpoint, p_t_set2$values, p_t_set2$endpoint, p_t_set2$values, lwd = 2)
      segments(p_t_set$startpoint, p_t_set$values, p_t_set$endpoint, p_t_set$values, col = "gray")
    } else {
      #If there are no intervals where the test rejects, we produce empty plots
      plot(NA, xlim=c(0, t_len),  ylim = c(0, 1), xlab="", ylab = "", xaxt = "n",
           mgp=c(2,0.5,0), yaxt = "n")
      axis(side = 1, at = ticks, labels = as.yearqtr(dates[ticks + 1], format = '%Y-Q%q'),
           cex.axis = 0.95, mgp=c(1, 0.5, 0))
      title(main = "(b) (minimal) intervals produced by our test", font.main = 1, line = 0.5)
    }
    dev.off()
  }
}



