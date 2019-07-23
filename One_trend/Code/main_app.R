# This is the main file for the analysis of both applications which is reported in Section 6.
rm(list=ls())

library(Rcpp)
library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

# The following file contains functions that are used to estimate the AR parameters, in particular, to compute
# the estimators of the coefficients, a_1, ... ,a_p, the estimator of the variance of the innovation term, sigma_eta^2 
# and the estimator of the long-run variance sigma^2.
source("functions/long_run_variance.r")

# The following file contains one main function that performs different multiscale tests on the data to determine where 
# the trend is increasing or decreasing. It also produces all the necessary plots such as SiZer map or the plot of minimal intervals.
# All the necessary arguments for this function are described in detail in the file.
source("functions/AnalyzeData.r")


###############################################
#Analysis of UK annual temperature time series#
###############################################

alpha   <- 0.05 # alpha for calculating quantiles
SimRuns <- 5000 # number of simulation runs to produce critical values

#Loading the real data for yearly temperature in England
temperature  <- read.table("data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
T_tempr      <- length(yearly_tempr)

#Order selection
q <- 10:35
r <- 10:15
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
    AR.struc      <- AR_lrv(data=yearly_tempr, q=criterion_matrix$q[[i]], r.bar=criterion_matrix$r[[i]], p=order)
    sigma_eta_hat <- sqrt(AR.struc$vareta)
    FPE <- c(FPE, (sigma_eta_hat^2 * (T_tempr + order)) / (T_tempr - order))
    AIC <- c(AIC, T_tempr * log(sigma_eta_hat^2) + 2 * order)
    AICC <- c(AICC, T_tempr * log(sigma_eta_hat^2) + T_tempr* (1 + order / T_tempr)/(1 - (order +2)/T_tempr))
    SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(T_tempr) / T_tempr)
    HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(T_tempr)) / T_tempr)
  }
  criterion_matrix$FPE[[i]]  <- which.min(FPE)
  criterion_matrix$AIC[[i]]  <- which.min(AIC)
  criterion_matrix$AICC[[i]] <- which.min(AICC)
  criterion_matrix$SIC[[i]]  <- which.min(SIC)
  criterion_matrix$HQ[[i]]   <- which.min(HQ)
}

write.csv(criterion_matrix,"plots/criterion_matrix_UK.csv", row.names = FALSE)

#Setting tuning parameters for testing
p <- 2
q <- 25
r <- 10

#And finally testing itself!
AnalyzeData(data = yearly_tempr, ts_start= 1659, filename_table = "plots/minimal_intervals_UK.tex",
              filename_plot = "plots/temperature_UK.pdf", axis_at = seq(17/T_tempr, 367/T_tempr, by = 50/T_tempr),
              axis_labels = seq(1675, 2025, by = 50), xaxp = c(1675, 2025, 7), yaxp = c(1.9, 2.2, 3),
              alpha = alpha, SimRuns = SimRuns, order = p, q = q, r.bar = r)
  

#########################################
#Analysis of the global temperature data#
#########################################
alpha   <- 0.05 # alpha for calculating quantiles
SimRuns <- 5000 # number of simulation runs to produce critical values

#Loading the real data for global yearly temperature
temperature_global  <- read.table("data/global_temp.txt", header = TRUE, skip = 16)
yearly_tempr_global <- temperature_global[(temperature_global$ANNUAL > -99) &(temperature_global$YEAR <= 1998), 'ANNUAL']
T_tempr_global      <- length(yearly_tempr_global)

#And finally testing itself!
AnalyzeData(data = yearly_tempr_global, ts_start= 1850, filename_table = "plots/minimal_intervals_global.tex",
              filename_plot = "plots/temperature_global.pdf", axis_at = seq(25/T_tempr_global, 125/T_tempr_global, by = 50/T_tempr_global),
              axis_labels = seq(1875, 1975, by = 50), xaxp = c(1875, 1975, 2), yaxp = c(1.8, 2.4, 2), 
              alpha = alpha,  SimRuns = SimRuns, plot_SiZer = "no", sigma_supplied = 'yes', sigma = sqrt(0.01558))