rm(list=ls())

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/grid_construction.r")
source("functions/kernel_weights.r")
source("functions/multiscale_statistics.r")
source("functions/multiscale_quantiles.r")
source("functions/multiscale_testing.r")
source("functions/long_run_variance.r")
source("functions/sim.r")
source("functions/minimal_intervals.r")
source("functions/SiZer_functions.r")
source("functions/data_analysis.r")

if (Sys.info()[[1]] == 'Windows'){
  dyn.load("functions/SiZer_functions.dll")
  dyn.load("functions/kernel_weights.dll")
} else {
  dyn.load("functions/SiZer_functions.so")
  dyn.load("functions/kernel_weights.so") 
}


###############################################
#Analysis of UK annual temperature time series#
###############################################

alpha   <- 0.05     #alpha for calculating quantiles
h       <- c(0.05, 0.1, 0.15, 0.2) #Different bandwidth for plotting. Number must be <=6 in order for the plot to be readable
SimRuns <- 5000        # number of simulation runs to produce critical values


#Loading the real data for yearly temperature in England
temperature  <- read.table("data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
T_tempr      <- length(yearly_tempr)

#Order selection for England
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

write.csv(criterion_matrix,"plots/criterion_matrix.csv", row.names = FALSE)

#Setting tuning parameters for testing
p <- 2
q <- 15
r <- 10

# pdf("Plots/temperature.pdf", width=10, height=3, paper="special")
# par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
# par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
# data <- ts(yearly_tempr, start=1659, end=2017, frequency=1)
# plot(data, ylab="", xlab = "", yaxp  = c(7, 11, 4), xaxp = c(1675, 2025, 7), type = 'l', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
# dev.off()


data_analysis(data = yearly_tempr, filename_table = "plots/minimal_intervals_UK.tex",
              filename_plot = "plots/UK_temperature.pdf", axis_at = seq(17/T_tempr, 367/T_tempr, by = 50/T_tempr),
              axis_labels = seq(1675, 2025, by = 50), xaxp = c(1675, 2025, 7), yaxp = c(1.75, 2.5, 3),
              ts_start= 1659, ts_end = 2017, alpha = alpha, SimRuns = SimRuns, order = p, q = q, r.bar = r)
  

#########################################
#Analysis of the global temperature data#
#########################################

#Loading the real data for global yearly temperature
temperature_global  <- read.table("data/global_temp.txt", header = TRUE, skip = 16)
yearly_tempr_global <- temperature_global[temperature_global$ANNUAL > -99, 'ANNUAL']
T_tempr_global      <- length(yearly_tempr_global)

#Order selection for global
q <- 10:35
r <- 5:15
criterion_matrix_global <- expand.grid(q = q, r = r)

criterion_matrix_global$FPE <- numeric(length = nrow(criterion_matrix_global))
criterion_matrix_global$AIC <- numeric(length = nrow(criterion_matrix_global))
criterion_matrix_global$AICC <- numeric(length = nrow(criterion_matrix_global))
criterion_matrix_global$SIC <- numeric(length = nrow(criterion_matrix_global))
criterion_matrix_global$HQ  <- numeric(length = nrow(criterion_matrix_global))

for (i in 1:nrow(criterion_matrix_global)){
  FPE <- c()
  AIC <- c()
  AICC <- c()
  SIC <- c()
  HQ <- c()
  different_orders <- (1:9)
  for (order in different_orders){
    AR.struc      <- AR_lrv(data=yearly_tempr_global, q=criterion_matrix_global$q[[i]], r.bar=criterion_matrix_global$r[[i]], p=order)
    sigma_eta_hat <- sqrt(AR.struc$vareta)
    FPE <- c(FPE, (sigma_eta_hat^2 * (T_tempr_global + order)) / (T_tempr_global - order))
    AIC <- c(AIC, T_tempr_global * log(sigma_eta_hat^2) + 2 * order)
    AICC <- c(AICC, T_tempr_global * log(sigma_eta_hat^2) + T_tempr_global* (1 + order / T_tempr_global)/(1 - (order +2)/T_tempr_global))
    SIC <- c(SIC, log(sigma_eta_hat^2) + order * log(T_tempr_global) / T_tempr_global)
    HQ <- c(HQ, log(sigma_eta_hat^2) + 2 * order * log(log(T_tempr_global)) / T_tempr_global)
  }
  criterion_matrix_global$FPE[[i]] <- which.min(FPE)
  criterion_matrix_global$AIC[[i]] <- which.min(AIC)
  criterion_matrix_global$AICC[[i]] <- which.min(AICC)
  criterion_matrix_global$SIC[[i]] <- which.min(SIC)
  criterion_matrix_global$HQ[[i]]  <- which.min(HQ)
}

write.csv(criterion_matrix_global,"plots/criterion_matrix_global.csv", row.names = FALSE)



#Setting tuning parameters for testing global temperature
p_global <- 4
q_global <- 15
r_global <- 10

data_analysis(data = yearly_tempr_global, filename_table = "plots/minimal_intervals_global.tex",
              filename_plot = "plots/global_temperature.pdf", axis_at = seq(25/T_tempr_global, 125/T_tempr_global, by = 50/T_tempr_global),
              axis_labels = seq(1875, 1975, by = 50), xaxp = c(1875, 1975, 2), yaxp = c(1.75, 3.75, 4), 
              ts_start= 1850, ts_end = 2014, alpha = alpha,
              SimRuns = SimRuns, order = p_global, q = q_global, r.bar = r_global,
              plot_SiZer = "no", sigma_supplied = 'yes', sigma = sqrt(0.01558))