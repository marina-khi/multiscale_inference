rm(list=ls())

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

source("functions/grid_construction.r")
dyn.load("functions/kernel_weights.dll")
source("functions/kernel_weights.r")
source("functions/multiscale_statistics.r")
source("functions/multiscale_quantiles.r")
source("functions/multiscale_testing.r")
source("functions/long_run_variance.r")
source("functions/sim.r")
source("functions/minimal_intervals.r")
source("functions/SiZer_functions.r")
dyn.load("functions/SiZer_functions.dll")


###############################
#Defining necessary parameters#
###############################

alpha   <- 0.05 #alpha for calculating quantiles
h       <- c(0.05, 0.1, 0.15, 0.2) #Different bandwidth for plotting. Number must be <=6 in order for the plot to be readable
SimRuns <- 5000


#########################################################
#Loading the real data for yearly temperature in England#
#########################################################

temperature  <- read.table("data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
T_tempr      <- length(yearly_tempr)


#############################
#Order selection for England#
#############################
q <- 10:35
r <- 10:15
criterion_matrix <- expand.grid(q = q, r = r)

criterion_matrix$FPE <- numeric(length = nrow(criterion_matrix))
criterion_matrix$AIC <- numeric(length = nrow(criterion_matrix))
criterion_matrix$AICC <- numeric(length = nrow(criterion_matrix))

criterion_matrix$SIC <- numeric(length = nrow(criterion_matrix))
criterion_matrix$HQ  <- numeric(length = nrow(criterion_matrix))

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
  criterion_matrix$FPE[[i]] <- which.min(FPE)
  criterion_matrix$AIC[[i]] <- which.min(AIC)
  criterion_matrix$AICC[[i]] <- which.min(AICC)
  criterion_matrix$SIC[[i]] <- which.min(SIC)
  criterion_matrix$HQ[[i]]  <- which.min(HQ)
}

#######################################
#Setting tuning parameters for testing#
#######################################
p <- 2
q <- 15
r <- 10


###################################################################
#Calculating smoothed curve for the data using Epanechnikov kernel#
###################################################################
# grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for plotting and estimating
# 
# for (i in 1:length(h)){
#   smoothed_curve <- mapply(local_linear_smoothing, grid_points, MoreArgs = list(yearly_tempr, grid_points, h[i]))
#   cat("End point:", smoothed_curve[T_tempr], "\n")
# }
#
# pdf("Plots/temperature.pdf", width=10, height=3, paper="special")
# par(mar = c(0, 0.5, 0, 0)) #Margins for each plot
# par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
# data <- ts(yearly_tempr, start=1659, end=2017, frequency=1)
# plot(data, ylab="", xlab = "", yaxp  = c(7, 11, 4), xaxp = c(1675, 2025, 7), type = 'l', mgp=c(2,0.5,0), cex = 1.2, tck = -0.025)
# dev.off()
# 

###############
#Data analysis#
###############
AR.struc  <- AR_lrv(data=yearly_tempr, q=q, r.bar=r, p=p)
sigma_hat <- sqrt(AR.struc$lrv)


#Construct grid
grid    <- grid_construction(T_tempr)
gset    <- grid$gset
u.grid  <- sort(unique(gset[,1]))
h.grid  <- sort(unique(gset[,2]))
correct <- sqrt(2*log(1/(2*gset[,2])))


# Compute kernel weights and critical value for multiscale test
wghts  <- kernel_weights(T=T_tempr, grid=grid)
stats  <- multiscale_statistics(data=yearly_tempr, weights=wghts, sigmahat=sigma_hat, grid=grid)
vals   <- stats$values

filename = paste0("quantiles/distr_T_", T_tempr,".RData")
if(!file.exists(filename)) {
  quants <- multiscale_quantiles(T=T_tempr, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
  save(quants, file = filename)
} else {
  load(filename)
}

 


# Select (1-alpha) quantile of the multiscale statistic under the null
probs.ms <- as.vector(quants$quant_ms[1,])
quant.ms <- as.vector(quants$quant_ms[2,])

if(sum(probs.ms == (1-alpha)) == 0){
  pos.ms <- which.min(abs(probs.ms-(1-alpha)))
} else {
  pos.ms <- which.max(probs.ms == (1-alpha)) 
}


quant.ms    <- quant.ms[pos.ms]


#Compute test results
vals2            <- abs(vals) - correct
gset_result      <- cbind(gset, vals, vals2)
gset_result$test <- (gset_result$vals2 > quant.ms) * sign(gset_result$vals)

a_t_set <- subset(gset_result, test == 1, select = c(u, h, vals2))
p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*T_tempr + 1659, 'endpoint' = (a_t_set$u + a_t_set$h)*T_tempr + 1659, 'values' = a_t_set$vals2)
p_t_set <- subset(p_t_set, endpoint <= 2017, select = c(startpoint, endpoint, values)) 
p_t_set <- choosing_minimal_intervals(p_t_set)

print.xtable(xtable(subset(p_t_set, select = c(startpoint, endpoint)), digits = c(0)), type="latex", file="plots/minimal_intervals.tex")


gamma = c()
for (k in 0:(T_tempr-1)){
  gamma_temp <- gamma
  gamma = c(gamma, autocovariance_function_AR2(k, AR.struc$ahat[[1]], AR.struc$ahat[[2]], sqrt(AR.struc$vareta), gamma_temp))  #Note that gamma[i] := \gamma(i-1)
}
rm(gamma_temp)

sizer.wghts   <- SiZer_weights(T=T_tempr, grid=grid)
sizer.std     <- SiZer_std(weights=sizer.wghts, autocov=gamma, T_tempr)
sizer.quants  <- SiZer_quantiles(alpha=alpha, T=T_tempr, grid=grid, autocov=gamma)
SiZer.values  <- sizer.wghts %*% yearly_tempr
sizer.vals    <- as.vector(SiZer.values)
SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)

test.res      <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)



#Parameters for plotting
grid_points <- seq(from = 1/T_tempr, to = 1, length.out = T_tempr) #grid points for estimating
grid_time <- seq(from = 1659, to = 2017, length.out = T_tempr) #grid points for plotting 

pdffilename <- paste0("plots/application.pdf")
pdf(pdffilename, width=8, height=10, paper="special")

par(mfrow = c(4,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

# Plotting the real data
plot(grid_time, xlim = c(1659, 2019), yearly_tempr, type = "l", mgp=c(1,0.5,0), xaxp = c(1675, 2025, 7)) 


#Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
ymaxlim = max(p_t_set$values)
yminlim = min(p_t_set$values)
plot(NA, xlim=c(1659,2019), xaxt = "n",  ylim = c(yminlim - 0.2, ymaxlim + 0.2), yaxp  = c(1.75, 2.5, 3), mgp=c(2,0.5,0))
segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set$endpoint, p_t_set[['values']])
abline(h = quant.ms, lty = 2)

SiZermap(u.grid, h.grid, test.res$test_ms, plot.title = expression(italic(T)[MS]))
axis(1, at=seq(17/T_tempr, 367/T_tempr, by = 50/T_tempr), labels = FALSE, mgp=c(1,0.5,0))
    
SiZermap(u.grid, h.grid, SiZer_results$test, plot.title = expression(italic(T)[SiZer]))
    
axis(1, at=seq(17/T_tempr, 367/T_tempr, by = 50/T_tempr), labels = seq(1675, 2025, by = 50), mgp=c(1,0.5,0))
dev.off()



#############################
#Point 3 in Referee Report 2#
#############################

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

#Setting tuning parameters for testing global temperature
p_global <- 4
q_global <- 15
r_global <- 10

#Data analysis
AR.struc  <- AR_lrv(data=yearly_tempr_global, q=q_global, r.bar=r_global, p=p_global)
sigma_hat <- sqrt(AR.struc$lrv)


#Construct grid
grid    <- grid_construction(T_tempr_global)
gset    <- grid$gset
u.grid  <- sort(unique(gset[,1]))
h.grid  <- sort(unique(gset[,2]))
correct <- sqrt(2*log(1/(2*gset[,2])))


# Compute kernel weights and critical value for multiscale test
wghts  <- kernel_weights(T=T_tempr_global, grid=grid)
stats  <- multiscale_statistics(data=yearly_tempr_global, weights=wghts, sigmahat=sigma_hat, grid=grid)
vals   <- stats$values

filename = paste0("quantiles/distr_T_", T_tempr_global,".RData")
if(!file.exists(filename)) {
  quants <- multiscale_quantiles(T=T_tempr_global, grid=grid, weights=wghts, kappa=0.1, SimRuns=SimRuns)
  save(quants, file = filename)
} else {
  load(filename)
}




# Select (1-alpha) quantile of the multiscale statistic under the null
probs.ms <- as.vector(quants$quant_ms[1,])
quant.ms <- as.vector(quants$quant_ms[2,])

if(sum(probs.ms == (1-alpha)) == 0){
  pos.ms <- which.min(abs(probs.ms-(1-alpha)))
} else {
  pos.ms <- which.max(probs.ms == (1-alpha)) 
}


quant.ms    <- quant.ms[pos.ms]


#Compute test results
vals2            <- abs(vals) - correct
gset_result      <- cbind(gset, vals, vals2)
gset_result$test <- (gset_result$vals2 > quant.ms) * sign(gset_result$vals)

a_t_set <- subset(gset_result, test == 1, select = c(u, h, vals2))
p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*T_tempr_global + 1850, 'endpoint' = (a_t_set$u + a_t_set$h)*T_tempr_global + 1850, 'values' = a_t_set$vals2)
#p_t_set <- subset(p_t_set, endpoint <= 2015, select = c(startpoint, endpoint, values)) 
p_t_set <- choosing_minimal_intervals(p_t_set)

print.xtable(xtable(subset(p_t_set, select = c(startpoint, endpoint)), digits = c(0)), type="latex", file="plots/minimal_intervals_global.tex")


gamma <- AR_acf(AR.struc$ahat, AR.struc$vareta, T_tempr_global)

sizer.wghts   <- SiZer_weights(T=T_tempr_global, grid=grid)
sizer.std     <- SiZer_std(weights=sizer.wghts, autocov=gamma, T_tempr_global)
sizer.quants  <- SiZer_quantiles(alpha=alpha, T=T_tempr_global, grid=grid, autocov=gamma)
SiZer.values  <- sizer.wghts %*% yearly_tempr_global
sizer.vals    <- as.vector(SiZer.values)
SiZer_results <- SiZer_test(values=sizer.vals, std.devs=sizer.std, quants=sizer.quants, grid=grid)

test.res      <- multiscale_testing(alpha=alpha, quantiles=quants, values=vals, grid=grid)



#Parameters for plotting
grid_points <- seq(from = 1/T_tempr_global, to = 1, length.out = T_tempr_global) #grid points for estimating
grid_time <- seq(from = 1850, to = 2014, length.out = T_tempr_global) #grid points for plotting 

pdffilename <- paste0("plots/application_global.pdf")
pdf(pdffilename, width=8, height=10, paper="special")

par(mfrow = c(4,1), cex = 1.1, tck = -0.025) #Setting the layout of the graphs
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

# Plotting the real data
plot(grid_time, xlim = c(1850, 2014), yearly_tempr_global, type = "l", mgp=c(1,0.5,0), xaxp = c(1900, 2000, 4)) 


#Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
ymaxlim = max(p_t_set$values)
yminlim = min(p_t_set$values)
plot(NA, xlim = c(1850, 2014), xaxt = "n",  ylim = c(yminlim - 0.2, ymaxlim + 0.2), mgp=c(2,0.5,0))
segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set$endpoint, p_t_set[['values']])
abline(h = quant.ms, lty = 2)

SiZermap(u.grid, h.grid, test.res$test_ms, plot.title = expression(italic(T)[MS]))
#axis(1, at=seq(17/T_tempr, 367/T_tempr, by = 50/T_tempr), labels = FALSE, mgp=c(1,0.5,0))

SiZermap(u.grid, h.grid, SiZer_results$test, plot.title = expression(italic(T)[SiZer]))

#axis(1, at=seq(17/T_tempr, 367/T_tempr, by = 50/T_tempr), labels = seq(1675, 2025, by = 50), mgp=c(1,0.5,0))
dev.off()
