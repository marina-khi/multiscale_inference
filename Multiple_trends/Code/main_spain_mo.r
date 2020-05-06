#########################################
#Analysis of Spain mortality time series#
#########################################
rm(list=ls())

library(Rcpp)
library(tidyr)

# The following file contains functions that are used to estimate the AR parameters, in particular, to compute
# the estimators of the coefficients, a_1, ... ,a_p, the estimator of the variance of the innovation term, sigma_eta^2 
# and the estimator of the long-run variance sigma^2.
source("functions/estimate_lrv.R")
source("functions/construct_grid.R")
source("functions/plot_SiZer_map.R")
source("functions/compute_quantiles.R")
source("functions/multiscale_test.R")
source("functions/compute_minimal_intervals.R")
source("functions/functions.R")
sourceCpp("functions/multiscale_stats.cpp")


alpha   <- 0.05 # alpha for calculating quantiles
SimRuns <- 100 # number of simulation runs to produce critical values

#Loading the real data for mortality in Spain
spain_mo <- read.csv("data/mortality_spain.csv", sep = ",", dec = ".", stringsAsFactors = FALSE, na.strings = "N/A")
spain_mo <- spain_mo[spain_mo$ambito == 'nacional' & spain_mo$cod_sexo == 'all' & spain_mo$cod_gedad =='all', 'defunciones_observadas']

Tlen      <- length(spain_mo)

#Order selection
source("functions/select_order.R")
orders <- select_order(as.matrix(spain_mo[1:700]), q = 30:50, r = 10:15)


#Setting tuning parameters for testing
p <- 9
q <- 35
r <- 15

ts_start = as.Date('2018/4/4')  
ts_end   = ts_start + Tlen - 1 #the last point of time series

#Estimate the long-run variance
sigmahat_vector_order_9 <- c()
AR.struc  <- AR_lrv(data=spain_mo[1:690], q = q, r.bar = r, p = p)
sigma_hat <- sqrt(AR.struc$lrv)
sigmahat_vector_order_9 <- c(sigmahat_vector_order_9, sigma_hat)

set.seed(1)
#Compute test results for the multiscale method
results <- multiscale_test(alpha = alpha, data = spain_mo, sigma_vec = sigmahat_vector_order_9, SimRuns = SimRuns, N_ts = 1, deriv_order = 1)
quant   <- results$quant


#Produce minimal intervals (Here - only the increases! But this is because we do not have decreases for our applications)
a_t_set <- subset(results$gset, test == 1, select = c(u, h, vals))
p_t_set <- data.frame('startpoint' = (a_t_set$u - a_t_set$h)*Tlen, 'endpoint' = (a_t_set$u + a_t_set$h)*Tlen, 'values' = a_t_set$vals)
p_t_set$endpoint[p_t_set$endpoint > ts_end] <- ts_end
#p_t_set <- subset(p_t_set, endpoint <= ts_end, select = c(startpoint, endpoint, values)) 
p_t_set <- compute_minimal_intervals(p_t_set)


pdf("plots/spain_mo.pdf", width=7, height=6.6, paper="special")
layout(matrix(c(1,2,3),ncol=1), widths=c(3,3,3), heights=c(1,0.8,1), TRUE) # Setting the layout of the graphs

#Parameters for plotting
grid_time <- seq(from = ts_start, to = ts_end, length.out = Tlen) #grid points for plotting 
p_t_set$plottingindex <- (1:nrow(p_t_set))/nrow(p_t_set)

par(cex = 1.1, tck = -0.025)
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

# Plotting the real data
plot(x = grid_time, xlim = c(ts_start, ts_end), y = spain_mo, type = "l", mgp=c(1,0.5,0), xaxs='i')
title(main = expression((a) ~ observed ~ temperature ~ time ~ series), line = 1)

par(mar = c(0.5, 0.5, 3.5, 0)) #Margins for each plot

#Plotting the minimal intervals. Do not have any negative minimal intervals, so plotting all (positive) ones
#ymaxlim = max(p_t_set$values)
#yminlim = min(min(p_t_set$values), quant.ms)
plot(NA, xlim=c(0,Tlen), xaxt = "n",  ylim = c(0, 1 + 1/nrow(p_t_set)), mgp=c(2,0.5,0))
title(main = expression((b) ~ minimal ~ intervals ~ produced ~ by ~ italic(T)[MS]), line = 1)
segments(p_t_set[['startpoint']], p_t_set$plottingindex, p_t_set$endpoint, p_t_set$plottingindex, lwd = 2)
#abline(h = quant.ms, lty = 2)

#SiZer 
par(mar = c(0.5, 0.5, 2, 0)) #Margins for each plot
plot_SiZer_map(sort(unique(results$gset[,1])), sort(unique(results$gset[,2])), results$test_matrix, plot.title = expression((c) ~ SiZer ~ map ~ 'for' ~ italic(T)[MS]))

axis_at = seq(4/Tlen, 724/Tlen, by = 30/Tlen)
axis_labels = seq(as.Date("2018/4/4"), as.Date("2020/4/4"), by = 30)
axis(1, at=axis_at, labels = axis_labels, mgp=c(1,0.5,0))

dev.off()