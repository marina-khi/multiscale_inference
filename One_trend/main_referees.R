library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
source("Shape/functions.R")
source("Shape/estimating_sigma_new.R")
dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_ll.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/data_analysis.R")


###############################
#Defining necessary parameters#
###############################
alpha <- 0.05 #alpha for calculating quantiles
h     <- c(0.05, 0.1, 0.15, 0.2) #Different bandwidth for plotting. Number must be <=6 in order for the plot to be readable

test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 
pdffilename_global = paste0("Paper/Plots/temperature_data_global.pdf") #Filename for the graph

#####################################################
#Loading the real data for global yearly temperature#
#####################################################
temperature_global  <- read.table("Shape/data/global_temp.txt", header = TRUE, skip = 16)
yearly_tempr_global <- temperature_global[temperature_global$ANNUAL > -99, 'ANNUAL']
T_tempr_global      <- length(yearly_tempr_global)


############################
#Order selection for global#
############################
q <- 10:20
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
    sigma_eta_hat_method2_global <- estimating_variance_new(yearly_tempr_global, criterion_matrix_global$q[[i]], order, criterion_matrix_global$r[[i]])[[3]]
    FPE <- c(FPE, (sigma_eta_hat_method2_global^2 * (T_tempr_global + order)) / (T_tempr_global - order))
    AIC <- c(AIC, T_tempr_global * log(sigma_eta_hat_method2_global^2) + 2 * order)
    AICC <- c(AICC, T_tempr_global * log(sigma_eta_hat_method2_global^2) + T_tempr_global* (1 + order / T_tempr_global)/(1 - (order +2)/T_tempr_global))
    SIC <- c(SIC, log(sigma_eta_hat_method2_global^2) + order * log(T_tempr_global) / T_tempr_global)
    HQ <- c(HQ, log(sigma_eta_hat_method2_global^2) + 2 * order * log(log(T_tempr_global)) / T_tempr_global)
  }
  criterion_matrix_global$FPE[[i]] <- which.min(FPE)
  criterion_matrix_global$AIC[[i]] <- which.min(AIC)
  
  criterion_matrix_global$AICC[[i]] <- which.min(AICC)
  criterion_matrix_global$SIC[[i]] <- which.min(SIC)
  criterion_matrix_global$HQ[[i]]  <- which.min(HQ)
}


##########################################################
#Setting tuning parameters for testing global temperature#
##########################################################
p_global <- 4
q_global <- 18
r_global <- 10


###############
#Data analysis#
###############

sigma_hat_global <- estimating_variance_new(yearly_tempr_global, q_global, order = p_global, r_global)[[1]]
data_analysis_global(alpha, yearly_tempr_global, test_problem, sigma_hat_global, pdffilename_global)



pdf("Paper/Plots/fitting_different_AR.pdf")
par(mfrow = c(9,2), cex = 0.8, tck = -0.025) #Setting the layout of the graphs
par(mar = c(0, 0.5, 1.0, 0)) #Margins for each plot
par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins

for (order_AR in 1:9){
  coefficient <- estimating_variance_new(yearly_tempr_global, q_global, order = order_AR, r_global)
  ts.sim <- arima.sim(list(c(order_AR, 0, 0), ar = coefficient[[2]]), n=T_tempr_global, rand.gen=function(n){rnorm(n, sd=coefficient[[3]])} )
  
  plot(ts.sim, ylab="", xlab = "", ylim = c(-1, 1), type = 'l',axes=FALSE, frame.plot=TRUE, cex = 1.2, tck = -0.025)
  Axis(side=2, at  = c(-0.75,-0.25, 0.25, 0.75))
  Axis(side=1, labels=FALSE)
  plot(yearly_tempr_global - ts.sim, ylab="", xlab = "", ylim = c(-1, 1),yaxp  = c(-0.75, 0.75, 3), type = 'l', axes=FALSE, frame.plot=TRUE, cex = 1.2, tck = -0.025)
  Axis(side=1, labels=FALSE)
  Axis(side=2, labels = FALSE)
}

dev.off()