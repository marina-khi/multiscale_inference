source("Shape/functions_for_referees.R")
source("Shape/estimating_sigma_new.R")
dyn.load("Shape/C_code/psihat_statistic.dll")
dyn.load("Shape/C_code/psihat_statistic_without_lambda.dll")
source("Shape/C_code/psihat_statistic.R")
source("Shape/data_analysis.R")


#############################
#Point 3 in Referee Report 2#
#############################

# #Defining necessary parameters
# alpha <- 0.05 #alpha for calculating quantiles
# h     <- c(0.05, 0.1, 0.15, 0.2) #Different bandwidth for plotting. Number must be <=6 in order for the plot to be readable
# 
# test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported. 
# pdffilename_global = paste0("Paper/Plots/temperature_data_global.pdf") #Filename for the graph
# 
# #Loading the real data for global yearly temperature
# temperature_global  <- read.table("Shape/data/global_temp.txt", header = TRUE, skip = 16)
# yearly_tempr_global <- temperature_global[temperature_global$ANNUAL > -99, 'ANNUAL']
# T_tempr_global      <- length(yearly_tempr_global)
# 
# #Order selection for global
# q <- 10:20
# r <- 5:15
# criterion_matrix_global <- expand.grid(q = q, r = r)
# 
# criterion_matrix_global$FPE <- numeric(length = nrow(criterion_matrix_global))
# criterion_matrix_global$AIC <- numeric(length = nrow(criterion_matrix_global))
# criterion_matrix_global$AICC <- numeric(length = nrow(criterion_matrix_global))
# criterion_matrix_global$SIC <- numeric(length = nrow(criterion_matrix_global))
# criterion_matrix_global$HQ  <- numeric(length = nrow(criterion_matrix_global))
# 
# for (i in 1:nrow(criterion_matrix_global)){
#   FPE <- c()
#   AIC <- c()
#   AICC <- c()
#   SIC <- c()
#   HQ <- c()
#   different_orders <- (1:9)
#   for (order in different_orders){
#     sigma_eta_hat_method2_global <- estimating_variance_new(yearly_tempr_global, criterion_matrix_global$q[[i]], order, criterion_matrix_global$r[[i]])[[3]]
#     FPE <- c(FPE, (sigma_eta_hat_method2_global^2 * (T_tempr_global + order)) / (T_tempr_global - order))
#     AIC <- c(AIC, T_tempr_global * log(sigma_eta_hat_method2_global^2) + 2 * order)
#     AICC <- c(AICC, T_tempr_global * log(sigma_eta_hat_method2_global^2) + T_tempr_global* (1 + order / T_tempr_global)/(1 - (order +2)/T_tempr_global))
#     SIC <- c(SIC, log(sigma_eta_hat_method2_global^2) + order * log(T_tempr_global) / T_tempr_global)
#     HQ <- c(HQ, log(sigma_eta_hat_method2_global^2) + 2 * order * log(log(T_tempr_global)) / T_tempr_global)
#   }
#   criterion_matrix_global$FPE[[i]] <- which.min(FPE)
#   criterion_matrix_global$AIC[[i]] <- which.min(AIC)
#   criterion_matrix_global$AICC[[i]] <- which.min(AICC)
#   criterion_matrix_global$SIC[[i]] <- which.min(SIC)
#   criterion_matrix_global$HQ[[i]]  <- which.min(HQ)
# }
# 
# #Setting tuning parameters for testing global temperature
# p_global <- 4
# q_global <- 18
# r_global <- 10
# 
# #Data analysis
# 
# sigma_hat_global <- estimating_variance_new(yearly_tempr_global, q_global, order = p_global, r_global)[[1]]
# data_analysis_global(alpha, yearly_tempr_global, test_problem, sigma_hat_global, pdffilename_global)
# 
# #plotting fitted data
# pdf("Paper/Plots/fitting_different_AR.pdf")
# par(mfrow = c(9,2), cex = 0.8, tck = -0.025) #Setting the layout of the graphs
# par(mar = c(0, 0.5, 1.0, 0)) #Margins for each plot
# par(oma = c(1.5, 1.5, 0.2, 0.2)) #Outer margins
# 
# for (order_AR in 1:9){
#   coefficient <- estimating_variance_new(yearly_tempr_global, q_global, order = order_AR, r_global)
#   ts.sim <- arima.sim(list(c(order_AR, 0, 0), ar = coefficient[[2]]), n=T_tempr_global, rand.gen=function(n){rnorm(n, sd=coefficient[[3]])} )
#   
#   plot(ts.sim, ylab="", xlab = "", ylim = c(-1, 1), type = 'l',axes=FALSE, frame.plot=TRUE, cex = 1.2, tck = -0.025)
#   Axis(side=2, at  = c(-0.75,-0.25, 0.25, 0.75))
#   Axis(side=1, labels=FALSE)
#   plot(yearly_tempr_global - ts.sim, ylab="", xlab = "", ylim = c(-1, 1),yaxp  = c(-0.75, 0.75, 3), type = 'l', axes=FALSE, frame.plot=TRUE, cex = 1.2, tck = -0.025)
#   Axis(side=1, labels=FALSE)
#   Axis(side=2, labels = FALSE)
# }
# 
# dev.off()



#############################
#Point 7 in Referee Report 1#
#############################


source("Shape/functions_for_referees.R")

#Defining necessary parameters
alpha <- 0.05 #alpha for calculating quantiles

test_problem  <- "constant" #Only "zero" (H_0: m = 0) or "constant" (H_0: m = const) testing problems are currently supported.
pdffilename = paste0("JRSSB_submission/Plots/temperature_data_SiZer_maps.pdf") #Filename for the graph
colorlist   <- c('red', 'purple', 'blue', 'grey')

#Recoding testing problem and type of kernel estimator 
if (test_problem == "zero"){
  kernel_ind = 1
} else if (test_problem == "constant"){
  kernel_ind = 2
} else {
  print('Given testing problem is currently not supported')
}

#Loading the real data for yearly temperature in England
temperature  <- read.table("Shape/data/cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']
T_tempr      <- length(yearly_tempr)

different_i <- seq(from = 5/T_tempr, to = 1, by = 5/T_tempr)
different_h <- seq(from = 3/T_tempr, to = 1/4+3/T_tempr, by = 5/T_tempr)

#Setting tuning parameters for testing
p <- 2
q <- 25
r <- 10


#Data analysis
parameters <- estimating_variance_new(yearly_tempr, q, order = p, r) 

sigma_hat <- parameters[[1]]
a_hat     <- parameters[[2]]
sigma_eta <- parameters[[3]]

gamma = c()
for (k in 0:(T_tempr-1)){
  gamma_temp <- gamma
  gamma = c(gamma, autocovariance_function_AR2(k, a_hat[[1]], a_hat[[2]], sigma_eta, gamma_temp))  #Note that gamma[i] := \gamma(i-1)
}
rm(gamma_temp)


#Calculating \Var(\bar{Y}) based on the true values of gamma(k)
true_var <- gamma[1] / T_tempr
for (k in 1:(T_tempr-1)){true_var = true_var + (2/T_tempr) * (1 - k/T_tempr) * gamma[k+1]}

T_star   <- gamma[1]/true_var

SiZer_matrix      <- calculating_SiZer_matrix(different_i, different_h, T_tempr, T_star, alpha, gamma)  

g_t_set <- psihat_statistic(yearly_tempr, SiZer_matrix, kernel_ind, sigma_hat)[[1]]
gaussian_quantile <- calculating_gaussian_quantile(T_tempr, SiZer_matrix, "data", kernel_ind, alpha)
g_t_set$gaussian_quantile <- gaussian_quantile

g_t_set_temp <- NULL

for (bandwidth in different_h){
  SiZer_matrix_temp                           <- subset(SiZer_matrix, h == bandwidth, select = c(u, h, lambda))
  if (nrow(SiZer_matrix_temp)>0){
    gaussian_quantile_rowwise                   <- calculating_gaussian_quantile(T_tempr, SiZer_matrix_temp, paste0("data_h_", bandwidth*1000), kernel_ind, alpha)
    g_t_set_temp_temp                           <- psihat_statistic(yearly_tempr, SiZer_matrix_temp, kernel_ind, sigma_hat)[[1]]
    g_t_set_temp_temp$gaussian_quantile_rowwise <- gaussian_quantile_rowwise
    g_t_set_temp                                <- rbind(g_t_set_temp, g_t_set_temp_temp)
  }
}

g_t_set_rowwise <- merge(SiZer_matrix, g_t_set_temp, by = c('h', 'u', 'lambda'))


for (row in 1:nrow(g_t_set)){
  i              = g_t_set[row, 'u']
  h              = g_t_set[row, 'h']
  q_h            = g_t_set[row, 'q_h']
  sd_m_hat_prime = g_t_set[row, 'sd']
  
  XtWX_inverse_XtW = g_t_set$XtWX_inv_XtW[[row]]
  
  if (!is.null(XtWX_inverse_XtW)) {
    m_hat_prime <- (XtWX_inverse_XtW %*% yearly_tempr)[2]
    if (m_hat_prime - q_h * sd_m_hat_prime > 0){
      g_t_set$results_their[[row]] = 1
    } else if (m_hat_prime + q_h * sd_m_hat_prime < 0) {
      g_t_set$results_their[[row]] = -1
    } else {
      g_t_set$results_their[[row]] = 0
    }
  } else {
    g_t_set$results_their[[row]] = 2
  }
  
  if (g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + g_t_set[row, 'gaussian_quantile']){
    g_t_set$results_our[[row]] = 1
  } else if (-g_t_set[row, 'values_with_sign'] > g_t_set[row, 'lambda'] + g_t_set[row, 'gaussian_quantile']){
    g_t_set$results_our[[row]] = -1
  } else {
    g_t_set$results_our[[row]] = 0
  }
  
  if (g_t_set_rowwise[row, 'values_with_sign'] > g_t_set_rowwise[row, 'lambda'] + g_t_set_rowwise[row, 'gaussian_quantile_rowwise']){
    g_t_set_rowwise$results_our_rowwise[[row]] = 1
  } else if (-g_t_set_rowwise[row, 'values_with_sign'] > g_t_set_rowwise[row, 'lambda'] + g_t_set_rowwise[row, 'gaussian_quantile_rowwise']){
    g_t_set_rowwise$results_our_rowwise[[row]] = -1
  } else {
    g_t_set_rowwise$results_our_rowwise[[row]] = 0
  }
  
}
result_SiZer <- subset(g_t_set, select = c(u, h, results_their))
result_our   <- subset(g_t_set, select = c(u, h, results_our))
result_our_rowwise        <- subset(g_t_set_rowwise, select = c(u, h, results_our_rowwise))

pdf(pdffilename, width = 7, height = 3)
       
par(mfrow = c(1, 3), cex = 0.5, tck = -0.025) #Setting the layout of the graphs
par(mar = c(1, 2, 3, 0.5)) #Margins for each plot
par(oma = c(1.5, 1.5, 3, 0.2)) #Outer margins

plot.SiZer(result_SiZer$results_their, different_i, different_h, ylab=expression(log[10](h)),
           colorlist=colorlist, title = "SiZer results for the data")

plot.SiZer(result_our$results_our, different_i, different_h, ylab=expression(log[10](h)),
           colorlist=colorlist, title = "Our global results for the data")

plot.SiZer(result_our_rowwise$results_our_rowwise, different_i, different_h, ylab=expression(log[10](h)),
           colorlist=colorlist, title = "Our rowwise results for the data")
dev.off()
