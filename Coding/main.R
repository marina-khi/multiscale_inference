source("functions.R")
source("estimating_sigma.R")
dyn.load("psihat_statistic.dll")
source("psihat_statistic.R")

##############################
#Defining necessary constants#
##############################

sigma <- 2 #square root of long run-variance
T<-1000 #length 
alpha <-0.05
noise_to_signal <- 1
p <- 1 #Order of AR(p) process of the error terms
sigmahat <- sigma #Here will be the estimate of the square root of the long-run variance sigma^2
kernel_f = "biweight" #Only "epanechnikov" and "biweight" kernel functions are currently supported

if (kernel_f == "epanechnikov"){
  kernel_ind = 1
} else if (kernel_f == "biweight"){
  kernel_ind = 2
} else {
  print('Currently only Epanechnikov and Biweight kernel functions are supported')
}

#####################################
#Calculating gaussian quantile for T#
#####################################

g_t_set <- creating_g_set(T)
gaussian_quantile <- calculating_gaussian_quantile(T, g_t_set, sigmahat, alpha)


#################################################
#Defining the data for simulations Y = m + noise#
#################################################
set.seed(16) #For reproducibility

#Adding a function that is 0 on the first half and linear on the second half
a <- sqrt(6*sigma*sigma/noise_to_signal)#This is the value of m(1), it depends on Var(e) and Noise to Signal Ratio
b <-0.5 * T/a #Constant needed just for calculation
m <- numeric(T)
for (i in 1:T) {
  if (i/T < 0.5) {m[i] <- 0} else {m[i] <- (i - 0.5*T)/b}
} 
y_data_1 <- m + rnorm(T, 0, sigma)

#Adding a function that is const on the whole interval
const <- sqrt(sigma*sigma/noise_to_signal)
y_data_2 <- const + rnorm(T, 0, sigma)

#Not adding anything, here the null hypothesis is true
y_data_3 <- rnorm(T, 0, sigma)


###########################################
#Calculating the statistic for simulations#
###########################################


#Function that takes the data as argument and plots the regions where
#the corrected test statistic is bigger than the critical value of the 
#gaussian version.
#Returns the set with regions and value of the statistic in this regions.
plotting_all_rejected_intervals <-function(data, g_t_set, kernel_ind, sigm, dir, plotname){
  
  #Calculating our statistic
  result <-psihat_statistic(data, g_t_set, kernel_ind, sigm)
  g_t_set_with_values <- result[[1]]
  psihat_statistic_value <- result[[2]]
  
  #And now the testing itself
  if (psihat_statistic_value > gaussian_quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")

    #The collection of intervals where the corrected test statistic lies above the critical value (as in (2.6))
    a_t_set <- subset(g_t_set_with_values, values > gaussian_quantile, select = c(u, h, values))
    p_t_set <- data.frame('startpoint' = a_t_set$u - a_t_set$h,
                          'endpoint' = a_t_set$u + a_t_set$h,
                          'values' = a_t_set$values)
    
    p_t_set <- choosing_minimal_intervals(p_t_set)
    
    a_t_set_plus <- subset(g_t_set_with_values, values_with_sign > gaussian_quantile + lambda, select = c(u, h, values_with_sign))
    p_t_set_plus <- data.frame('startpoint' = a_t_set_plus$u - a_t_set_plus$h,
                          'endpoint' = a_t_set_plus$u + a_t_set_plus$h,
                          'values' = a_t_set_plus$values_with_sign)
    
    p_t_set_plus <- choosing_minimal_intervals(p_t_set_plus)
    
    a_t_set_minus <- subset(g_t_set_with_values, -values_with_sign > gaussian_quantile + lambda, select = c(u, h, values_with_sign))
    p_t_set_minus <- data.frame('startpoint' = a_t_set_minus$u - a_t_set_minus$h,
                               'endpoint' = a_t_set_minus$u + a_t_set_minus$h,
                               'values' = a_t_set_minus$values_with_sign)
    
    p_t_set_minus <- choosing_minimal_intervals(p_t_set_minus)
      
    #The plotting itself
    plotname1 = paste(dir, plotname, sep = "")
    jpeg(filename=plotname1)
    ymaxlim = max(p_t_set$values)
    yminlim = min(p_t_set$values)
    plot(NA, xlim=c(0,1), ylim = c(yminlim - 1, ymaxlim + 1), xlab="x", ylab="y")
    segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set[['endpoint']], p_t_set[['values']])
    dev.off()

    #The plotting everything on one level
    plotname2 = paste(dir, "level_", plotname, sep = "")
    jpeg(filename=plotname2)
    plot(NA, xlim=c(0,1), ylim = c(yminlim - 1, ymaxlim + 1), xlab="x", ylab="y")
    segments(p_t_set[['startpoint']], (yminlim + ymaxlim)/2, p_t_set[['endpoint']], (yminlim + ymaxlim)/2)
    dev.off()
     
    #Plotting the set A_plus
    plotname3 = paste(dir, "plus_", plotname, sep = "")
    jpeg(filename = plotname3)
    plot(NA, xlim=c(0,1), ylim = c(yminlim - 1, ymaxlim + 1), xlab="x", ylab="y")
    segments(p_t_set_plus[['startpoint']], p_t_set_plus[['values']], p_t_set_plus[['endpoint']], p_t_set_plus[['values']])
    dev.off()  
    
    #Plotting the set A_plus on one level
    plotname4 = paste(dir, "level_plus_", plotname, sep = "")
    plot(NA, xlim=c(0,1), ylim = c(yminlim, ymaxlim + 1), xlab="x", ylab="y")
    segments(p_t_set_plus[['startpoint']], (yminlim + ymaxlim)/2, p_t_set_plus[['endpoint']], (yminlim + ymaxlim)/2)
    dev.off()
    
    if (nrow(p_t_set_minus) > 0)
    #Plotting the set A_minus
    #plotname5 = paste(dir, "minus_", plotname, sep = "")
    #jpeg(filename=plotname5)
    #plot(NA, xlim=c(0,1), ylim = c(yminlim, ymaxlim + 1), xlab="x", ylab="y")
    #segments(p_t_set_minus[['startpoint']], p_t_set_minus[['values']], p_t_set_minus[['endpoint']], p_t_set_minus[['values']])
    #dev.off()  
    
    #Plotting the set A_minus
    #plotname6 = paste(dir, "level_minus_", plotname, sep = "")
    #jpeg(filename=plotname6)
    #plot(NA, xlim=c(0,1), ylim = c(yminlim, ymaxlim + 1), xlab="x", ylab="y")
    #segments(p_t_set_minus[['startpoint']], (yminlim + ymaxlim)/2, p_t_set_minus[['endpoint']], (yminlim + ymaxlim)/2)
    #dev.off()  
    
    return(list(p_t_set, p_t_set_plus, p_t_set_minus))
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    return(NULL)
  }
}

#a_t_set_1 <- plotting_all_rejected_intervals(y_data_1, g_t_set, kernel_ind, sigmahat, "Output/", "rightplot.jpg")
#a_t_set_2 <- plotting_all_rejected_intervals(y_data_2, g_t_set, kernel_ind, sigmahat, "Output/", "constantplot.jpg") #We expect to reject H_0 and the plot is everywhere
#a_t_set_3 <- plotting_all_rejected_intervals(y_data_3, g_t_set, kernel_ind, sigmahat, "Output/", "nullplot.jpg") #We expect to fail to reject H_0


###################################
#Loading the real data for England#
###################################

temperature <- read.table("cetml1659on.dat", header = TRUE, skip = 6)
yearly_tempr <- temperature[temperature$YEAR > -99, 'YEAR']

T_tempr <- length(yearly_tempr)

#Tuning parameters
L1 <-  floor(sqrt(T_tempr))
L2 <- floor(2 * sqrt(T_tempr))

sigmahat_tempr2 <- estimating_sigma_for_AR1(yearly_tempr, L1, L2)

g_t_set_tempr <- creating_g_set(T_tempr)
gaussian_quantile <- calculating_gaussian_quantile(T_tempr, g_t_set_tempr, kernel_ind, sqrt(sigmahat_tempr2), alpha)
results <- plotting_all_rejected_intervals(yearly_tempr, g_t_set_tempr, kernel_ind, sqrt(sigmahat_tempr2), "Output/", "tempr.jpg")
p_t_set <- results[[1]]
p_t_set_plus <- results[[2]]
p_t_set_minus <- results[[3]]