library(microbenchmark)
source("functions.R")
dyn.load("psihat_statistic.dll")
source("psihat_statistic.R")

##############################
#Defining necessary constants#
##############################

sigma <- 2 #square root of long run-variance
T<-100 #length 
alpha <-0.05
noise_to_signal <- 1
sigmahat <- sigma #Here will be the estimate of the square root of the long-run variance sigma^2
kernel_f = "epanechnikov" #Only "epanechnikov" and "biweight" kernel functions are currently supported
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

#If we have already calculated quantiles and stored them in a file 'distribution.RData'
#then no need to calculate them once more, we just load them from this file.
#Ohterwise simulate the statistic 1000 times in order to calculate the qunatiles
calculating_gaussian_quantile <- function(T, g_t_set, sigmahat, alpha = 0.05){
  filename = paste("distribution_with_T_equal_to_", T, ".RData", sep = "")
  if(!file.exists(filename)) {
    gaussian_statistic_distribution <- replicate(1000, {
      z = rnorm(T, 0, 1)
      z_temp = sigmahat * z
      psistar_statistic(z_temp, g_t_set, kernel_ind, sigmahat)
    })
    save(gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }
  #Calculate the quantiles for gaussian statistic defined in the previous step
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)
  return(gaussian_quantile)
}

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


##################################
#Calculating the statistic itself#
##################################

#Auxiliary lines for checking function time
#microbenchmark(psihat_statistic(y_data_1, g_t_set, kernel_ind, sigmahat))
#microbenchmark(psihat_statistic_old(y_data_1, g_t_set, epanechnikov_kernel, sigmahat))


#Comparing C function with R function 
result <-psihat_statistic(y_data_1, g_t_set, kernel_ind, sigmahat)
g_t_set_with_values <- result[[1]]
psihat_statistic_value <- result[[2]]

result_old <-psihat_statistic_old(y_data_1, g_t_set, epanechnikov_kernel, sigmahat)
g_t_set_with_values <- result_old[[1]]
psihat_statistic_value <- result_old[[2]]

#Function that takes the data as argument and plots the regions where
#the corrected test statistic is bigger than the critical value of the 
#gaussian version. This function needs to be called only in environment
#with g_t_set and sigmahat defined.
#Returns the set with regions and value of the statistic in this regions.
plotting_all_rejected_intervals <-function(data, plotname, second_plotname){
  
  #Calculating our statistic
  result <-psihat_statistic(data, g_t_set, kernel_ind, sigmahat)
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
      
    #The plotting itself
#    jpeg(filename=plotname)
#    ymaxlim = max(p_t_set$values)
#    yminlim = min(p_t_set$values)
#    plot(NA, xlim=c(0,1), ylim = c(yminlim, ymaxlim + 1), xlab="x", ylab="y")
#    segments(p_t_set[['startpoint']], p_t_set[['values']], p_t_set[['endpoint']], p_t_set[['values']])
#    dev.off()  
    
    #The plotting itself
#    second_plotname <- cat("level_", plotname, sep="")
#    jpeg(filename=second_plotname)
#    ymaxlim = max(p_t_set$values)
#    yminlim = min(p_t_set$values)
#    plot(NA, xlim=c(0,1), ylim = c(yminlim, ymaxlim + 1), xlab="x", ylab="y")
#    segments(p_t_set[['startpoint']], (yminlim + ymaxlim)/2, p_t_set[['endpoint']], (yminlim + ymaxlim)/2)
#    dev.off()  

    return(list(a_t_set, p_t_set))
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    return(NULL)
  }
}

#a_t_set_1 <- plotting_all_rejected_intervals(y_data_1, "rightplot.jpg", "level_rightplot.jpg")[[1]]
#p_t_set_1 <- plotting_all_rejected_intervals(y_data_1, "rightplot.jpg", "level_rightplot.jpg")[[2]]
#a_t_set_2 <- plotting_all_rejected_intervals(y_data_2, "constantplot.jpg", "level_constantplot.jpg")[[1]] #We expect to reject H_0 and the plot is everywhere
#a_t_set_3 <- plotting_all_rejected_intervals(y_data_3, "nullplot.jpg", "level_nullplot.jpg")[[1]] #We expect to fail to reject H_0

#a_t_set_1[which.min(a_t_set_1$u),]