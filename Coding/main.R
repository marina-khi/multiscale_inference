source("functions.R")
library(microbenchmark)
library(gsubfn)

#Defining constants 
sigma <- 2 #square root of long run-variance
T<-100
alpha <-0.05
noise_to_signal <- 2


#Creating g_t_set over which we are taking the maximum (from Section 2.1)
creating_g_set <- function(T){
  u <- seq(4/T, 1, length.out = T/4)
  h <- seq(3/T, 1/4+3/T, length.out = T/20)

  g_t_set_temp <- expand.grid(u = u, h = h) #Creating a dataframe with all possible combination of u and h
  g_t_set_temp$values <-numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero

  g_t_set <- subset(g_t_set_temp, u - h >= 0 & u + h <= 1, select = c(u, h, values)) #Subsetting u and h such that [u-h, u+h] lies in [0,1]
  g_t_set$lambda <- lambda(g_t_set[['h']]) #Calculating the lambda(h) in order to speed up the function psistar_statistic
  return(g_t_set)
}

#If we have already calculated quantiles and stored them in a file 'distribution.RData'
#then no need to calculate them once more, we just load them from this file.
calculating_gaussian_quantile <- function(T){
  filename = paste("distribution_with_T_equal_to_", T, ".RData", sep = "")
  if(!file.exists(filename)) {
    gaussian_statistic_distribution <- replicate(1000, {
      z = rnorm(T, 0, 1)
      z_temp = sigma * z
      psistar_statistic(z_temp, g_t_set, epanechnikov_kernel, sigma)
    })
    save(gaussian_statistic_distribution, file = filename)
    } else {
    load(filename)
    }
  #Calculate the quantiles for gaussian statistic defined in the previous step
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)
  return(gaussian_quantile)
}

############################################
#Calculating quantiles for T=100, 500, 1000#
############################################
for (t in c(100)) {
  g_t_set <- creating_g_set(t)
  gaussian_quantile <- calculating_gaussian_quantile(t)
}

##################################
#Calculating the statistic itself#
##################################


#Defining the data Y = m + noise
set.seed(1) #For reproducibility

a <- sqrt(4*48/(5*noise_to_signal))
b <-0.5 * T/a
m <- numeric(T)
for (i in 1:T) {
  if (i/T < 0.5) {m[i] <- 0} else {m[i] <- (i - 0.5*T)/b}
} 
y_data_1 <- m + rnorm(T, 0, sigma)#Adding to noise a function that is 0 on the first half and linear on the second half
y_data_2 <- 1 + rnorm(T, 0, sigma)#Adding constant function m=1
y_data_3 <- rnorm(T, 0, sigma)#Here the null hypothesis is true, m=0

sigmahat <- sigma #Here will be the estimate of the square root of the long-run variance sigma^2

g_t_set <- creating_g_set(T)

#Auxiliary line for checking function time
#microbenchmark(psihat_statistic(y_data, g_t_set, epanechnikov_kernel, sigmahat))


#Function that takes the data as argument and plots the regions where
#the corrected test statistic is bigger than the critical value of the 
#gaussian version. This function needs to be called only in environment
#with g_t_set and sigmahat defined.
#Returns the set with regions and value of the statistic in this regions.
plotting_all_rejected_intervals <-function(data, plotname){
  
  #Calculating our statistic
  result <-psihat_statistic(data, g_t_set, epanechnikov_kernel, sigmahat)
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
                          'value' = a_t_set$values)
    
#    p_t_set <- choosing_minimal_intervals(p_t_set)
      
    #The plotting itself
    jpeg(filename=plotname)
    ymaxlim = max(a_t_set$values)
    yminlim = min(a_t_set$values)
    plot(NA, xlim=c(0,1), ylim = c(yminlim, ymaxlim + 1), xlab="x", ylab="y")
    segments(a_t_set[['u']] - a_t_set[['h']], a_t_set[['values']], a_t_set[['u']] + a_t_set[['h']], a_t_set[['values']])
    dev.off()  

    return(list(a_t_set, p_t_set))
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    return(NULL)
  }
}

a_t_set_1 <- plotting_all_rejected_intervals(y_data_1, "rightplot.jpg")[[1]]
p_t_set_1 <- plotting_all_rejected_intervals(y_data_1, "rightplot.jpg")[[2]]#We expect to reject H_0 and the plot is mostly to the right
a_t_set_2 <- plotting_all_rejected_intervals(y_data_2, "constantplot.jpg")[[1]] #We expect to reject H_0 and the plot is everywhere
a_t_set_3 <- plotting_all_rejected_intervals(y_data_3, "nullplot.jpg")[[1]] #We expect to fail to reject H_0

a_t_set_1[which.min(a_t_set_1$u),]