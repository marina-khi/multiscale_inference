# Epanechnikov kernel function, which is defined as f(x) = 3/4(1-x^2)
# for |x|<=1 and 0 elsewhere
epanechnikov_kernel <- function(x)
{
  if (abs(x)<=1)
  {
    result = 3/4 * (1 - x*x)
  } else {
    result = 0
  }
  return(result)
}


#Nadaraya-Watson estimator using the epanechnikov kernel. 
epanechnikov_smoothing <-function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result = 0
    norm = 0
    for (i in 1:length(data_p)){
      result = result + epanechnikov_kernel((u - grid_p[i])/bw) * data_p[i]
      norm = norm + epanechnikov_kernel((u - grid_p[i])/bw) 
    }
    return(result/norm)
  }
}


#Additive correction tern \lambda(h) that depends only on the bandwidth h
lambda <- function(h)
{
  result = tryCatch(sqrt(2*log(1/(2*h))), warning = function(w) print("h is exceeding h_max"))
  return(result)
}


#Creating g_t_set over which we are taking the maximum (from Section 2.1)
creating_g_set <- function(T, kernel_method){
  u <- seq(from = 5/T, to = 1, by = 5/T)
  h <- seq(from = 3/T, to = 1/4+3/T, by = 5/T)
  
  g_t_set_temp                  <- expand.grid(u = u, h = h) #Creating a dataframe with all possible combination of u and h
  g_t_set_temp$values           <- numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero
  g_t_set_temp$values_with_sign <- numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero
  
  if (kernel_method == "nw"){
    g_t_set <- subset(g_t_set_temp, u - h >= 0 & u + h <= 1, select = c(u, h, values, values_with_sign)) #Subsetting u and h such that [u-h, u+h] lies in [0,1]
  } else if (kernel_method == "ll"){
    g_t_set <- subset(g_t_set_temp, u >= 0 & u <= 1, select = c(u, h, values, values_with_sign)) #Subsetting u and h such that [u-h, u+h] lies in [0,1]
  } else {
    print('Given method is currently not supported')
  }
  
  g_t_set$lambda <- lambda(g_t_set[['h']]) #Calculating the lambda(h) in order to speed up the function psistar_statistic
  return(g_t_set)
}

#If we have already calculated quantiles and stored them in a file 'distr_T_....RData'
#then no need to calculate them once more, we just load them from this file.
#Ohterwise simulate the \Psi^star statistic 1000 times in order to calculate the quantiles
calculating_gaussian_quantile <- function(T, N_ts, g_t_set, kernel_method, alpha = 0.05){
  filename = paste("distribution/distr_T_", T, "_and_N_", N_ts, "_and_method_", kernel_method, ".RData", sep = "")
  if(!file.exists(filename)) {
      gaussian_statistic_distribution <- replicate(1000, {
        z_matrix <- matrix(rnorm(T * N_ts, 0, 1), T, N_ts)
        psistar_statistic(z_matrix, N_ts, g_t_set, sqrt(2), kernel_method)
        #statistic_vector <- c()
        # for (i in (1:(N_ts - 1))){
        #   for (j in ((i + 1):N_ts)){
        #     statistic_vector = c(statistic_vector, )
        #   }
        # }
        # max(statistic_vector)
      })
    save(gaussian_statistic_distribution, file = filename)
  } else {
    load(filename)
  }
  #Calculate the quantiles for gaussian statistic defined in the previous step
  gaussian_quantile <- quantile(gaussian_statistic_distribution, probs = (1 - alpha), type = 1)
  cat("Gaussian quantile is", gaussian_quantile, "\n")
  return(gaussian_quantile)
}


#This functions finds minimal intervals as described in Duembgen(2002)
choosing_minimal_intervals <- function(dataset){
  set_cardinality <- nrow(dataset) 
  if (set_cardinality > 1) {
    dataset <- dataset[order(dataset$startpoint, -dataset$endpoint),] #Ordering such that we do not have 
    rownames(dataset) <- 1:nrow(dataset) #restoring the indices after ordering
    dataset[['contains']] <- numeric(set_cardinality)
    for (i in 1:(set_cardinality-1)){
      for (j in (i+1):set_cardinality){
        if ((dataset$startpoint[i] <= dataset$startpoint[j]) & (dataset$endpoint[i] >= dataset$endpoint[j])) {
          dataset[['contains']][i] <- 1 #We are marking all the intervals that contain at least one another interval
          break
        }
      }
    }
    p_t_set <- subset(dataset, contains == 0, select = c(startpoint, endpoint, values))#Subsetting everything not marked
  }
}






#Function that plot the intervals (either minimal or not) to a file with name = "name"
#with respect to the value of the statistic on this interval. The set with the intervals should
#have the column named named "startpoint", "endpoint", and "values"
plotting <- function(set, name){
  jpeg(filename=name)
  ymaxlim = max(set$values)
  yminlim = min(set$values)
  plot(NA, xlim=c(0,1), ylim = c(yminlim - 1, ymaxlim + 1), xlab="x", ylab="y")
  segments(set[['startpoint']], set[['values']], set[['endpoint']], set[['values']])
  dev.off()
}


#Function that plot all the intervals (either minimal or not) to a file with name = "name"
#on one level. The set with the intervals should have the columns named "startpoint", "endpoint", and "values"
plotting_one_level <- function(set, name){
  jpeg(filename=name)
  ymaxlim = max(set$values)
  yminlim = min(set$values)
  plot(NA, xlim=c(0,1), ylim = c(yminlim - 1, ymaxlim + 1), xlab="x", ylab="y")
  segments(set[['startpoint']], (yminlim + ymaxlim)/2, set[['endpoint']], (yminlim + ymaxlim)/2)
  dev.off()
}


#Function that takes the data as argument and plots the regions where
#the corrected test statistic is bigger than the critical value of the 
#gaussian version.
#Returns the set with regions and value of the statistic in this regions.
plotting_all_rejected_intervals <-function(data, g_t_set, quantile, kernel_ind, sigma, dir, plotname){
  
  #Calculating our statistic
  result                 <- psihat_statistic(data, g_t_set, kernel_ind, sigma)
  g_t_set_with_values    <- result[[1]]
  psihat_statistic_value <- result[[2]]
  
  #And now the testing itself
  if (psihat_statistic_value > quantile) {
    cat("We reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    
    #The collection of intervals where the corrected test statistic lies above the critical value (as in (2.6))
    a_t_set <- subset(g_t_set_with_values, values > gaussian_quantile, select = c(u, h, values))
    p_t_set <- data.frame('startpoint' = a_t_set$u - a_t_set$h,
                          'endpoint' = a_t_set$u + a_t_set$h,
                          'values' = a_t_set$values)
    p_t_set <- choosing_minimal_intervals(p_t_set)
    
    #The collection of intervals where the values_with_sign > gaussian_quantile + lambda (as in (2.6))
    a_t_set_plus <- subset(g_t_set_with_values, values_with_sign > gaussian_quantile + lambda, select = c(u, h, values_with_sign))
    p_t_set_plus <- data.frame('startpoint' = a_t_set_plus$u - a_t_set_plus$h,
                               'endpoint' = a_t_set_plus$u + a_t_set_plus$h,
                               'values' = a_t_set_plus$values_with_sign)
    p_t_set_plus <- choosing_minimal_intervals(p_t_set_plus)
    
    #The collection of intervals where the -values_with_sign > gaussian_quantile + lambda (as in (2.6))
    a_t_set_minus <- subset(g_t_set_with_values, -values_with_sign > gaussian_quantile + lambda, select = c(u, h, values_with_sign))
    p_t_set_minus <- data.frame('startpoint' = a_t_set_minus$u - a_t_set_minus$h,
                                'endpoint' = a_t_set_minus$u + a_t_set_minus$h,
                                'values' = a_t_set_minus$values_with_sign)
    p_t_set_minus <- choosing_minimal_intervals(p_t_set_minus)
    
    #The plotting itself
    plotname1 = paste(dir, plotname, sep = "")
    plotname2 = paste(dir, "level_", plotname, sep = "")
    
    plotting(p_t_set, plotname1)
    plotting_one_level(p_t_set, plotname2)
    

    return(list(p_t_set, p_t_set_plus, p_t_set_minus))
  } else {
    cat("We fail to reject H_0 with probability", alpha, "Psihat_statistic = ", psihat_statistic_value,
        "Gaussian quantile value = ", gaussian_quantile, "\n")
    return(NULL)
  }
}

creating_matrix_and_texing <- function(vect, vect_T, vect_alpha, filename){
  matrix_ <- matrix(vect, nrow = length(vect_T), ncol = length(vect_alpha), byrow = TRUE)
  rownames(matrix_) <- vect_T
  colnames(matrix_) <- vect_alpha
  print.xtable(xtable(matrix_), type="latex", file=filename)
  }