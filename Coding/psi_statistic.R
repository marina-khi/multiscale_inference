# Epanechnikov kernel function,
# which is defined as f(x) = 3/4(1-x^2)
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


# Additive correction tern \lambda(h) that depends only on the bandwidth h
lambda <- function(h)
{
  result = tryCatch(sqrt(2*log(1/(2*h))), warning = function(w) print("h is exceeding h_max"))
  return(result)
}


psihat_statistic <- function(y_data, g_t_set, kernel_ind = 1, sigmahat){
  #Wrapper of a C function from psihat_statistic.C
  #Arguments:
  # y_data          vector of length T 
  # g_t_set         dataframe with columns "u", "h", "values" (equal to 0), "lambda"
  # kernel_function type of kernel function used. Currently only epanechnikov and biweight kernels are supported
  # sigmahat        appropriate estimate of sqrt(long-run variance)
  
  T <- as.integer(length(y_data))
  N <- as.integer(nrow(g_t_set))

  #Passing an indicator to C functions that determines a necessary kernel function
  #if (kernel_function == "epanechnikov") {kernel_ind = 1}
  #else {kernel_ind = 2}

  storage.mode(kernel_ind)  <- "integer"
  storage.mode(sigmahat)    <- "double"
  
  g_t_set_vec <- unlist(g_t_set[c('u', 'h', 'lambda')])
  values      <- vector(mode = "double",length = N)
  maximum     <- vector(mode = "double", length = 1)
  
  result <- .C("psihat_statistic", y_data, T, g_t_set_vec, N, kernel_ind, sigmahat, maximum, values)

#  cat("kernel = ", result[[5]], ", T=", result[[2]], ", N=", result[[4]], ", sigmahat = ", result[[6]],
#      "maximum = ", result[[7]], "len of values = ", length(result[[8]]), sep=" ")
  
  g_t_set_result  <- c(g_t_set$u, g_t_set$h, result[[8]], g_t_set$lambda)
  statistic       <- result[[7]]

  return(list(g_t_set_result, statistic)) 
}

# Function that calculates the auxiliary statistic \Psi^star_T.
# The only difference with the previus function is in the return values.
psistar_statistic <- function(y_data, g_t_set, kernel_function = epanechnikov_kernel, sigmahat) {
  g_t_set$values <- abs(mapply(psi_average, u = g_t_set$u, h = g_t_set$h,
                               MoreArgs = list(data = y_data, k_function = kernel_function))/sigmahat) - .subset2(g_t_set, 'lambda')
  result = max(g_t_set$values)
  return(result)
}

#This functions chooses minimal intervals as described in Dumbgen()
choosing_minimal_intervals <- function(dataset){
  set_cardinality <- nrow(dataset) 
  if (set_cardinality > 1) {
    dataset <- dataset[order(dataset$startpoint, -dataset$endpoint),]
    rownames(dataset) <- 1:nrow(dataset)
    dataset[['contains']] <- numeric(set_cardinality)
    for (i in 1:(set_cardinality-1)){
      for (j in (i+1):set_cardinality){
        if ((dataset$startpoint[i] <= dataset$startpoint[j])&(dataset$endpoint[i] >= dataset$endpoint[j])) {
          #          cat("[",dataset$startpoint[i], ", ", dataset$endpoint[i],"], [", dataset$startpoint[j], ", ", dataset$endpoint[j], "]", sep="")
          dataset[['contains']][i] <- 1
          break
        }
      }
    }
    p_t_set <- subset(dataset, contains == 0, select = c(startpoint, endpoint, values))
  }
}


#Creating g_t_set over which we are taking the maximum (from Section 2.1)
creating_g_set <- function(T){
  u <- seq(4/T, 1, length.out = T/4)
  h <- seq(3/T, 1/4+3/T, length.out = T/20)

  #u <- seq(1/T, 1, length.out = T)#This is useful for testing small tme series, for example T=3 or 10
  #h <- seq(1/T, 1, length.out = T)
  
  g_t_set_temp <- expand.grid(u = u, h = h) #Creating a dataframe with all possible combination of u and h
  g_t_set_temp$values <-numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero
  
  g_t_set <- subset(g_t_set_temp, u - h >= 0 & u + h <= 1, select = c(u, h, values)) #Subsetting u and h such that [u-h, u+h] lies in [0,1]
  g_t_set$lambda <- lambda(g_t_set[['h']]) #Calculating the lambda(h) beforehand in order to speed up the function psistar_statistic
  return(g_t_set)
}