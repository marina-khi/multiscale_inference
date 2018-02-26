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


# Kernel weight function that calculates
# the vector of (\omega_{1,T}(u,h), \ldots, \omega_{T,T}(u,h)).
# Meanwhile, the function also calculates
# the norm ||K||_{u,h,T} of the kernel function.
omega <- function(T, u, h, k_function)
{
  result_temp = c()
  K_norm_temp = 0 
  for (i in 1:T) {
    x = (u - i/T)/h
    k = k_function(x)
    result_temp[i] = k
    K_norm_temp = K_norm_temp + k^2
  }
  K_norm = sqrt(K_norm_temp)
  result = result_temp / K_norm
  #  cat("K_norm:", K_norm, " omega_t:", result)
  return(result)
}


# Kernel average function \psi_T(u,h) that takes u, h, data
# and the type of kernel function as arguments.
# The data can be y_data for \hat{\Psi}_T or independent gaussian rv
# z_temp = sigma*z for \Psi^star_T.
# The output is one value for each u and h.
psi_average <- function(data, u, h, k_function)
{
  T = length(data)
  result = sum(omega(T, u, h, k_function)*data)
  return(result)
}


# Function that calculates the multiscale statistic \hat{\Psi}_T.
# It takes the following entities as arguments.
#   y-data: the data
#   g_t_set: range of different locations u and bandwidths h
#   k_function: type of kernel function
#   sigmahat: the estimator of the square root of the long-run error variance \sigma^2
# It produces the value of the test statistic as an output
psihat_statistic <- function(y_data, g_t_set, kernel_function = epanechnikov_kernel, sigmahat) {
  g_t_set$values <- abs(mapply(psi_average, u = g_t_set$u, h = g_t_set$h,
                               MoreArgs = list(data = y_data, k_function = kernel_function))/sigmahat) - .subset2(g_t_set, 'lambda')
  result = max(g_t_set$values)
  return(list(g_t_set, result))
}

#old version
#psihat_statistic <- function(y_data, g_t_set, k_function = epanechnikov_kernel, sigmahat) {
#  g_t_set_card = nrow(g_t_set)
#  for (i in 1:g_t_set_card) {
#    g_t_set[['values']][i] <- abs(psi_average(y_data, .subset2(g_t_set,'u')[i], .subset2(g_t_set,'h')[i], k_function)/sigmahat) - .subset2(g_t_set, 'lambda')[i]
#  }
#  result = max(g_t_set$values)
#  return(result)
#}


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
  
  g_t_set_temp <- expand.grid(u = u, h = h) #Creating a dataframe with all possible combination of u and h
  g_t_set_temp$values <-numeric(nrow(g_t_set_temp)) # Setting the values of the statistic to be zero
  
  g_t_set <- subset(g_t_set_temp, u - h >= 0 & u + h <= 1, select = c(u, h, values)) #Subsetting u and h such that [u-h, u+h] lies in [0,1]
  g_t_set$lambda <- lambda(g_t_set[['h']]) #Calculating the lambda(h) in order to speed up the function psistar_statistic
  return(g_t_set)
}