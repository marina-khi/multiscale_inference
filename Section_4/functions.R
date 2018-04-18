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

#Functions needed for local linear smoothing
s_t_1 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

s_t_2 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h) * ((i/T_size - u) / h) * ((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}

s_t_0 <- function(u, h, T_size) {
  result = 0
  for (i in 1:T_size) {
    result = result + epanechnikov_kernel((i/T_size - u) / h)
  }
  return(result / (T_size * h));
}


#Local Linear estimator using the epanechnikov kernel. 
local_linear_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    s_t_2_value = s_t_2(u, bw, T_size)
    s_t_1_value = s_t_1(u, bw, T_size) 
    for (i in 1:T_size){
      k = (s_t_2_value - s_t_1_value * ((grid_p[i] - u) / bw)) * epanechnikov_kernel((grid_p[i] - u) / bw)
      result = result + k * data_p[i]
      norm = norm + k 
    }
    return(result/norm)
  }
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
  sigma_vector <- rep(1, N_ts)
  if(!file.exists(filename)) {
      gaussian_statistic_distribution <- replicate(1000, {
        z_matrix <- matrix(rnorm(T * N_ts, 0, 1), T, N_ts)
        psistar_statistic(z_matrix, N_ts, g_t_set, sigma_vector, kernel_method)
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


#Create a matrix (for size and power table for example) and write them in the tex file
creating_matrix_and_texing <- function(vect, vect_T, vect_alpha, filename){
  matrix_ <- matrix(vect, nrow = length(vect_T), ncol = length(vect_alpha), byrow = TRUE)
  rownames(matrix_) <- vect_T
  colnames(matrix_) <- vect_alpha
  
  addtorow <- list()
  addtorow$pos <- list(0, 0)
  addtorow$command <- c("& \\multicolumn{3}{c}{nominal size $\\alpha$} \\\\\n",
                        "$T$ & 0.01 & 0.05 & 0.1 \\\\\n") 
  print.xtable(xtable(matrix_, digits = c(3), align = "cccc"), type="latex",  file=filename, add.to.row = addtorow, include.colnames = FALSE)
}