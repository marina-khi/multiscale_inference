dyn.load("C_code/psihat_statistic_ij.dll")
source("C_code/psihat_statistic_ij.R")
source("C_code/estimating_sigma.R")
dyn.load("C_code/estimating_sigma.dll")
source("functions.R")

##############################
#Defining necessary constants#
##############################

N     <- 3 #number of different time series
N_rep <- 1000

different_T     <- c(250, 350, 500, 1000) #Different lengths of time series for which we calculate size and power
different_alpha <- c(0.01, 0.05, 0.1) #Different alpha for which we calculate size and power

kernel_f <- "epanechnikov" #Only "epanechnikov" and "biweight" kernel functions are currently supported
a_hat <- 0.5

if (kernel_f == "epanechnikov"){
  kernel_ind = 1
} else if (kernel_f == "biweight"){
  kernel_ind = 2
} else {
  print('Currently only Epanechnikov and Biweight kernel functions are supported')
}


######################################
#Simulating the data for the iid case#
######################################

size_iid <- c()

for (T_size in different_T){
  simulated_data <- matrix(NA, nrow = T_size, ncol = N)
  colnames(simulated_data) <- 1:N
  L1 <- floor(sqrt(T_size))
  L2 <- floor(2 * sqrt(T_size))
  g_t_set_sim = creating_g_set(T_size)
  
  for (alpha in different_alpha){
    #Calculating gaussian quantiles for given T and alpha
    gaussian_quantile = calculating_gaussian_quantile_ij(T_size, N, g_t_set_sim, kernel_ind, alpha)
    
    size_of_the_test = replicate(N_rep, {
      sigmahat_vector_2 <- c()
      for (i in 1:N){
        simulated_data[, i] <- rnorm(T_size, 0, 1)
        sigma_i = estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
        sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_i * sigma_i)
      }
      sigmahat <- sqrt(sum(sigmahat_vector_2)/N)
      
      matrix_of_statistic <- matrix(, nrow = N, ncol = N)
      for (i in (1 : (N - 1))){
        for (j in ((i + 1):N)){
          result = psihat_statistic_ij(simulated_data[, i], simulated_data[, j], g_t_set_sim, kernel_ind, sigmahat)
          matrix_of_statistic[i, j] = result[[2]]
          #      cat(matrix_of_statistic[i, j], "with i =", i, "and j = ",j, "\n")
        }
      }
      result_with_trend = max(as.vector(matrix_of_statistic), na.rm=TRUE)
      if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
      d
    })
    size_iid = c(size_iid, sum(size_of_the_test)/N_rep)
    cat(sum(size_of_the_test)/N_rep, "with alpha = ", alpha, "and T = ", T_size, "\n")
  }
}


######################################
#Simulating the data for the AR1 case#
######################################

size_ar1 <- c()

for (T_size in different_T){
  simulated_data <- matrix(NA, nrow = T_size, ncol = N)
  colnames(simulated_data) <- 1:N
  L1 <- floor(sqrt(T_size))
  L2 <- floor(2 * sqrt(T_size))
  g_t_set_sim = creating_g_set(T_size)
  
  for (alpha in different_alpha){
    #Calculating gaussian quantiles for given T and alpha
    gaussian_quantile = calculating_gaussian_quantile_ij(T_size, N, g_t_set_sim, kernel_ind, alpha)
    
    size_of_the_test = replicate(N_rep, {
      sigmahat_vector_2 <- c()
      for (i in 1:N){
        simulated_data[, i] <- arima.sim(model = list(ar = a_hat), innov = rnorm(T_size, 0, 1), n = T_size)
        sigma_i = estimating_sigma_for_AR1(simulated_data[, i], L1, L2)[[1]]
        sigmahat_vector_2 <- c(sigmahat_vector_2, sigma_i * sigma_i)
      }
      sigmahat <- sqrt(sum(sigmahat_vector_2)/N)
      
      matrix_of_statistic <- matrix(, nrow = N, ncol = N)
      for (i in (1 : (N - 1))){
        for (j in ((i + 1):N)){
          result = psihat_statistic_ij(simulated_data[, i], simulated_data[, j], g_t_set_sim, kernel_ind, sigmahat)
          matrix_of_statistic[i, j] = result[[2]]
          #      cat(matrix_of_statistic[i, j], "with i =", i, "and j = ",j, "\n")
        }
      }
      result_with_trend = max(as.vector(matrix_of_statistic), na.rm=TRUE)
      if (result_with_trend > gaussian_quantile) {d = 1} else {d = 0}
      d
    })
    size_ar1 = c(size_ar1, sum(size_of_the_test)/N_rep)
    cat(sum(size_of_the_test)/N_rep, "with alpha = ", alpha, "and T = ", T_size, "\n")
  }
}

