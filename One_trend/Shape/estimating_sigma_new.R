estimating_sigma_new_method <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model    

  T_size      <- length(y_data)
  y_data_diff <- c(0, diff(y_data))

  Q_matrix <- matrix(data = NA, nrow = (T_size - p - 1), ncol = (p + 1))
  
  for (t in 1:(T_size - p - 1)){
    for (i in 1:(p+1)){
      Q_matrix[t, i] = Q(t + p + 1, i - 1, p, y_data_diff)  
    }
  }
  
  gamma_hat        <- c()
  Gamma_matrix_hat <- matrix(rep(0, p * p), nrow = p, ncol = p)
  
  for (i in 1:p){
    for (j in i:p){
      temp = sum(Q_matrix[, i + 1] * Q_matrix[, j + 1])
      Gamma_matrix_hat[i, j] = temp
      Gamma_matrix_hat[j, i] = temp
    }
    gamma_hat <- c(gamma_hat, sum(Q_matrix[, i + 1] * Q_matrix[, 1]))
  }
  
  # gamma_hat        <- c()
  # Gamma_matrix_hat <- matrix(rep(0, p * p), nrow = p, ncol = p)
  # 
  # for (i in 1:p){
  #   for (j in i:p){
  #     temp <- 0
  #     for (t in (p + 2):T_size){
  #       temp = temp + Q(t, i, p, y_data_diff) * Q(t, j, p, y_data_diff)
  #     }
  #     Gamma_matrix_hat[i, j] = temp
  #     Gamma_matrix_hat[j, i] = temp
  #   }
  # }
  # 
  # for (j in 1:p){
  #   temp <- 0
  #   for (t in (p + 2):T_size){
  #     temp = temp + Q(t, 0, p, y_data_diff) * Q(t, j, p, y_data_diff)
  #   }
  #   gamma_hat <- c(gamma_hat, temp)
  # }

  a_hat     <- solve(Gamma_matrix_hat) %*% gamma_hat
  
  sigma_eta_squared <- 0
  for (t in 1:(T_size - p - 1)){
    sigma_eta_squared <- sigma_eta_squared + (1 / (T_size - p - 1)) * (Q_matrix[t, 1] - sum(a_hat * Q_matrix[t, 2:(p + 1)]))^2
  }
  
  sigma_eta <- sqrt(sigma_eta_squared)
  sigma_hat <- sigma_eta / (abs(1 - sum(a_hat)))
  return(list(sigma_hat, a_hat, sigma_eta))
}