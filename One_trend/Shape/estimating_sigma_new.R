estimating_sigma_new <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model    

  T_size      <- length(y_data)
  y_data_diff <- c(0, diff(y_data))

  # Q_matrix <- matrix(data = NA, nrow = (T_size - p - 1), ncol = (p + 1))
  # 
  # for (t in 1:(T_size - p - 1)){
  #   for (i in 1:(p+1)){
  #     Q_matrix[t, i] = Q(t + p + 1, i - 1, p, y_data_diff)  
  #   }
  # }
  # 
  # gamma_hat        <- c()
  # Gamma_matrix_hat <- matrix(rep(0, p * p), nrow = p, ncol = p)
  # 
  # for (i in 1:p){
  #   for (j in i:p){
  #     temp = sum(Q_matrix[, i + 1] * Q_matrix[, j + 1])
  #     Gamma_matrix_hat[i, j] = temp
  #     Gamma_matrix_hat[j, i] = temp
  #   }
  #   gamma_hat <- c(gamma_hat, sum(Q_matrix[, i + 1] * Q_matrix[, 1]))
  # }

  gamma_hat        <- c()
  Gamma_matrix_hat <- matrix(rep(0, p * p), nrow = p, ncol = p)

  for (i in 1:p){
    for (j in i:p){
      temp <- 0
      for (t in (p + 2):T_size){
        temp = temp + Q(t, i, p, y_data_diff) * Q(t, j, p, y_data_diff)
      }
      Gamma_matrix_hat[i, j] = temp
      Gamma_matrix_hat[j, i] = temp
    }
  }

  for (j in 1:p){
    temp <- 0
    for (t in (p + 2):T_size){
      temp = temp + Q(t, 0, p, y_data_diff) * Q(t, j, p, y_data_diff)
    }
    gamma_hat <- c(gamma_hat, temp)
  }

  a_hat     <- solve(Gamma_matrix_hat) %*% gamma_hat
  
  # sigma_eta_squared <- 0
  # for (t in 1:(T_size - p - 1)){
  #   sigma_eta_squared <- sigma_eta_squared + (1 / (T_size - p - 1)) * (Q_matrix[t, 1] - sum(a_hat * Q_matrix[t, 2:(p + 1)]))^2
  # }
  
  sigma_eta_squared <- 0
  for (t in (p+2):T_size){
    sigma_eta_squared <- sigma_eta_squared + (1 / (T_size - p - 1)) * (Q(t, 0, p, y_data_diff) - a_hat * Q(t, 1, p, y_data_diff))^2
  }
  
  sigma_eta <- sqrt(sigma_eta_squared)
  sigma_hat <- sigma_eta / (abs(1 - sum(a_hat)))
  return(list(sigma_hat, a_hat, sigma_eta))
}

estimating_sigma_old_simple <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model
  
  y_data_diff <- c(0, diff(y_data))
  
  gamma_hat_temp <- sapply(0:(p+1), FUN = sample_autocovariance, y_data_diff)
  if (p == 1){
    Gamma_matrix_hat <- gamma_hat_temp[2]
    gamma_hat        <- gamma_hat_temp[3]
    a_hat            <- gamma_hat / Gamma_matrix_hat     
  } else {
    Gamma_matrix_hat <- matrix(data = NA, nrow = p, ncol = p)
    for (i in 1:p){
      for (j in 1:p){
        if (i >= j) {
          Gamma_matrix_hat[i, j] = gamma_hat[i - j + 2]
        } else {
          Gamma_matrix_hat[i, j] = gamma_hat[j - i]
        }
      }
    }
    gamma_hat <-  gamma_hat_temp[3:(p+2)]
    a_hat     <- inv(Gamma_matrix_hat) %*% gamma_hat
  }
#  sigma_eta_squared <- sum(a_hat * gamma_hat_temp[2:(p+1)]) - gamma_hat_temp[1]
#  sigma_eta <- sqrt(sigma_eta_squared)
#  sigma_hat <- sigma_eta / (abs(1 - sum(a_hat)))
  return(a_hat)
}  



estimating_sigma_overidentified <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model
  
  y_data_diff <- c(0, diff(y_data))
  
  gamma_hat <- sapply(1:(p+2), FUN = sample_autocovariance, y_data_diff)

  Kappa_matrix_hat <- matrix(data = NA, nrow = p, ncol = p)
  Kappa_vector     <- c()
  for (i in 1:p){
    for (j in 1:p){
        Kappa_matrix_hat[i, j] = sum(gamma_hat[(2 - i):(p + 2 - i)] * gamma_hat[(2 - j):(p + 2 - j)])
    }
  Kappa_vector <- c(Kappa_vector, sum(gamma_hat[2:(p + 2)] * gamma_hat[(2 - i):(p + 2 - i)]))
  }
  a_hat     <- solve(Kappa_matrix_hat) %*% Kappa_vector
  return(a_hat)
}

estimating_sigma_overidentified_twice <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model
  
  y_data_diff <- c(0, diff(y_data))
  
  gamma_hat <- sapply(1:(p+3), FUN = sample_autocovariance, y_data_diff)
  
  Kappa_matrix_hat <- matrix(data = NA, nrow = p, ncol = p)
  Kappa_vector     <- c()
  for (i in 1:p){
    for (j in 1:p){
      Kappa_matrix_hat[i, j] = sum(gamma_hat[(2 - i):(p + 3 - i)] * gamma_hat[(2 - j):(p + 3 - j)])
    }
    Kappa_vector <- c(Kappa_vector, sum(gamma_hat[2:(p + 3)] * gamma_hat[(2 - i):(p + 3 - i)]))
  }
  a_hat     <- solve(Kappa_matrix_hat) %*% Kappa_vector
  return(a_hat)
}


estimating_sigma_overidentified_3 <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model
  
  y_data_diff <- c(0, diff(y_data))
  
  gamma_hat <- sapply(1:(p+4), FUN = sample_autocovariance, y_data_diff)
  
  Kappa_matrix_hat <- matrix(data = NA, nrow = p, ncol = p)
  Kappa_vector     <- c()
  for (i in 1:p){
    for (j in 1:p){
      Kappa_matrix_hat[i, j] = sum(gamma_hat[(2 - i):(p + 4 - i)] * gamma_hat[(2 - j):(p + 4 - j)])
    }
    Kappa_vector <- c(Kappa_vector, sum(gamma_hat[2:(p + 4)] * gamma_hat[(2 - i):(p + 4 - i)]))
  }
  a_hat     <- solve(Kappa_matrix_hat) %*% Kappa_vector
  return(a_hat)
}


estimating_sigma_oracle <- function(y_data, p){
  #Function that calculates the estimate of the square root of long-run variance sigma^2
  #Arguments:
  #y_data           vector of data length T 
  #p                order of the AR(p) model
  
  T_size <- length(y_data)
  X_matrix <- matrix(data = NA, nrow = T_size - p, ncol = p)
  y_vector <- y_data[(p+1):T_size]
  for (i in 1:(T_size - p)){
      X_matrix[i, ] <- rev(y_data[i:(i + p - 1)])
  }
  a_hat     <- solve(t(X_matrix) %*%X_matrix) %*% t(X_matrix) %*% y_vector
  return(a_hat)
}