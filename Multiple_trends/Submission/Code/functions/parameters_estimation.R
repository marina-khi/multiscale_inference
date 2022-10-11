#############################
#Estimating the coefficients#
#############################

parameters_hp <- function(X_mat, n_ts, t_len, countries){
  hp_log           <- matrix(NA, ncol = n_ts, nrow = t_len)
  colnames(hp_log) <- countries
  
  hp_log_augm           <- matrix(NA, ncol = n_ts, nrow = t_len)
  colnames(hp_log_augm) <- countries
  
  beta      <- matrix(data = NA, nrow = 4, ncol = n_ts)
  alpha_vec <- c()

  i <- 1
  
  for (country in countries){
    tmp <- X_mat[X_mat$iso == country, ]
    tmp <- tmp[order(tmp$year), ]
    
    #Calculating first difference of the growth rate
    tmp <-
      tmp %>%
      mutate(delta_log_hp = coalesce(delta_log_hp, 0)) %>%
      mutate(delta_log_gdp = coalesce(delta_log_gdp, 0)) %>%
      mutate(delta_ltrate = coalesce(delta_ltrate, 0)) %>%
      mutate(delta_log_pop = coalesce(delta_log_pop, 0)) %>%
      mutate(delta_infl = coalesce(delta_infl, 0))
    
    #Estimating beta_i
    y_vec_tmp <- as.matrix(tmp[-1, 'delta_log_hp'])
    x_mat_tmp <- as.matrix(tmp[-1, c('delta_log_gdp', 'delta_log_pop',
                                     'delta_ltrate', 'delta_infl')])
    
    beta_tmp  <- solve(t(x_mat_tmp) %*% x_mat_tmp) %*% t(x_mat_tmp) %*% y_vec_tmp
    beta[, i] <- beta_tmp
    
    #Estimating alpha_i
    alpha_tmp     <- mean(tmp$log_hp - as.vector(as.matrix(tmp[, c('log_gdp',
                                                                   'log_pop',
                                                                   'ltrate',
                                                                   'infl')]) %*% beta_tmp))
    alpha_vec[i]  <- alpha_tmp
    
    
    #Working with adjusted time series and storing the original one
    y_vec_adj        <- tmp$log_hp - as.vector(as.matrix(tmp[, c('log_gdp',
                                                                 'log_pop',
                                                                 'ltrate',
                                                                 'infl')]) %*% beta_tmp) - alpha_tmp
    hp_log[, i]      <- tmp$log_hp
    hp_log_augm[, i] <- y_vec_adj
    i = i + 1
  }
  return(list("hp_log" = hp_log, "hp_log_augm" = hp_log_augm,
              "beta" = beta, "alpha" = alpha_vec))
}