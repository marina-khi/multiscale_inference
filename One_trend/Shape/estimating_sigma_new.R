estimating_variance_new <- function(y_data, L1, L2, order, K1 = order + 1, K2 = 10){
  
  a_hat_method1         <- AR_coefficients(y_data, L1, L2, rep(0,L2), order)
  sigma_eta_hat_method1 <- calculating_sigma_eta(y_data, a_hat_method1, order)
  sigma_hat_method1     <- sqrt(sigma_eta_hat_method1^2 / (1 - sum(a_hat_method1))^2) 
  
  corrections_value     <- corrections(a_hat_method1, sigma_eta_hat_method1, K2+1)
  a_hat_method2         <- AR_coefficients(y_data, K1, K2, corrections_value, order)
  sigma_eta_hat_method2 <- calculating_sigma_eta(y_data, a_hat_method2, order)
  sigma_hat_method2     <- sqrt(sigma_eta_hat_method2^2 / (1 - sum(a_hat_method2))^2) 
  
  return(list(sigma_hat_method2, a_hat_method2, sigma_eta_hat_method2))
}