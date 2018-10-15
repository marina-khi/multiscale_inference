ptm <- proc.time()
SiZer_matrix <- calculating_SiZer_matrix(different_i, different_h, T_size, T_star, alpha, gamma)  
proc.time() - ptm

ptm <- proc.time()
SiZer_matrix_1 <- calculating_SiZer_matrix_1(different_i, different_h, T_size, T_star, alpha, gamma)  
proc.time() - ptm

y_data_ar_1 <- arima.sim(model = list(ar = a_1), n = T_size, innov = rnorm(T_size, 0, sigma_eta))

for (i in 1:nrow(SiZer_matrix)){
  SiZer_matrix$oldmethod[[i]] <- (SiZer_matrix$XtWX_inv_XtW[[i]] %*% y_data_ar_1)[2]
}

SiZer_matrix_1 <- SiZer_matrix_calculations(y_data_ar_1, SiZer_matrix_1)

