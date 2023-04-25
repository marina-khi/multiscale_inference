#####################
#Estimating variance#
#####################

sigma <- function(matrix_, n_ts, q_ = 15, r_ = 10){
  #Calculating each sigma_hat_i separately
  lrv_vector            <- c()
  ahat_vector           <- c()
  innovation_var_vector <- c()
  
  for (i in 1:n_ts){
    AR.struc              <- estimate_lrv(data = matrix_[, i], q = q_,
                                          r_bar = r_, p = 1)
    lrv_vector            <- c(lrv_vector, AR.struc$lrv)
    ahat_vector           <- c(ahat_vector, AR.struc$ahat)
    innovation_var_vector <- c(innovation_var_vector, AR.struc$vareta)
  }
  return(lrv_vec = lrv_vector, a_vec = ahat_vector,
         innovation_var_vec = innovation_var_vector)
}