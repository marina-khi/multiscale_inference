SiZer_matrix_calculations <- function(y_data, SiZer_mat){
  #Wrapper of a C function from SiZer_matrix.C
  #Function that calculates the multiscale statistic \hat{\Psi}_T  under local linear kernel assumption.
  #Arguments:
  # y_data          vector of length T 
  # SiZer_mat       dataframe with columns "u", "h", "values_with_sign" (equal to 0)

  T        <- as.integer(length(y_data))
  N        <- as.integer(nrow(SiZer_mat))

  SiZer_mat_vec     <- unlist(SiZer_mat[c('u', 'h')])
  values_with_sign  <- vector(mode = "double",length = N)

  result <- .C("SiZer_matrix", y_data, T, SiZer_mat_vec, N, values_with_sign)
  
  SiZer_mat$values_with_sign <- result[[5]]

  return(SiZer_mat) 
}
