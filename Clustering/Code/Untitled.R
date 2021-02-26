nadaraya_watson_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    result = sum((abs((grid_p - u) / bw) <= 1) * data_p)
    norm = sum((abs((grid_p - u) / bw) <= 1))
    return(result/norm)
  }
}