nadaraya_watson_smoothing <- function(u, data_p, grid_p, bw){
  if (length(data_p) != length(grid_p)){
    cat("Dimensions of the grid and the data do not match, please check the arguments")
    return(NULL)
  } else {
    result      = 0
    norm        = 0
    T_size      = length(data_p)
    for (i in 1:T_size){
      if (abs((grid_p[i] - u) / bw) <= 1){
        result = result + data_p[i]
        norm = norm + 1
      }
    }
    return(result/norm)
  }
}