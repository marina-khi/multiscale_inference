#This functions finds minimal intervals as described in Duembgen(2002)
choosing_minimal_intervals <- function(dataset){
  set_cardinality <- nrow(dataset) 
  if (set_cardinality > 1) {
    dataset <- dataset[order(dataset$startpoint, -dataset$endpoint),] #Ordering such that we don't need to check previous intervals
    rownames(dataset) <- 1:nrow(dataset) #restoring the indices after ordering
    dataset[['contains']] <- numeric(set_cardinality)
    for (i in 1:(set_cardinality-1)){
      for (j in (i+1):set_cardinality){
        if ((dataset$startpoint[i] <= dataset$startpoint[j]) & (dataset$endpoint[i] >= dataset$endpoint[j])) {
          dataset[['contains']][i] <- 1 #We are marking all the intervals that contain at least one another interval
          break
        }
      }
    }
    p_t_set <- subset(dataset, contains == 0, select = c(startpoint, endpoint, values))#Subsetting everything not marked
  } else {p_t_set <- dataset}
  return(p_t_set) 
}
