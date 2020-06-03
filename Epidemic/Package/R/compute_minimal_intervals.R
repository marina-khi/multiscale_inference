#' Compute the set of minimal intervals as described in Duembgen (2002)
#'
#' @description The result of our multiscale test is the set of all intervals
#'              that have a corresponding test statistic bigger than
#'              the respective critical value. In order to make understandable
#'              statistic statements about these intervals, we need to find
#'              so-called minimal intervals: for a given set of intervals K,
#'              all intervals J such that K does not contain a proper subset
#'              of J are called minimal. Given K, this function computes the
#'              set of minimal intervals.
#'              Procedure is described in Duembgen (2002).
#' @export
#'
#'
#' @param dataset    Set of all intervals that have a corresponding test
#'                   statistic bigger than the respective critical value.
#' @return p_t_set   Set of minimal intervals

compute_minimal_intervals <- function(dataset) {
  set_cardinality <- nrow(dataset)
  if (set_cardinality > 1) {
    #Ordering the dataset such that we don't need to check previous intervals
    dataset <- dataset[order(dataset$startpoint, -dataset$endpoint), ]
    #We restore the indices after ordering
    rownames(dataset) <- seq_len(nrow(dataset))
    dataset$contains  <- numeric(set_cardinality)
    for (i in seq_len(set_cardinality - 1)) {
      for (j in (i + 1):set_cardinality) {
        if ((dataset$startpoint[i] <= dataset$startpoint[j]) &
            (dataset$endpoint[i] >= dataset$endpoint[j])) {
          dataset[["contains"]][i] <- 1
          #We mark all the intervals that contain at least one another interval
          break
        }
      }
    }
    #No we subset everything that is not marked
    p_t_set <- subset(dataset, contains == 0, select = -contains)
  } else {
    p_t_set <- dataset
  }
  return(p_t_set)
}
