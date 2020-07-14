#' Computes quantiles of the gaussian multiscale statistics.
#'
#' @description        Quantiles from the gaussian version of the test
#'                     statistics which are used to approximate
#'                     the critical values for the multiscale test.
#' @param t_len        An integer. Sample size.
#' @param grid         Grid of location-bandwidth points as produced by
#'                     the function \code{\link{construct_grid}} or
#'                     \code{\link{construct_weekly_grid}}, list with
#'                     the elements 'gset', 'bws', 'gtype'. If not provided,
#'                     then the defalt grid is produced and used.
#'                     For the construction of the default grid,
#'                     see \code{\link{construct_grid}}.
#' @param n_ts         An integer. Number of time series analyzed. Default is 1.
#' @param ijset        A matrix of integers. In case of multiple time series,
#'                     we need to know which pairwise comparisons to perform.
#'                     This matrix consists of all pairs of indices \eqn{(i, j)}
#'                     that we want to compare. If not provided, then all
#'                     possible pairwise comparison are performed.
#' @param sigma        A double that is equal to \eqn{\sqrt{long-run varaince}}
#'                     in case of n_ts = 1, or the overdispersion in case of
#'                     n_ts > 1.If not given, then the default is 1.
#' @param deriv_order  An integer. Order of the derivative of the trend
#'                     that is being investigated. Default is 0.
#' @param sim_runs     Number of simulation runs to produce quantiles.
#'                     Default is 1000.
#' @param probs        A numeric vector of probability levels (1-alpha)
#'                     for which the quantiles are computed.
#'                     Default is probs=seq(0.5,0.995,by=0.005).
#'
#' @return quant       Matrix with 2 rows where the first row contains
#'                     the vector of probabilities and the second contains
#'                     corresponding quantile of the gaussian statistics
#'                     distribution.
#' @export
#'
#' @examples
#' compute_quantiles(100)

compute_quantiles <- function(t_len, grid = NULL, n_ts = 1,
                              ijset = NULL, sigma = 1,
                              deriv_order = 0, sim_runs=1000,
                              probs=seq(0.5, 0.995, by = 0.005)) {

  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }

  if (is.null(ijset)) {
    ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
    ijset <- ijset[ijset$i < ijset$j, ]
  }

  ijset_cpp               <- as.matrix(ijset)
  ijset_cpp               <- as.vector(ijset_cpp)
  storage.mode(ijset_cpp) <- "integer"

  gset                   <- grid$gset
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp)
  storage.mode(gset_cpp) <- "double"
  storage.mode(sigma)    <- "double"

  phi <- simulate_gaussian(t_len = t_len, n_ts = n_ts, sim_runs = sim_runs,
                           gset = gset_cpp, ijset = ijset_cpp,
                           sigma = sigma, deriv_order = deriv_order)
  quant  <- as.vector(quantile(phi, probs = probs))
  quant  <- rbind(probs, quant)

  colnames(quant) <- NULL
  rownames(quant) <- NULL

  return(list("quant" = quant, "phi" = phi))
}
