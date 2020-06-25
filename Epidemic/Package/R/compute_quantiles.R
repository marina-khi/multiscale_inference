#' Computes quantiles of the gaussian multiscale statistics.
#'
#' @description        Quantiles from this distribution are used to approximate
#'                     the quantiles for the multiscale test.
#' @param t_len        An integer. Sample size.
#' @param grid         Grid of location-bandwidth points as produced by
#'                     the function \code{\link{construct_grid}}, list with
#'                     the elements 'gset', 'bws', 'gtype'. If not provided,
#'                     then the defalt grid is produced and used.
#' @param n_ts         An integer. Number of time series analyzed. Default is 1
#' @param sigma_vector A numeric vector of estimated sqrt(long-run variance)
#'                     for each time series. If not given, then the default is
#'                     a vector of ones of length n_ts (1, ..., 1).
#' @param deriv_order  An integer. Order of the derivative of the trend
#'                     that is being investigated. Default is 0.
#' @param sim_runs     Number of simulation runs to produce quantiles.
#'                     Default is 1000.
#' @param probs        A numeric vector of probability levels (1-alpha)
#'                     for which the quantiles are computed.
#'                     Default is probs=seq(0.5,0.995,by=0.005).
#' @param epidem       A logical parameter. True if we are investigating
#'                     epidemiological model. Default is FALSE.
#'
#' @return quant       Matrix with 2 rows where the first row contains
#'                     the vector of probabilities and the second contains
#'                     corresponding quantile of the gaussian statistics
#'                     distribution.
#' @export
#'
#' @examples
#' compute_quantiles(100)

compute_quantiles <- function(t_len, grid = NULL, n_ts = 1, sigma_vector = NULL,
                              deriv_order = 0, sim_runs=1000,
                              probs=seq(0.5, 0.995, by = 0.005),
                              epidem = FALSE) {

  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }

  gset                   <- grid$gset
  n                      <- as.integer(dim(grid$gset)[1])
  gset_cpp               <- as.matrix(gset)
  gset_cpp               <- as.vector(gset_cpp)
  storage.mode(gset_cpp) <- "double"
  storage.mode(epidem)   <- "logical"

  if (is.null(sigma_vector)) {
    sigma_vector <- rep(1, n_ts)
  }

  #filename = paste0("quantiles/distr_T_", t_len, "_n_", n_ts, ".RData")
  #if(!file.exists(filename)) {
  phi <- simulate_gaussian(t_len = t_len, n_ts = n_ts, sim_runs = sim_runs,
                           gset = gset_cpp, n = n, sigma_vec = sigma_vector,
                           deriv_order = deriv_order, epidem = epidem)
  quant  <- as.vector(quantile(phi, probs = probs))
  quant  <- rbind(probs, quant)

  colnames(quant) <- NULL
  rownames(quant) <- NULL

  #  save(quant, file = filename)
  #} else {
  #  load(filename)
  #}
  return(list("quant" = quant, "phi" = phi))
}
