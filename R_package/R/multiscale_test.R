#' Carries out the multiscale test given that the values the estimatates of
#' long-run variance have already been computed.
#'
#' @param n_ts           Number of time series analysed. Default is 1.
#' @param data           Vector (in case of n_ts = 1) or matrix (in case of
#'                       n_ts > 1) that contains (a number of) time series
#'                       that needs to be analyzed.
#' @param sigma          A double. Estimator of the square root of the long-run
#'                       variance in case of one time series, or estimator of
#'                       the overdispersion parameter in case of multiple time
#'                       series.
#' @param n_ts           Number of time series analysed. Default is 1.
#' @param grid           Grid of location-bandwidth points as produced by
#'                       the function \code{\link{construct_grid}}.
#' @param ijset          Matrix of all pairs of indices (i, j) that we want
#'                       to compare.
#' @param alpha          Significance level. Default is 0.05.
#' @param sim_runs       Number of simulation runs to produce quantiles.
#'                       Default is 1000.
#' @param deriv_order    An integer. Order of the derivative of the trend
#'                       that is being investigated. Default is 0.
#' @export
#'
#' @return quant         Quantile that wasused for testing calculated from
#'                       the gaussian distribution.
#' @return statistics    Value of the multiscale statistics.
#' @return test_matrix   Return in case of n_ts = 1. Matrix of test results
#'                       for the multiscale test defined in
#'                       Khismatullina and Vogt (2019).
#'                       test_matrix[i,j] = -1: test rejects the null for the
#'                                j-th location u and the i-th bandwidth h and
#'                                indicates a decrease in the trend
#'                       test_matrix[i,j] = 0:  test does not reject the null
#'                                for the j-th location u and the i-th
#'                                bandwidth h
#'                       test_matrix[i,j] = 1:  test rejects the null for the
#'                                j-th location u and the i-th bandwidth h and
#'                                indicates an increase in the trend
#'                       test_matrix[i,j] = 2:  no test is carried out at j-th
#'                                location u and i-th bandwidth h (because
#'                                the point (u, h) is excluded from the grid
#'                                as specified by the 'deletions' option
#'                                in the function \code{\link{construct_grid}})
#' @return test_matrices  Return in case of n_ts > 1. List of matrices,
#'                        each matrix contains test results for the pairwise
#'                        comparison between time series.
#'                        Each matrix is coded exactly as in case of n_ts = 1.
#' @return gset_with_vals Either a matrix (in case of n_ts = 1) or a list of
#'                        matrices (in case of n_ts > 1) that contains test
#'                        results together with location-bandwidth points.
#'
multiscale_test <- function(data, sigma, n_ts = 1, grid = NULL,
                            ijset = NULL, alpha = 0.05, sim_runs = 1000,
                            deriv_order = 0) {

  if (n_ts == 1) {
    t_len <- length(data)
  } else {
    t_len <- nrow(data)
  }

  #If grid is not supplied, we construct it by default
  if (is.null(grid)) {
    grid <- construct_grid(t_len)
  }

  #If ijset is not supplied, we compare all
  #possible pairs of time series.
  if (is.null(ijset)) {
    ijset <- expand.grid(i = 1:n_ts, j = 1:n_ts)
    ijset <- ijset[ijset$i < ijset$j, ]
  }


  #These full grids we need for plotting the SiZer maps
  gset_full    <- grid$gset_full
  u_grid_full  <- unique(gset_full[, 1])
  pos_full     <- grid$pos_full

  test_res     <- rep(2, length(pos_full))

  # Select (1-alpha) quantile of the multiscale statistic under the null
  quantiles <- compute_quantiles(t_len = t_len, grid = grid, n_ts = n_ts,
                                 ijset = ijset, sigma = sigma,
                                 sim_runs = sim_runs,
                                 deriv_order = deriv_order)

  probs <- as.vector(quantiles$quant[1, ])
  quant <- as.vector(quantiles$quant[2, ])

  if (sum(probs == (1 - alpha)) == 0)
    pos <- which.min(abs(probs - (1 - alpha)))
  if (sum(probs == (1 - alpha)) != 0)
    pos <- which.max(probs == (1 - alpha))

  quant <- quant[pos]

  psi   <- compute_statistics(data = data, sigma = sigma, n_ts = n_ts,
                              grid = grid, deriv_order = deriv_order)
  stat  <- psi$stat

  if (n_ts == 1) {
    gset_with_vals             <- psi$gset_with_vals
    test_results               <- (gset_with_vals$vals_cor > quant) * sign(gset_with_vals$vals)
    test_res[!is.na(pos_full)] <- test_results
    test                       <- matrix(test_res, ncol = length(u_grid_full),
                                         byrow = TRUE)

    if (stat > quant) {
      cat("For the given time series given we reject H_0 with probability",
          alpha, ". Psihat_statistic = ", stat,
          ". Gaussian quantile value = ", quant, "\n")
    } else {
      cat("For the given time series we fail to reject H_0 with probability",
          alpha, ". Psihat_statistic = ", stat,
          ". Gaussian quantile value = ", quant, "\n")
    }

    gset_with_vals$test <- test_results

    return(list(quant = quant, stat = psi$stat, test_matrix = test,
                gset_with_vals = gset_with_vals))

  } else {
    if (stat > quant) {
      cat("We reject H_0 with probability", alpha, ". Psihat_statistic = ",
          stat, ". Number of pairwise rejections = ",
          sum(stat > quant, na.rm = TRUE),
          ". Gaussian quantile value = ", quant, "\n")
    } else {
      cat("We fail to reject H_0 with probability", alpha,
          ". Psihat_statistic = ", stat,
          ". Gaussian quantile value = ", quant, "\n")
    }

    gset_with_values <- psi$gset_with_values

    for (i in seq_len(nrow(ijset))) {
      test_results               <- gset_with_values[[i]]$vals > quant
      gset_with_values[[i]]$test <- test_results
    }

    return(list(quant = quant, stat = stat, stat_pairwise = psi$stat_pairwise,
                gset_with_values = gset_with_values, ijset = ijset))
  }
}
