#' Computes autocovariances at lags 0 to p for the ell-th differences of data.
#'
#' @keywords internal
#' @param data       Time series for which we calculate autocovariances.
#' @param ell        Order of differences used.
#' @param p          Maximum lag for the autocovariances
#' @return autocovs  Vector of length (p + 1) that consists of empirical
#'                   autocoavariances fr the corresponding lag - 1.

emp_acf <- function(data, ell, p) {
  y_diff   <- diff(data, ell)
  len      <- length(y_diff)
  autocovs <- rep(0, p + 1)
  for (k in seq_len(p + 1)) {
    autocovs[k] <- sum(y_diff[k:len] * y_diff[1:(len - k + 1)]) / len
  }
  return(autocovs)
}

#' Computes variance of AR(p) innovation terms eta.
#'
#' @keywords internal
#' @param data      Time series.
#' @param coefs     Estimated coefficients.
#' @param p         AR order of the time series.
#' @return var.eta  Variance of the innovation term

variance_eta <- function(data, coefs, p) {
  y_diff <- diff(data)
  len    <- length(y_diff)
  x_diff <- matrix(0, nrow = len - p, ncol = p)
  for (i in seq_len(p)) {
     x_diff[, i] <- y_diff[(p + 1 - i):(len - i)]
  }
  y_diff  <- y_diff[(p + 1):len]
  resid   <- y_diff - as.vector(x_diff %*% coefs)
  var_eta <- mean(resid^2) / 2
  return(var_eta)
}

#' Computes vector of correction terms for second-stage estimator
#' of AR parameters as described in Khismatullina, Vogt (2019).
#'
#' @keywords internal
#' @param coefs     Given (estimated) coefficients of the AR time series.
#' @param var.eta   Variance of the innovation term
#' @param len       Length of the vector of the corrections
#' @return c.vec    Vector of the corrections terms of length len.
corrections <- function(coefs, var_eta, len) {
  p        <- length(coefs)
  c_vec    <- rep(0, len)
  c_vec[1] <- 1
  for (j in 2:len) {
    lags     <- (j - 1):max(j - p, 1)
    c_vec[j] <- sum(coefs[1:length(lags)] * c_vec[lags])
  }
  c_vec <- c_vec * var_eta
  return(c_vec)
}

#' Computes estimator of the AR(p) coefficients by the procedure from
#' Khismatullina and Vogt (2019).
#'
#' @keywords internal
#' @param data      Time series.
#' @param l1,l2     Tuning parameters.
#' @param correct   Vector of the corrections, either zero or calculated by
#'                  the function \code{\link{corrections}}.
#' @param p         AR order of the time series.
#' @return a_hat    Vector of length p  of estimated AR coefficients.
ar_coef <- function(data, l1, l2, correct, p) {
  pos <- 0
  a_mat <- matrix(0, nrow = l2 - l1 + 1, ncol = p)
  for (ell in l1:l2) {
    pos <- pos + 1
    autocovs <- emp_acf(data, ell, p)
    cov_mat  <- matrix(0, ncol = p, nrow = p)
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
         cov_mat[i, j] <- autocovs[abs(i - j) + 1]
      }
    }
    if (ell >= p) {
      correct_vec <- correct[(ell - 1 + 1):(ell - p + 1)]
    } else {
      correct_vec <- c(correct[(ell - 1 + 1):1], rep(0, p - ell))
    }
    cov_vec <- autocovs[2:(p + 1)] + correct_vec
    a_mat[pos, ] <- solve(cov_mat) %*% cov_vec
  }
  a_hat <- colMeans(a_mat)
  a_hat <- as.vector(a_hat)
  return(a_hat)
}

#' Computes estimator of the long-run variance of the error terms.
#'
#' @description     A difference based estimator for the coefficients and
#'                  long-run variance in case of a nonparametric regression
#'                  model \eqn{Y(t) = m(t/T) + \epsilon(t)} where the errors are AR(p).
#'                  The procedure was first introduced in
#'                  Khismatullina and Vogt (2019).
#'
#' @param data      A vector of Y(1), Y(2), ... Y(T).
#' @param q,r_bar   Integers, tuning parameters.
#' @param p         AR order of the error terms.
#' @return lrv      Estimator of the long run variance of the error terms.
#' @return ahat     Vector of length p of estimated AR coefficients.
#' @return vareta   Estimator of the variance of the innovation term
#'
#' @export
estimate_lrv <- function(data, q, r_bar, p) {
    a_tilde       <- ar_coef(data = data, l1 = q, l2 = q,
                             correct = rep(0, max(p, q) + 1), p = p)
    sig_eta_tilde <- variance_eta(data = data, coefs = a_tilde, p = p)
    correct       <- corrections(coefs = a_tilde, var_eta = sig_eta_tilde,
                                 len = r_bar + 1)
    a_hat         <- ar_coef(data = data, l1 = 1, l2 = r_bar,
                             correct = correct, p = p)
    sig_eta_hat   <- variance_eta(data = data, coefs = a_hat, p = p)
    lrv_hat       <- sig_eta_hat / (1 - sum(a_hat))^2
    return(list(lrv = lrv_hat, ahat = a_hat, vareta = sig_eta_hat))
}
