#' Multiscale Inference for Nonparametric Regression(s)
#'
#' This package performs a multiscale analysis of a nonparametric regression
#' or nonparametric regressions with time series errors.
#'
#' @description This package performs a multiscale analysis of a nonparametric
#'   regression or nonparametric regressions with time series errors.
#'
#'   In case of a single nonparametric regression, the multiscale method to
#'   test qualitative hypotheses about the nonparametric time trend \eqn{m}
#'   in the model \eqn{Y_t = m(t/T) + \epsilon_t} with time series errors
#'   \eqn{\epsilon_t} is provided. The method was first proposed in
#'   Khismatullina and Vogt (2019). It allows to test for shape properties
#'   (areas of monotonic decrease or increase) of the trend \eqn{m}.
#'
#'   This method require an estimator of the long-run error variance
#'   \eqn{\sigma^2 = \sum_{l=-\infty}^{\infty} Cov(\epsilon_0, \epsilon_l)}.
#'   Hence, the package also provides the difference-based
#'   estimator for the case that the errors belong to the class of
#'   \eqn{AR(\infty)} processes. The estimator was first proposed in
#'   Khismatullina and Vogt (2019).
#'
#'   In case of multiple nonparametric regressions, the multiscale method
#'   to test qualitative hypotheses about the nonparametric time trends
#'   in the context of epidemic modelling is provided. Specifically,
#'   we assume that the we observe a sample of the count data
#'   \eqn{\{\mathcal{X}_i = \{ X_{it}: 1 \le 1 \le T \}\}}, where \eqn{X_{it}}
#'   are quasi-Poisson distributed with time-varying intensity parameter
#'   \eqn{\lambda_i(t/T)}. The multiscale method allows to test whether
#'   intenisty parameters are different or not, and if they are, it detects
#'   the regions where the trends significantly differ from each other.
#'   The method was first proposed in Khismatullina and Vogt (2020) and can be
#'   used for comparing the rates of infection of COVID-19 across countries.
#'
#' @docType package
#' @name multiscale-package
#' @aliases multiscale
#' @seealso
#'    \url{https://rss.onlinelibrary.wiley.com/doi/full/10.1111/rssb.12347}
#' @references Khismatullina M., Vogt M. Multiscale inference and long-run
#'             variance estimation in non-parametric regression with
#'             time series errors //Journal of the Royal Statistical Society:
#'             Series B (Statistical Methodology). - 2019.
NULL
