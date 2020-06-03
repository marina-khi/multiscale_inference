#' Computes the location-bandwidth grid for the multiscale test.
#'
#' @param t          Sample size.
#' @param u_grid     Vector of location points in the unit interval [0,1].
#'                   If u.grid=NULL, a default grid is used.
#' @param h_grid     Vector of bandwidths, each bandwidth is supposed to lie
#'                   in (0,0.5). If h.grid=NULL, a default grid is used.
#' @param deletions  Logical vector of the length len(u.grid) * len(h.grid).
#'                   Each element is either TRUE, which means that
#'                   the corresponding location-bandwidth point (u, h)
#'                   is NOT deleted from the grid, or FALSE, which means that
#'                   the corresponding location-bandwidth point (u, h) IS
#'                   deleted from the grid. Default is deletions = NULL
#'                   in which case nothing is deleted.
#'                   See vignette for the use.
#' @return gset      Matrix of location-bandwidth points (u,h) that remains
#'                   after deletions, the i-th row gset[i,] corresponds to
#'                   the i-th point (u,h).
#' @return bws       Vector of bandwidths (after deletions).
#' @return lens      Vector of length=length(bws), lens[i] gives the number of
#'                   locations in the grid for the i-th bandwidth level.
#' @return gtype     Type of grid that is used, either 'default' or
#'                   'non-default'.
#' @return gset_full Matrix of all location-bandwidth pairs (u,h) including
#'                   deleted ones.
#' @return pos_full  Logical vector indicating which points (u,h) have been
#'                   deleted.
#'
#' @export
#'
#' @examples
#' construct_grid(100)
#' construct_grid(100, u_grid = seq(from = 0.05, to = 1, by = 0.05),
#'                h_grid = c(0.1, 0.2, 0.3, 0.4))


construct_grid <- function(t, u_grid = NULL, h_grid = NULL, deletions = NULL) {

  grid_type <- "non-default"
  if (is.null(u_grid) & is.null(h_grid)) {
    grid_type <- "default"
    t_len  <- t
    u_grid <- seq(from = 5 / t_len, to = 1, by = 5 / t_len)
    h_grid <- seq(from = 5 / t_len, to = 1 / 4, by = 5 / t_len)
    h_grid <- h_grid[h_grid > log(t_len) / t_len]
  }

  gset      <- expand.grid(u = u_grid, h = h_grid)
  gset_full <- gset
  n         <- dim(gset)[1]
  pos_full  <- rep(TRUE, n)

  if (!is.null(deletions))
    pos_full <- deletions

  gset <- gset[pos_full, ]
  bws  <- unique(gset[, 2])
  lengths_u <- rep(NA, length(bws))
  for (i in seq_len(length(bws)))
     lengths_u[i] <- sum(gset[, 2] == bws[i])

  return(list(gset = gset, bws = bws, lens = lengths_u, gtype = grid_type,
              gset_full = gset_full, pos_full = pos_full))
}
