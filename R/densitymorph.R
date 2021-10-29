#' Calculates the quasi-distance DensityMorph between two set of points of identical dimension and number
#' @title Quasi-distance DensityMorph
#' @param X1 A matrix of points, with each row corresponding to a point's coordinates in a multi-dimensional real space.
#' @param X2 A second matrix of points, to be compared to X1. The number of columns and rows must be equal to those of X1.
#' @param Xref A pre-calculated reference self-PDLR distribution
#' @param K The number of the first K nearest neighbours (NN) to average over, defaults to K = 1
#' @return A pair of real numbers measuring the quasi-distances from X1 to X2, and from X2 to X1
#' @export

densitymorph <- function(X1, X2, Xref, K) { # stands for transport distance, returns a 2*r matrix - take median of rows, then take the maximum
  if (nrow(X1) != nrow(X2) | ncol(X1) != ncol(X2)) stop("Dimensions of X1 and X2 must be identical")
  if (nrow(X1) != length(Xref)) stop("Length of Xref must have the same number of rows as X1, X2")
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  d1 <- PDLR(X1, X2, K)
  d2 <- PDLR(X2, X1, K)
  w1 <- wasserstein1d(Xref, d1)
  w2 <- wasserstein1d(Xref, d2)
  return(c(w1, w2))
  return(R)
}
