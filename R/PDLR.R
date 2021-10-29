#' Calculates the paired density logratio (PDLR) between two multi-dimensional sets of points
#' @title Paired density logratio
#' @param X1 A matrix of points, with each row corresponding to a point's coordinates in a multi-dimensional real space.
#' @param X2 A second matrix of points, to be compared to X1. The number of columns and rows must be equal to those of X1.
#' @param K The number of the first K nearest neighbours (NN) to average over, defaults to K = 1
#' @return A set of real numbers corresponding to the log-10 ratio of the cross NN distances to the self NN distances
#' @export

PDLR <- function(X1, X2, K = 1) {

  if (nrow(X1) != nrow(X2) | ncol(X1) != ncol(X2)) stop("Dimensions of X1 and X2 must be identical")

  if (K == 1) {
    n1 <- nn2(X1, eps = 0, k = 2)$nn.dists[, 2]
    n2 <- nn2(X2, X1, eps = 0, k = 2)$nn.dists[, 1]
  } else {
    n1 <- rowMeans(nn2(X1, eps = 0, k = K+1)$nn.dists[, 2:(K+1)])
    n2 <- rowMeans(nn2(X2, X1, eps = 0, k = K)$nn.dists[, 1:K])

  }
  ratio <- log10(n2/n1)
  return(ratio)
}
