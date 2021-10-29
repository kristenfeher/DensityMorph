#' A function to pre-calculated the reference self-PDLR distribution
#' @title Self-PDLR distribution
#' @param dim The dimension of the points to be compared
#' @param Nin The number of points used to calculate the initial self-PDLR values, e.g. Nin = 100*Nout
#' @param Nout The number of points used to calculate the final self-PDLR values, based on Nout equally spaced quantiles. Nout should be equal to the number of points in the datasets to be compared.
#' @return A set of real numbers corresponding to a self-PDLR distribution
#' @export


reference_distribution <- function(dim, Nin, Nout) {
  G1 <- rmvnorm(Nin, mean = rep(0, dim))
  G2 <- rmvnorm(Nin, mean = rep(0, dim))
  Xref <- PDLR(G1, G2, K = 1)
  Xref <- quantile(Xref, seq(0, 1, length.out = Nout))
  return(Xref)
}
