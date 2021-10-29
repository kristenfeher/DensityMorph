#' Given a set of matrices with idential dimension, calculate a pairwise matrix of quasi-distances between them using DensityMorph
#' @title Calculate latent space coordinates
#' @param D A list of matrices of identical dimension. The points are stored as rows in each matrix.
#' @param Xref A pre-calculated self-PDLR distribution
#' @param K The number of the first K nearest neighbours (NN) to average over, defaults to K = 1
#' @return A named list containing a symmetrical matrix containing the quasi-distances calculated using DensityMorph (the dimension is equal to the number of matrices in D) and the latent space with coordinates stored in rows
#' @examples
#' p <- 3
#' Ntype <- 10

#' set.seed(21)
#' library(mvtnorm)
#' experiment <- replicate(Ntype, {
#'   nclus <- 5;
#'  clus_means <- t(sapply(1:nclus, function(i) runif(p, 0, 5)));
#'  clus_var <- runif(nclus, 1, 3);
#'  clus_cor <- runif(nclus, -0.3, 0.3);
#'  clus_count <- runif(nclus, 1, 10);
#'  N <- rmultinom(1, 1000, clus_count);
#'
#'  sample <- lapply(1:nclus, function(i) {
#'    M <- as.vector(rmvnorm(1, mean = clus_means[i, ], sigma = diag(rep(0.1, p))))
#'    S <- rnorm(1, mean = clus_var[i], sd = 0.05)
#'    C <- rnorm(1, mean = clus_cor[i], sd = 0.05)
#'    sig <- matrix(C, nrow = p, ncol = p)
#'    diag(sig) <- S
#'    pop <- rmvnorm(N[i], M, sig)
#'    return(list(pop = pop, M = M, S = S, C = C, N = N[i]))
#'  }
#'  );
#'  s <- mapply(function(x) x$pop, sample, SIMPLIFY = FALSE);
#'  s <- do.call(rbind, s);
#'  params <- data.frame(t(mapply(function(x) unlist(x[2:5]), sample)));
#'  reps <- list(data = s, params = params)
#'
#'}, simplify = FALSE
#')
#'expr_data <- mapply(function(x) x$data, experiment, SIMPLIFY = FALSE)

#'Xref <- reference_distribution(p, 10^5, 10^3)
#'# Dimension reduction should be performed at this step if necessary.
#'# Dimension should preferably be reduced to 3,4 or 5 to take maximum advantage of the time saving.
#'# The same transformation should be applied to all datasets, i.e. calculated over a pooled dataset.
#'# One possibility is to use spectral map analysis and it is implemented in the function SMA.
#'LS <- latent_space(expr_data, Xref, 1)
#'bulk_data <- bulk(expr_data, c(3, 3, 3), c('A', 'B', 'C'))
#'EC_A <- EC(LS$latent_space[, 1:3], bulk_data[[1]][, 1])
#'
#'plot(LS$latent_space[, 1:2], pch = 16, cex = 1)
#'arrows(0, 0, EC_A$proj[1]/10, EC_A$proj[2]/10)

latent_space <- function(D, Xref, K = 1) {
  distmtx <- lapply(1:(length(D) - 1), function(i) sapply((i+1):length(D), function(j) {
    max(densitymorph(D[[i]], D[[j]], Xref, K))
  } ))
  dm <- mapply(function(x) c(rep(0, length(D) - length(x)), x), distmtx)
  dm <- cbind(dm, 0)
  dm <- dm + t(dm)
  LS <- ade4::dudi.pco(as.dist(dm^0.5), nf = 3, scannf = FALSE)$tab

  return(list(dist_mtx = dm, latent_space = LS))
}
