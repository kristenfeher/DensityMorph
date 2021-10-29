#' Takes a set of single cell samples and extracts three types of univariate 'bulk' measurements from them: proportion of cells positive for each marker, average MFI/expression of the positive cells for each marker, and the log-odds ratio of double positive cells for each marker pair.
#' @title Calculate bulk quantities
#' @param D A list of single cell datasets of equal dimension
#' @param threshold A vector of expression thresholds, one for each marker, above which a cell is deemed to be positive for that marker. Currently you need to supply a threshold for each marker, there is no functionality to select a subset of markers
#' @param marker_names A vector of marker names, of equal length to threshold.
#' @return A named list of matrices containing 'bulk' measurements for each sample and marker/marker pair.
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


bulk <- function(D, threshold, marker_names) {

  markers <- mapply(function(x) sign(sweep(x, 2, threshold, '-')), D, SIMPLIFY = FALSE)

  markers_1 <- t(mapply(function(x) apply(x, 2, function(y) length(which(y > 0)))/nrow(x), markers))
  markers_expr <- do.call(rbind, lapply(1:length(markers), function(i) colSums((as.matrix(markers[[i]])+1) * as.matrix(D[[i]])) / apply(markers[[i]], 2, function(q) length(which(q == 1)))))

  colnames(markers_1) <- marker_names
  colnames(markers_expr) <- marker_names

  markers_2 <- t(mapply(function(x) unlist(lapply(1:(ncol(x) -1), function(i) lapply((i+1):ncol(x), function(j) {
    tmp <- table(x[, i], x[, j]); (tmp[1, 1]*tmp[2, 2])/(tmp[1, 2] * tmp[2, 1])
  }))), markers))
  markers_2[which(is.infinite(markers_2), arr.ind = TRUE)] <- max(markers_2[!is.infinite(markers_2)])*2
  markers_2[which(markers_2 == 0, arr.ind = TRUE)] <- min(markers_2[markers_2 > 0])/2
  markers_2 <- log10(markers_2)

  L <- unlist(sapply(1:(length(marker_names) - 1), function(i) sapply((i+1):length(marker_names), function(j) {
    paste(marker_names[i], marker_names[j], sep = '_')
  })))
  colnames(markers_2) <- L

  return(list(proportion = markers_1, mean_expr = markers_expr, LOR = markers_2))

}
