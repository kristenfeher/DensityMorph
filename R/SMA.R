#' A function to perform Spectral Map Analysis (SMA), used for dimension reduction with an emphasis on difference between cell expression profiles, as opposed to each marker's dynamic range.
#' @title Spectral Map Analysis
#' @param X The log or logicle transformed data matrix, with rows as samples/observations and columns as features/marker expression
#' @return The transformed matrix. The first column is the first component, the second column is the second component, etc.

SMA <- function(X) {
  print("Reminder: did you log/logicle transform X?")

  SMA_rn <- t(apply(X, 1, function(x) x - mean(x)))
  SMA_pr <- prcomp(SMA_rn, scale. = FALSE, center = TRUE)
  return(SMA_pr)
}
