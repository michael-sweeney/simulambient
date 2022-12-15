#' Simulate snRNA-seq Contamination with Multinomial Distribution
#'
#' simulate_contamination simulates snRNA-seq ambient RNA contamination given parameters for the multinomial distribution.
#'
#' @param geneMeans Matrix-like object of ambient mRNA gene expression means for each metadata group. (Rows = genes, columns = metadata groups)
#'
#' @param numCounts A vector of length k corresponding to the contamination counts needed for each nucleus.
#'
#' @param batch A vector of length k corresponding to batch information of the ambient mRNA droplets in the nucleus \code{dataset}.
#'
#' @param seed Seed for random number generation. Can be set by the user. Defaults to 615.
#'
#' @return A sparse matrix of size n x k corresponding to the contamination counts to be added to each simulated nucleus.
#'
#' @export
#'
#' @import stats
#' @import Matrix

simulate_contamination <- function(geneMeans, numCounts, batch = NULL, seed) {

  set.seed(seed)
  numNuclei <- length(numCounts)

  outputMatrix <- Matrix::sparseMatrix(i = c(), j = c(), x = 0, dims = c(nrow(geneMeans), numNuclei))

  if (is.null(batch)) {
    for (i in 1:numNuclei) {
      outputMatrix[,i] <- stats::rmultinom(1, size = numCounts[i], prob = geneMeans[,1])
    }
  } else {
    batches <- colnames(geneMeans)
    for (i in 1:numNuclei) {
      ambient_batch <- which(batches == batch[i])
      outputMatrix[,i] <- stats::rmultinom(1, size = numCounts[i], prob = geneMeans[,ambient_batch])
    }
  }

  return(outputMatrix)
}
