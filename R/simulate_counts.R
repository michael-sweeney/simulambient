#' Simulate snRNA-seq Counts with Negative Binomial
#'
#' simulate_counts simulates snRNA-seq counts given parameters for the negative binomial.
#'
#' @param geneMeans Matrix-like object of gene expression means for each metadata group. (Rows = genes, columns = metadata groups)
#'
#' @param geneSDs Matrix-like object of gene expression standard devitaitons for each metadata group. (Rows = genes, column = metadata groups)
#'
#' @param numNucleiPerGroup Vector of integers correpsonding to how many nuclei need to be simulated within each metadata group.
#'
#' @param seed Seed for random number generation. Can be set by the user. Defaults to 615.
#'
#' @export
#'
#' @import stats
#' @import Matrix

simulate_counts <- function(geneMeans, geneSDs, numNucleiPerGroup, seed) {

  set.seed(seed)

  tol <- 1e-10 # Small value to add to avoid NAs/NaNs

  num_genes <- nrow(geneMeans)
  num_groups <- ncol(geneMeans)

  geneCVs <- ((geneSDs + tol) / geneMeans)^2 # Calculate

  outputMatrix <- NULL

  for (i in 1:num_groups) {
    df <- data.frame(geneMeans[,i], geneCVs[,i])
    function_simulate <- function(x, y) stats::rnbinom(numNucleiPerGroup[i], mu = x, size = y)
    if (numNucleiPerGroup[i] == 1) {
      counts <- Matrix(mapply(function_simulate, df[,1], df[,2]), sparse = TRUE)
    } else {
      counts <- Matrix(t(mapply(function_simulate, df[,1], df[,2])), sparse = TRUE)
    }
    if (is.null(outputMatrix)) {
      outputMatrix <- counts
    } else {
      outputMatrix <- cbind(outputMatrix, counts)
    }
  }
  rownames(outputMatrix) <- rownames(geneMeans)
  return(outputMatrix)
}
