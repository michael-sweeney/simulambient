#'
#'
#'
#'
#'
#'
#' @export
#'
#' @import stats
#' @import Matrix

simulate_counts <- function(geneMeans, geneSDs, numNucleiPerGroup, seed) {

  set.seed(seed)

  tol <- 1e-10 # Small value to add to avoid NAs

  num_genes <- nrow(geneMeans)
  num_groups <- ncol(geneMeans)

  geneCVs <- ((geneSDs + tol) / geneMeans)^2

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
