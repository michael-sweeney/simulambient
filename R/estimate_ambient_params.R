#' Estimate Ambient Parameters from Contamination
#'
#' estimate_ambient_params will take a snRNA-seq contamination dataset as input and estimate parameters for downstream analysis.
#'
#' @param contamination An n x j matrix-like object where n = number of genes and j = number of ambient mRNA droplets. The droplets in this matrix are believed to be by the user to be composed of ambient mRNA. Including droplets beneath a specific UMI cutoff (i.e. 100) is a basic idea to get ambient mRNA droplets.
#'
#' @param bgBatch A vector of length j corresponding to batch information of the ambient mRNA droplets in \code{contamination}.
#'
#' @return A list containing two values. The first value is a matrix \code{geneMeans} that lists the mean counts value for each gene (for that particular metadata, if applicable). The second value is a matrix \code{geneMeans} that lists the standad deviation for each gene (for that particular metadata, if applicable).
#'
#' @export
#'
#' @import Matrix
#' @import sparseMatrixStats
#' @import methods

estimate_ambient_params <- function(contamination, bgBatch = NULL) {

  num_genes <- nrow(contamination) # How many genes did the user provide in their dataset?

  if (is.null(bgBatch)) { # Metadata-free case. Get rowMeans and UMIs per droplet and return.
    genedf <- sparseMatrixStats::rowSums2(contamination) # Mean expression per gene
    names(genedf) <- rownames(contamination)
    return(matrix(genedf, nrow = num_genes))
  } else {
    batches <- unique(bgBatch)
    num_batches <- length(batches)
    genedf <- matrix(NA, ncol = num_batches, nrow = nrow(contamination))
    for (i in batches) {
      genedf[,i] <- sparseMatrixStats::rowSums2(contamination[,which(bgBatch == i)])
    }
    result <- matrix(genedf, nrow = num_genes)
    colnames(result) <- as.character(batches)
    return(result)
  }
}
