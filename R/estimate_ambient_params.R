#' Estimate Ambient Parameters from Contamination
#'
#' estimate_ambient_params will take a snRNA-seq contamination dataset as input and estimate parameters for downstream analysis.
#'
#' @param contamination An n x j matrix-like object where n = number of genes and j = number of ambient mRNA droplets. The droplets in this matrix are believed to be by the user to be composed of ambient mRNA. Including droplets beneath a specific UMI cutoff (i.e. 100) is a basic idea to get ambient mRNA droplets.
#'
#' @param bgBatch A vector of length j corresponding to batch information of the ambient mRNA droplets in \code{contamination}.
#'
#' @return A matrix with the unnormalized psuedobulk of the ambient profile to be used by the Multinomial distribution to simulate ambient counts.
#'
#' @export
#'
#' @import Matrix
#' @import sparseMatrixStats

estimate_ambient_params <- function(contamination, bgBatch = NULL) {

  num_genes <- nrow(contamination) # How many genes did the user provide in their dataset?

  if (is.null(bgBatch)) { # Metadata-free case. Get rowMeans and UMIs per droplet and return.
    genedf <- sparseMatrixStats::rowSums2(contamination) # Mean expression per gene
    names(genedf) <- rownames(contamination)
    return(matrix(genedf, nrow = num_genes))
  } else { # We want to have different ambient profiles for different batches, should the user provide that information.
    batches <- unique(bgBatch)
    num_batches <- length(batches)
    genedf <- matrix(NA, ncol = num_batches, nrow = nrow(contamination))
    rownames(genedf) <- rownames(contamination)
    for (i in 1:length(batches)) {
      genedf[,i] <- sparseMatrixStats::rowSums2(contamination[,which(bgBatch == batches[i])])
    }
    result <- matrix(genedf, nrow = num_genes)
    colnames(result) <- as.character(batches)
    return(result)
  }
}
