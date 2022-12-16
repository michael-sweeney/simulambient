#' Estimate Parameters from Dataset
#'
#' estimate_params will take a snRNA-seq dataset as input and estimate parameters for downstream analysis.
#'
#' @param dataset An n x k matrix-like object where n = number of genes and k = number of nuclei. Assumed to be demultiplexed and cell type-annotated.
#'
#' @param metadata A matrix containing metadata information provided by the user. Has k rows.
#'
#' @param metadata_indices A vector informing which nuclei belong to which metadata group.
#'
#' @return A list containing two values: a matrix \code{geneMeans} with mean gene expression for each gene for each metadata group, and a matrix \code{geneSDs} with the standard deviation gene expression for each gene in each metadata group.
#'
#' @export
#'
#' @import Matrix
#' @import matrixStats

estimate_params <- function(dataset, metadata = NULL, metadata_indices = NULL) {

  num_genes <- nrow(dataset) # How many genes did the user provide in their dataset?
  num_nuclei <- ncol(dataset) # How many nuclei did the user provide in their dataset?

  l <- ncol(metadata) # How many metadata parameters did the user include?

  if (is.null(metadata)) { # Metadata-free case. Get rowMeans and UMIs per droplet and return.
    genedf <- sparseMatrixStats::rowMeans2(dataset) # Mean expression per gene
    sddf <- sparseMatrixStats::rowSds(dataset) # Mean expression per gene
    names(genedf) <- rownames(dataset)
    return(list(geneMeans = matrix(genedf, nrow = num_genes), geneSDs = matrix(sddf, nrow = num_genes)))
  } else { # We want different mean and standard deviation gene expression for each metadata group.
    permutations <- unique(metadata_indices)
    num_permutations <- length(permutations)
    genedf <- matrix(NA, ncol = num_permutations, nrow = nrow(dataset))
    sddf <- matrix(NA, ncol = num_permutations, nrow = nrow(dataset))
    for (i in permutations) {
      if (length(which(metadata_indices == i)) == 1) {
        genedf[,i] <- dataset[,which(metadata_indices == i)]
        sddf[,i] <- 0
      } else {
        genedf[,i] <- sparseMatrixStats::rowMeans2(dataset[,which(metadata_indices == i)])
        sddf[,i] <- sparseMatrixStats::rowSds(dataset[,which(metadata_indices == i)])
      }
    }
    return(list(geneMeans = matrix(genedf, nrow = num_genes), geneSDs = matrix(sddf, nrow = num_genes)))
  }
}
