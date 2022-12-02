#' Simulambient DecontX Wrapper
#'
#' decontaminate will perform snRNA-seq decontamination on an snRNA-seq dataset using a contamination reference and cell types.
#'
#' @param dataset An n x k matrix-like object where n = number of genes and k = number of nuclei to be decontaminated.
#'
#' @param contamination An n x j matrix-like object where n = number of genes and j = number of ambient mRNA droplets. The droplets in this matrix are believed to be by the user to be composed of ambient mRNA. Including droplets beneath a specific UMI cutoff (i.e. 100) is a basic idea to get ambient mRNA droplets.
#'
#' @param cellTypes A vector of length k correpsonding to the cell types of the nuclei in \code{dataset}.
#'
#' @param batch A vector of length k corresponding to batch information of the nuclei in \code{dataset}.
#'
#' @param bgBatch A vector of length j corresponding to batch information of the ambient mRNA droplets in \code{contamination}.
#'
#' @param ... Any extra parameters passed to DecontX.
#'
#' @return A sparse matrix of class dgCMatrix with the decontaminated counts for the experiment.
#'
#' @export
#'
#' @import celda

decontaminate <- function(dataset, contamination = NULL, cellTypes = NULL, batch = NULL, bgBatch = NULL, ...) {
  decontaminated_object <- celda::decontX(x = dataset, background = contamination, z = cellTypes, batch = batch, bgBatch = bgBatch, verbose = FALSE, ...)
  return(decontaminated_object$decontXcounts)
}



