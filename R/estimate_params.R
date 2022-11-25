#' Estimate Parameters from Dataset
#'
#' estimate_params will take a snRNA-seq dataset as input and estimate parameters for downstream analysis.
#'
#' @param dataset An n x k dataset where n = number of genes and k = number of nuclei. Assumed to be demultiplexed and cell type-annotated.
#'
#' @param cellTypes A vector of length k that provides the cell type of each nucleus in \code{dataset}. If all nuclei are from one type or cell type is not to be accounted for, parameter can be left blank.
#'
#' @param batch A vector of length k that provides the batch of the nucleus in \code{dataset}. If there is no batch effect, parameter can be left blank.
#'
#' @param sample A vector of length k that provides the associated sex of each droplet in \code{dataset}.
#'
#' @return Idk yet
#'
#' @export
#'

estimate_params <- function(dataset, cellTypes = NULL, batch = NULL, sample = NULL) {

  num_genes <- nrow(dataset)
  num_nuclei <- ncol(dataset)
  metadata <- c()

  if (is.null(rownames(dataset))) {
    stop("Single-nucleus data should have row names to identify genes. Please attribute gene names to data before continuing.")
  }

  if (!is.null(cellTypes)) {
    if (length(cellTypes) != num_nuclei) {
      stop("Each nucleus should be assigned to a cell type. The cell type metadata does not match the number of nuclei.")
    }
    metadata <- c(metadata, "cellTypes")
  }

  if (!is.null(batch)) {
    if (length(batch) != num_nuclei) {
      stop("Each nucleus should be assigned to a batch. The batch metadata does not match the number of nuclei.")
    }
    metadata <- c(metadata, "batch")
  }

  if (!is.null(sample)) {
    if (length(sample) != num_nuclei) {
      stop("Each nucleus should be assigned to a sample. The sample metadata does not match the number of nuclei.")
    }
    metadata <- c(metadata, "sample")
  }

  umis <- colSums()
  if (length(metadata) == 0) {
    genedf <- rowMeans(dataset)
  } else {
    metadata_df <- data.frame(mget(metadata))
    metadata_string <- paste(metadata, collapse = "+")
    counts_with_metadata <- t(rbind(dataset, metadata_df))
    counts_with_metadata <- aggregate(eval(parse(text = paste(".~", metadata_string, sep = ""))),
                                      data = counts_with_metadata, FUN = mean)
  }
}
