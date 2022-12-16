#' Check User Arguments
#'
#' check_args will take any cellTypes, batch, and sex arguments provided by the user and check them for validity.
#'
#' @param dataset An n x k matrix-like object where n = number of genes and k = number of nuclei to be decontaminated.
#'
#' @param contamination An n x j matrix-like object where n = number of genes and j = number of ambient mRNA droplets. The droplets in this matrix are believed to be by the user to be composed of ambient mRNA. Including droplets beneath a specific UMI cutoff (i.e. 100) is a basic idea to get ambient mRNA droplets.
#'
#' @param metadata A matrix with k rows that corresponds to metadata of the k nuclei provided in \code{dataset}. Each column should represent a different aspect of metadata, i.e. sex, cell types, etc.

#' @param cellTypes A vector of length k should cell type information not be provided in \code{metadata}. If it is provided, this parameter should be a string of the column name that represents the cell type information in \code{metadata}. If cell type information is not desired to be included, default is NULL.
#'
#' @param batch A vector of length k that provides the batch of the nucleus in \code{dataset}. If there is no batch effect, parameter can be left blank.
#'
#' @param bgBatch A vector of length j corresponding to batch information of the ambient mRNA droplets in \code{contamination}.
#'
#' @return Validated metadata provided by the user.
#'
#' @export

check_args <- function(dataset, contamination = NULL, metadata = NULL, cellTypes = NULL, batch = NULL, bgBatch = NULL) {

  num_nuclei <- ncol(dataset)
  num_ambient <- ncol(contamination)
  metadata <- metadata

  # Throw error if genes aren't named.
  if (is.null(rownames(dataset))) {
    stop("Single-nucleus data should have row names to identify genes. Please attribute gene names to data before continuing.")
  }

  if (!is.null(contamination)) {
    if (is.null(rownames(contamination))) {
      stop("Ambient mRNA droplet data should have row names to identify genes. Please attribute gene names to data before continuing.")
    }
  }

  if (!all(rownames(dataset) == rownames(contamination))) {
    stop("Ambient mRNA droplet data should have the exact same genes as the single-nucleus data. Droplets from the same 10X experiment should have the same genes (in features.tsv). If there are differences between the single-nucleus data and the ambient mRNA profile, it should be manually updated by the user such that the row names of these two matrices are the same. Row names must be in the same order.")
  }

  # Throw error if user provides cell type information but metadata length differs from number of nuclei.
  if (!is.null(metadata)) {
    if (nrow(metadata) != num_nuclei) {
      stop("The rows of the metadata matrix does not match the number of nuclei.")
    }
    if (is.null(colnames(metadata))) {
      stop("Please include column names to the metadata.")
    }
  }

  if (!is.null(cellTypes)) {
    if (length(cellTypes) == 0) {
      stop("Cell types should either be NULL, a string, or a vector of length k.")
    } else if (length(cellTypes) == 1) { # String case
      if (cellTypes %in% colnames(metadata) == FALSE) {
        stop("cellTypes does not match a column name in metadata. Please make sure cellTypes perfectly corresponds to a column name of metadata.")
      } else {
        colnames(metadata)[which(colnames(metadata) == cellTypes)] <- "cellTypes"
      }
    } else { # Vector case
      if (length(cellTypes) != num_nuclei) {
        stop("The cell type metadata differs in length from the number of nuclei. Please make sure cellTypes is length of k.")
      } else {
        if (!is.null(metadata)) {
          metadata <- cbind(metadata, cellTypes)
          colnames(metadata)[ncol(metadata)] <- "cellTypes"
        } else {
          metadata <- as.data.frame(cellTypes)
          colnames(metadata) <- "cellTypes"
        }
      }
    }
  }

  # Throw error if user provides batch effect information but batch labels aren't same length as number of nuclei.
  if (!is.null(batch)) {
    if (length(batch) != num_nuclei) {
      stop("Each nucleus should be assigned to a batch. The batch metadata does not match the number of nuclei.")
    }
  }

  if ((is.null(bgBatch) == TRUE && is.null(batch) == FALSE) || (is.null(bgBatch) == FALSE && is.null(batch) == TRUE)) {
    stop("Batch metadata is required for the ambient mRNA profile if it is provided for the single-nucleus data and vice versa.")
  }

  if (all(bgBatch %in% batch) == FALSE || all(batch %in% bgBatch) == FALSE) {
    stop("Ambient batch information should be a subset of the nucleus batch information. Please do not include ambient information from batches that were not present in the experiment.")
  }

  if (!is.null(bgBatch)) {
    if (length(bgBatch) != num_ambient) {
      stop("Each ambient mRNA droplet should be assigned to a batch. The ambient mRNA profile batch metadata does not match the number of ambient mRNA droplets.")
    }
  }
  return(metadata)
}
