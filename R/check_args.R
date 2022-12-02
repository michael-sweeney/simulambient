#' Check User Arguments
#'
#' check_args will take any cellTypes, batch, and sex arguments provided by the user and check them for validity.
#'
#' @param dataset An n x k matrix-like object where n = number of genes and k = number of nuclei to be decontaminated.
#'
#' @param contamination An n x j matrix-like object where n = number of genes and j = number of ambient mRNA droplets. The droplets in this matrix are believed to be by the user to be composed of ambient mRNA. Including droplets beneath a specific UMI cutoff (i.e. 100) is a basic idea to get ambient mRNA droplets.
#'
#' @param cellTypes A vector of length k that provides the cell type of each nucleus in \code{dataset}. If all nuclei are from one type or cell type is not to be accounted for, parameter can be left blank.
#'
#' @param batch A vector of length k that provides the batch of the nucleus in \code{dataset}. If there is no batch effect, parameter can be left blank.
#'
#' @param bgBatch A vector of length j corresponding to batch information of the ambient mRNA droplets in \code{contamination}.
#'
#' @param sex A vector of length k that provides the associated sex of each droplet in \code{dataset}.
#'
#' @return A list of strings containing the valid metadata provided by the user.
#'
#' @export

check_args <- function(dataset, contamination = NULL, cellTypes = NULL, batch = NULL, bgBatch = NULL, sex = NULL) {

  num_nuclei <- ncol(dataset)
  num_ambient <- ncol(contamination)
  metadata <- c()

  # Throw error if genes aren't named.
  if (is.null(rownames(dataset))) {
    stop("Single-nucleus data should have row names to identify genes. Please attribute gene names to data before continuing.")
  }

  if (is.null(rownames(contamination))) {
    stop("Ambient mRNA droplet data should have row names to identify genes. Please attribute gene names to data before continuing.")
  }

  if (!all(rownames(dataset) == rownames(contamination))) {
    stop("Ambient mRNA droplet data should have the exact same genes as the single-nucleus data. Droplets from the same 10X experiment should have the same genes (in features.tsv). If there are differences between the single-nucleus data and the ambient mRNA profile, it should be manually updated by the user such that the row names of these two matrices are the same.")
  }

  # Throw error if user provides cell type information but cell type labels aren't same length as number of nuclei.
  if (!is.null(cellTypes)) {
    if (length(cellTypes) != num_nuclei) {
      stop("Each nucleus should be assigned to a cell type. The cell type metadata does not match the number of nuclei.")
    }
    metadata <- c(metadata, "cellTypes")
  }

  # Throw error if user provides batch effect information but batch labels aren't same length as number of nuclei.
  if (!is.null(batch)) {
    if (length(batch) != num_nuclei) {
      stop("Each nucleus should be assigned to a batch. The batch metadata does not match the number of nuclei.")
    }
    metadata <- c(metadata, "batch")
  }

  if (is.null(bgBatch) == TRUE && is.null(batch) == FALSE) {
    stop("Batch metadata is required for the ambient mRNA profile if it is provided for the single-nucleus data.")
  }

  if (!is.null(bgBatch)) {
    if (length(bgBatch) != num_ambient) {
      stop("Each ambient mRNA droplet should be assigned to a batch. The ambient mRNA profile batch metadata does not match the number of ambient mRNA droplets.")
    }
  }

  # Throw error if user provides sex information but labels aren't same length as number of nuclei.
  if (!is.null(sex)) {
    if (length(sex) != num_nuclei) {
      stop("Each nucleus should be assigned to a sex. The sex metadata does not match the number of nuclei.")
    }
    metadata <- c(metadata, "sex")
  }

  return(metadata)

}
