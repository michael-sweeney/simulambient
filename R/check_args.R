#' Check User Arguments
#'
#' check_args will take any cellTypes, batch, and sample arguments provided by the user and check them for validity.
#' 
#' @param cellTypes A vector of length k that provides the cell type of each nucleus in \code{dataset}. If all nuclei are from one type or cell type is not to be accounted for, parameter can be left blank.
#'
#' @param batch A vector of length k that provides the batch of the nucleus in \code{dataset}. If there is no batch effect, parameter can be left blank.
#'
#' @param sample A vector of length k that provides the associated sex of each droplet in \code{dataset}.
#' 
#' @return A list of strings containing the valid metadata provided by the user.
#' 
#' @export

check_args <- function(cellTypes = NULL, batch = NULL, sample = NULL) {
  
  # Throw error if genes aren't named.
  if (is.null(rownames(dataset))) {
    stop("Single-nucleus data should have row names to identify genes. Please attribute gene names to data before continuing.")
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
  
  # Throw error if user provides sample information but labels aren't same length as number of nuclei.
  if (!is.null(sample)) {
    if (length(sample) != num_nuclei) {
      stop("Each nucleus should be assigned to a sample. The sample metadata does not match the number of nuclei.")
    }
    metadata <- c(metadata, "sample")
  }
  
  return(metadata)
  
}