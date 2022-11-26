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
#' @return A list containing three values.
#'
#' @export
#'
#' @import Matrix
#' @import data.table

estimate_params <- function(dataset, cellTypes = NULL, batch = NULL, sample = NULL) {

  dataset <- as(dataset, "CsparseMatrix") # Convert data to sparse format

  num_genes <- nrow(dataset) # How many genes did the user provide in their dataset?
  num_nuclei <- ncol(dataset) # How many nuclei did the user provide in their dataset?
  metadata <- c() # Vector of strings. Could contain none, "cellTypes", "batch", "sex", or any combination depending on metadata provided by the user.

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

  umis <- colSums(dataset) # UMIs per sample
  if (length(metadata) == 0) {
    genedf <- rowMeans(dataset)
    names(genedf) <- rownames(dataset)
    return(list(genes = genedf, tbl = NULL, umis = umis))
  } else {
    metadata_df <- data.frame(mget(metadata))
    for (i in 1:ncol(metadata_df)) {
      metadata_df[,i] <- as.factor(metadata_df[,i])
      levels(metadata_df[,i]) <- 1:length(levels(metadata_df[,i]))
      metadata_df[,i] <- as.numeric(levels(metadata_df[,i]))[metadata_df[,i]]
    }
    metadata_df <- Matrix(as.matrix(metadata_df), sparse = TRUE)
    metadata_string <- paste(metadata, collapse = ",")
    counts_with_metadata <- t(rbind(dataset, metadata_df))
    counts_with_metadata <- aggregate(eval(parse(text = paste(".~", metadata_string, sep = ""))),
                                      data = counts_with_metadata, FUN = mean)

    colnames(counts_with_metadata) <- metadata # might be redundant

    if (!is.null(cellTypes)) {
      unique_cell_types <- unique(cellTypes)
      counts_with_metadata$cellTypes <- unique_cell_types[counts_with_metadata$cellTypes]
    }

    if (!is.null(batch)) {
      unique_batches <- unique(batch)
      counts_with_metadata$batch <- unique_batches[counts_with_metadata$batch]
    }

    if (!is.null(sample)) {
      unique_samples <- unique(sample)
      counts_with_metadata$sample <- unique_samples[counts_with_metadata$sample]
    }
  }

  tbl <- data.table(metadata_df)[, .N, by = metadata]

  return(genes = counts_with_metadata, tbl = tbl, umis = umis)
}
