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
#' @return A list containing two values.
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

  l <- length(metadata) # How many metadata parameters did the user include?

  if (l == 0) {
    genedf <- rowMeans(dataset)
    names(genedf) <- rownames(dataset)
    umis <- colSums(dataset)
    return(list(geneMeans = matrix(genedf, nrow = length(genedf)), umis = matrix(umis, nrow = length(umis))))
  } else {
    metadata_df <- data.frame(mget(metadata))

    tbl <- as.data.frame(data.table(metadata_df)[, .N, by = metadata])

    if (l == 1) {
      metadata_df_concat <- metadata_df[,1]
      tbl_concat <- tbl[,1]
    } else if (l == 2) {
      metadata_df_concat <- paste(metadata_df[,1], metadata_df[,2], sep = "")
      tbl_concat <- paste(tbl[,1], tbl[,2], sep = "")
    } else {
      metadata_df_concat <- paste(metadata_df[,1], metadata_df[,2], metadata_df[,3], sep = "")
      tbl_concat <- paste(tbl[,1], tbl[,2], tbl[,3], sep = "")
    }

    genedf <- matrix(NA, ncol = nrow(tbl), nrow = eval(nrow(dataset) + l))
    umis <- matrix(NA, ncol = nrow(tbl), nrow = eval(l + max(tbl[,eval(l+1)])))

    genedf[1:l,] <- t(tbl[,1:l])

    umis[1:l,] <- t(tbl[,1:l])

    for (i in 1:nrow(tbl)) {
      genedf[eval(l+1):eval(nrow(dataset) + l),i] <- rowMeans(dataset[,which(metadata_df_concat == tbl_concat[i])])
      umiVector <- colSums(dataset[,which(metadata_df_concat == tbl_concat[i])])
      umis[eval(l+1):eval(l+length(umiVector)),i] <- umiVector
    }

    rownames(genedf) <- c(colnames(metadata_df), rownames(dataset))
    rownames(umis) <- c(colnames(metadata_df), as.character(1:eval(nrow(umis)-2)))

    return(list(geneMeans = genedf, umis = umis))
  }
}
