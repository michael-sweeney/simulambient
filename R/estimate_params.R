#' Estimate Parameters from Dataset
#'
#' estimate_params will take a snRNA-seq dataset as input and estimate parameters for downstream analysis.
#'
#' @param dataset An n x k matrix-like object where n = number of genes and k = number of nuclei. Assumed to be demultiplexed and cell type-annotated.
#'
#' @param metadata A vector of strings corresponding to the metadata made available by the user. This argument will be handled inside the main function \code{simulambient} and is the return value of \code{check_args}.
#'
#' @param cellTypes A vector of length k that provides the cell type of each nucleus in \code{dataset}. If all nuclei are from one type or cell type is not to be accounted for, parameter can be left blank.
#'
#' @param batch A vector of length k that provides the batch of the nucleus in \code{dataset}. If there is no batch effect, parameter can be left blank.
#'
#' @param sample A vector of length k that provides the associated sex of each droplet in \code{dataset}.
#'
#' @return A list containing two values. The first value is a matrix \code{geneMeans} that lists the mean counts value for each gene (for that particular metadata, if applicable). The second value is a matrix \code{umis} that lists the UMI counts for each droplet (for that particular metadata, if applicable).
#'
#' @export
#'
#' @import Matrix
#' @import data.table
#' @import methods

estimate_params <- function(dataset, metadata = NULL, cellTypes = NULL, batch = NULL, sample = NULL) {

  dataset <- methods::as(dataset, "CsparseMatrix") # Convert data to sparse format (although it already should be... might take this out)

  num_genes <- nrow(dataset) # How many genes did the user provide in their dataset?
  num_nuclei <- ncol(dataset) # How many nuclei did the user provide in their dataset?

  l <- length(metadata) # How many metadata parameters did the user include?

  if (l == 0) { # Metadata-free case. Get rowMeans and UMIs per droplet and return.
    genedf <- rowMeans(dataset)
    names(genedf) <- rownames(dataset)
    umis <- colSums(dataset)
    return(list(geneMeans = matrix(genedf, nrow = length(genedf)), umis = matrix(umis, nrow = length(umis))))
  } else {

    metadata_df <- data.frame(mget(metadata)) # Make dataframe of metadata objects

    tbl <- as.data.frame(data.table::data.table(metadata_df)[, .N, by = metadata]) # Going to iterate over all unique permutations of metadata combinations

    # Point of this if/else switch is to find droplets in dataset corresponding to particular metadata permutation
    if (l == 1) { # One of cellTypes, batch, sample provided
      metadata_df_concat <- metadata_df[,1]
      tbl_concat <- tbl[,1]
    } else if (l == 2) { # Two provided
      metadata_df_concat <- paste(metadata_df[,1], metadata_df[,2], sep = "")
      tbl_concat <- paste(tbl[,1], tbl[,2], sep = "")
    } else {# Three provided
      metadata_df_concat <- paste(metadata_df[,1], metadata_df[,2], metadata_df[,3], sep = "")
      tbl_concat <- paste(tbl[,1], tbl[,2], tbl[,3], sep = "")
    }

    # Initialize matrices for function result
    genedf <- matrix(NA, ncol = nrow(tbl), nrow = eval(nrow(dataset) + l))
    umis <- matrix(NA, ncol = nrow(tbl), nrow = eval(l + max(tbl[,eval(l+1)])))

    genedf[1:l,] <- t(tbl[,1:l]) # Fill genedf with metadata information

    umis[1:l,] <- t(tbl[,1:l]) # Finll umis with metadata information

    # Fill genedf and umis with gene means and UMI information
    for (i in 1:nrow(tbl)) {
      genedf[eval(l+1):eval(nrow(dataset) + l),i] <- rowMeans(dataset[,which(metadata_df_concat == tbl_concat[i])])
      umiVector <- colSums(dataset[,which(metadata_df_concat == tbl_concat[i])])
      umis[eval(l+1):eval(l+length(umiVector)),i] <- umiVector
    }

    # Give metadata rownames
    rownames(genedf) <- c(colnames(metadata_df), rownames(dataset))
    rownames(umis) <- c(colnames(metadata_df), as.character(1:eval(nrow(umis)-2)))

    return(list(geneMeans = genedf, umis = umis))
  }
}
