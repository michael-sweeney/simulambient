#' Calculate Pre-DecontX Counts
#'
#' get_umis will get the pre-DecontX counts for each nucleus so that the scaling factor for Simulambient can be computed.
#'
#' @param dataset An n x k matrix-like object where n = number of genes and k = number of nuclei. Assumed to be demultiplexed and cell type-annotated.
#'
#' @param metadata A vector of strings corresponding to the metadata made available by the user. This argument will be handled inside the main function \code{simulambient} and is the return value of \code{check_args}.
#'
#' @param cellTypes A vector of length k that provides the cell type of each nucleus in \code{dataset}. If all nuclei are from one type or cell type is not to be accounted for, parameter can be left blank.
#'
#' @param sex A vector of length k that provides the associated sex of each droplet in \code{dataset}.
#'
#' @return A list containing four values. The first value is a matrix \code{geneMeans} that lists the mean counts value for each gene (for that particular metadata, if applicable). The second value is a matrix \code{geneMeans} that lists the standad deviation for each gene (for that particular metadata, if applicable). The third value is a matrix \code{umis} that lists the UMI counts for each droplet (for that particular metadata, if applicable). The fourth value is a matrix-like object that holds the metadata correpsonding to the columns/rows of the other returned objects.
#'
#' @export
#'
#' @import Matrix
#' @import matrixStats
#' @import data.table
#' @import methods

get_umis <- function(dataset, metadata = NULL, cellTypes = NULL, sex = NULL) {
  dataset <- methods::as(dataset, "CsparseMatrix") # Convert data to sparse format (although it already should be... might take this out)

  num_genes <- nrow(dataset) # How many genes did the user provide in their dataset?
  num_nuclei <- ncol(dataset) # How many nuclei did the user provide in their dataset?

  l <- length(metadata) # How many metadata parameters did the user include?

  if (l == 0) { # Metadata-free case. Get rowMeans and UMIs per droplet and return.
    umis <- colSums(dataset)
    return(matrix(umis, nrow = length(umis)))
  } else {

    metadata_df <- data.frame(mget(metadata)) # Make dataframe of metadata objects

    tbl <- as.data.frame(data.table::data.table(metadata_df)[, .N, by = metadata]) # Going to iterate over all unique permutations of metadata combinations

    # Point of this if/else switch is to find droplets in dataset corresponding to particular metadata permutation
    if (l == 1) { # One of cellTypes,  sex provided
      metadata_df_concat <- metadata_df[,1]
      tbl_concat <- tbl[,1]
    } else if (l == 2) { # Two provided
      metadata_df_concat <- paste(metadata_df[,1], metadata_df[,2], sep = "")
      tbl_concat <- paste(tbl[,1], tbl[,2], sep = "")
    }

    # Initialize matrices for function result
    umis <- matrix(NA, ncol = nrow(tbl), nrow = eval(max(tbl[,eval(l+1)])))

    # Fill genedf and umis with gene means and UMI information
    for (i in 1:nrow(tbl)) {
      metadata_columns <- which(metadata_df_concat == tbl_concat[i])
      if (length(metadata_columns) > 1) {
        umiVector <- colSums(dataset[,metadata_columns])
      } else {
        umiVector <- sum(dataset[,metadata_columns])
      }
      umis[1:length(umiVector),i] <- umiVector
    }

    # Give metadata rownames
    rownames(umis) <- as.character(1:nrow(umis))

    return(umis)
  }
}
