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
