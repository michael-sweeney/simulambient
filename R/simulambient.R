#' Simulambient
#'
#' simulambient is the main function of the package. It calls the helper functions to simulate snRNA-seq counts with varying levels of contamination.
#'
#' @param dataset An n x k matrix-like object where n = number of genes and k = number of nuclei to serve as the basis of the simulation.
#'
#' @param contamination_levels A numeric value or a vector of numeric values indicating the desired levels of contamination that the user would like to simulate. A value of 1 indicates the user wants to simulate the level of contamination detected in the input dataset; a value of 0.5 indicates that the user wants to simulate half of the contamination detected in the original dataset, etc.
#'
#' @param contamination An n x j matrix-like object where n = number of genes and j = number of ambient mRNA droplets. The droplets in this matrix are believed to be by the user to be composed of ambient mRNA. Including droplets beneath a specific UMI cutoff (i.e. 100) is a basic idea to get ambient mRNA droplets.
#'
#' @param metadata A matrix with k rows that corresponds to metadata of the k nuclei provided in \code{dataset}. Please see \code{cellTypes}, \code{batch}, and \code{bgBatch}. Each column should represent a different aspect of metadata, i.e. sex. It should be noted that the user should be rather parsimonious when choosing which metadata to include. Including too many parameters will increase the number of metadata permutations, which will increase runtime and potentially "overfit" the simulations. Cell types is strongly recommended, anything else is up to the user's discretion.
#'
#' @param cellTypes A vector of length k should cell type information not be provided in \code{metadata}. If it is provided, this parameter should be a string of the column name that represents the cell type information in \code{metadata}. If cell type information is not desired to be included, default is NULL.
#'
#' @param batch A vector of length k that provides the batch of the nucleus in \code{dataset}. If there is no batch effect, parameter can be left blank. This parameter must be provided if \code{bgBatch} is provided (as required by DecontX). Please do not include batch information in the \code{metadata} parameter.
#'
#' @param bgBatch A vector of length j corresponding to batch information of the ambient mRNA droplets in \code{contamination}. This parameter must be provided if \code{batch} is provided (as required by DecontX).
#'
#' @param seed Seed for random number generation. Can be set by the user. Defaults to 615.
#'
#' @param ... Any extra parameters passed to DecontX from the \code{celda} package.
#'
#' @return A sparse matrix or a list of sparse matrices containing simulated data. Note that simulated matrices have the same metadata as the input matrix. For example, for an input matrix and metadata with the first column belonging to a female of cell type A, the first column of the output matrix/matrices will belong to a female of cell type A.
#'
#' @export
#'
#' @import stats


simulambient <- function(dataset, contamination_levels = 1, contamination = NULL, metadata = NULL, cellTypes = NULL, batch = NULL, bgBatch = NULL, seed = 615, ...) {

  ## STEPS 0A-0C: CHECK ARGUMENTS FOR VALIDITY, CONFIGURE METADATA GROUPS, REMOVE GENES WITHOUT EXPRESSION

  message("Checking arguments for validity...")

  metadata <- check_args(dataset = dataset, contamination = contamination, metadata = metadata, cellTypes = cellTypes, batch = batch, bgBatch = bgBatch)

  message("Handling metadata permutations...")

  # metadata_indices --> Based on the user-input dataset, which nucleus belongs to which metadata group?
  # metadata_indices_sorted --> Used to sort user-input dataset to organize it so metadata groups are next to each other
  # metadata_indices_original --> If we sort the user-input dataset, we use this to "unsort" it and get back to original order.

  metadata_tbl <- unique(metadata) # Get unique permutations
  metadata_indices <- match(do.call(paste, data.frame(metadata)), do.call(paste, data.frame(metadata_tbl))) # Assign each nucleus to its metadata permutation

  metadata_sort_df <- data.frame(metadata_indices, 1:length(metadata_indices))
  metadata_indices_sorted <- metadata_sort_df[order(metadata_sort_df[,1]),]
  metadata_indices_sorted_df <- data.frame(metadata_indices_sorted, 1:nrow(metadata_indices_sorted))
  metadata_indices_sorted <- metadata_indices_sorted[,2]

  metadata_indices_sorted_df <- metadata_indices_sorted_df[order(metadata_indices_sorted_df[,2]),]
  metadata_indices_original <- metadata_indices_sorted_df[,3]

  num_nuclei_per_metadata_group <- as.data.frame(table(metadata_indices))[,2]

  preDecontX_umis <- colSums(dataset)[metadata_indices_sorted] # NOW IN METADATA ORDER

  message("Removing genes with no expression...")

  expressed_genes <- which(rowSums(dataset) > 0)
  dataset <- dataset[expressed_genes,]
  if (!is.null(contamination)) {
    contamination <- contamination[expressed_genes,]
  }

  ## STEP 1: DECONTAMINATE

  message("Decontaminating snRNA-seq data with DecontX...")

  if (!is.null(cellTypes)) {
    cellTypes <- metadata[,"cellTypes"]
  }

  decontaminated_counts <- decontaminate(dataset, contamination, cellTypes, batch, bgBatch, ...)

  postDecontX_umis <- colSums(decontaminated_counts)[metadata_indices_sorted] # NOW IN METADATA ORDER

  ## STEP 2: ESTIMATE PARAMETERS

  message("Estimating snRNA-seq data gene expression profile...")

  params <- estimate_params(dataset = decontaminated_counts, metadata = metadata, metadata_indices = metadata_indices)

  ## STEP 3: SIMULATE GROUND TRUTH

  message("Simulating ground truth...")

  counts <- simulate_counts(geneMeans = params$geneMeans, geneSDs = params$geneSDs, numNucleiPerGroup = num_nuclei_per_metadata_group, seed = seed)

  rownames(counts) <- rownames(dataset)
  colnames(counts) <- as.character(1:ncol(counts))

  ground_truth_output <- round(counts)
  ground_truth_output <- ground_truth_output[,metadata_indices_original]
  output_list <- list(ground_truth_output)

  ## STEP 4: SCALE GROUND TRUTH

  message("Applying nuclei-specific scale factors...")

  scaling_factor <- postDecontX_umis / colSums(counts)
  scaling_factor[!is.finite(scaling_factor)] <- 0

  counts <- t(t(counts) * scaling_factor)

  ## STEP 5: GET AMBIENT MRNA PROFILE MEANS

  message("Estimating ambient mRNA droplets gene expression profile...")

  if (!is.null(contamination)) {
    contamination_params <- estimate_ambient_params(contamination = contamination, bgBatch = bgBatch)
  } else {
    contamination_params <- estimate_ambient_params(contamination = dataset, bgBatch = bgBatch)
    # If no contamination reference is provided, you're just going to pull random contamination counts proportional to the observed data, which is suboptimal.
  }

  ## STEP 6: ADD DESIRED CONTAMINATION LEVELS TO GROUND TRUTH DATA

  for (i in 1:length(contamination_levels)) {
    message(paste("Simulating contamination counts at factor ", contamination_levels[i], "...", sep = ""))

    contamination_to_add_vector <- preDecontX_umis - postDecontX_umis
    contamination_to_add_vector <- contamination_to_add_vector * contamination_levels[i]

    contamination_counts <- simulate_contamination(geneMeans = contamination_params, numCounts = contamination_to_add_vector, batch = batch, seed = seed)

    # Scale down ground truth counts to account for UMI displacement
    if (contamination_levels[i] > 1) {
      scaling_factor <- (preDecontX_umis - contamination_to_add_vector) / colSums(counts)
      scaling_factor[!is.finite(scaling_factor)] <- 0
      scaling_factor[scaling_factor < 0] <- 0
      counts_temp <- t(t(counts) * scaling_factor)
    }

    final_matrix <- counts_temp + contamination_counts
    final_matrix <- round(final_matrix)
    final_matrix <- final_matrix[,metadata_indices_original]
    rownames(final_matrix) <- rownames(dataset)
    colnames(final_matrix) <- as.character(1:ncol(final_matrix))
    output_list <- append(output_list, final_matrix)
  }

  names(output_list) <- c("ground_truth", as.character(contamination_levels))
  return(output_list)
}
