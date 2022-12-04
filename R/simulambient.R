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
#' @param cellTypes A vector of length k that provides the cell type of each nucleus in \code{dataset}. If all nuclei are from one type or cell type is not to be accounted for, parameter can be left blank.
#'
#' @param batch A vector of length k that provides the batch of the nucleus in \code{dataset}. If there is no batch effect, parameter can be left blank. This parameter must be provided if \code{bgBatch} is provided (as required by DecontX).
#'
#' @param bgBatch A vector of length j corresponding to batch information of the ambient mRNA droplets in \code{contamination}. This parameter must be provided if \code{batch} is provided (as required by DecontX).
#'
#' @param sex A vector of length k that provides the associated sex of each droplet in \code{dataset}.
#'
#' @param seed Seed for random number generation. Can be set by the user. Defaults to 615.
#'
#' @param ... Any extra parameters passed to DecontX from the \code{celda} package.
#'
#' @export
#'
#' @import stats


simulambient <- function(dataset, contamination_levels = 1, contamination = NULL, cellTypes = NULL, batch = NULL, bgBatch = NULL, sex = NULL, seed = 615, ...) {

  ## STEPS 0A-0C: CHECK ARGUMENTS FOR VALIDITY & SAVE UMI INFORMATION

  print("Checking arguments for validity...")

  metadata <- check_args(dataset = dataset, contamination = contamination, cellTypes = cellTypes, batch = batch, bgBatch = bgBatch, sex = sex)
  preDecontX_umis <- get_umis(dataset = dataset, metadata = metadata, cellTypes = cellTypes, sex = sex)

  print("Removing genes with no expression...")

  expressed_genes <- which(rowSums(dataset) > 0)
  dataset <- dataset[expressed_genes,]
  if (!is.null(contamination)) {
    contamination <- contamination[expressed_genes,]
  }

  ## STEP 1: DECONTAMINATE

  print("Decontaminating snRNA-seq data with DecontX...")

  decontaminated_counts <- decontaminate(dataset, contamination, cellTypes, batch, ...)

  return(decontaminated_counts)

  ## STEP 2: ESTIMATE PARAMETERS

  print("Estimating snRNA-seq data gene expression profile...")

  params <- estimate_params(dataset = decontaminated_counts, metadata = metadata, cellTypes = cellTypes, sex = sex)

  ## STEP 3: SIMULATE GROUND TRUTH

  print("Simulating ground truth...")

  umis <- params$umis
  num_nuclei_per_metadata <- apply(umis, MARGIN = 2, FUN = function(x) length(stats::na.omit(x)))

  counts <- simulate_counts(geneMeans = params$geneMeans, geneSDs = params$geneSDs, numNucleiPerGroup = num_nuclei_per_metadata, seed = seed)

  ## STEP 4: SCALE GROUND TRUTH

  print("Applying nuclei-specific scale factors...")

  umis_vector <- na.omit(as.vector(umis))
  scaling_factor <- umis_vector / colSums(counts)
  scaling_factor[!is.finite(scaling_factor)] <- 0

  counts <- sweep(counts, 2, scaling_factor, FUN = "/") # Ground truth counts
  counts <- round(counts)

  output_list <- list(counts)

  ## STEP 5: GET AMBIENT MRNA PROFILE MEANS

  print("Estimating ambient mRNA droplets gene expression profile...")

  if (!is.null(bgBatch)) {
    metadata2 <- "cellTypes" # using this parameter as batch substitute
  } else {
    metadata2 <- c()
  }
  contamination_params <- estimate_params(dataset = contamination, metadata = metadata2, cellTypes = bgBatch, sex = NULL)

  contamination_umis <- contamination_params$umis
  if (is.null(batch)) {
    num_droplets_per_metadata <- c(ncol(counts))
  } else {
    ## TODO: Fix this
    num_droplets_per_metadata <- as.data.frame(table(batch))[,2]
  }

  print("Simulating contamination counts...")

  contamination_counts <- simulate_counts(geneMeans = contamination_params$geneMeans * 100, geneSDs = contamination_params$geneSDs, numNucleiPerGroup = num_droplets_per_metadata, seed = seed)

  ## STEP 6: ADD DESIRED CONTAMINATION LEVELS TO GROUND TRUTH DATA

  preDecontX_umis_vector <- na.omit(as.vector(preDecontX_umis))

  contamination_to_add_vector <- preDecontX_umis_vector - umis_vector

  for (i in 1:length(contamination_levels)) {
    scaling_factor <- contamination_to_add_vector / colSums(contamination_counts)
    scaling_factor[!is.finite(scaling_factor)] <- 0

    contamination_counts <- sweep(contamination_counts, 2, scaling_factor, FUN = "/") # Contaminated counts
    contamination_counts <- contamination_counts * contamination_levels[i]
    contamination_counts <- round(contamination_counts)
    final_matrix <- counts + contamination_counts
    output_list <- append(output_list, final_matrix)
  }

  names(output_list) <- c("ground_truth", as.character(contamination_levels))
  return(output_list)
}
