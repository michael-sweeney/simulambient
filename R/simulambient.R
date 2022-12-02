#' Simulambient
#'
#' simulambient is the main function of the package.
#'
#'
#'
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
