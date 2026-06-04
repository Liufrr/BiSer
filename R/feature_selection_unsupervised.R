################################################################################
# Unsupervised Feature Selection for Gene Expression Data
#
# Standard bioinformatics approaches:
#   - hvg:      scran::modelGeneVar — standard scRNA-seq HVG selection,
#               decomposes variance into biological/technical components
#   - variance: simple variance ranking — standard for microarray/bulk data
#   - mad:      Median Absolute Deviation — robust alternative
#
# Dependencies: Matrix, scran, SingleCellExperiment, scuttle
################################################################################

library(Matrix)

#' Select top variable genes (unsupervised, standard bioinformatics methods)
#'
#' @param data Gene expression matrix (rows = genes, cols = cells/samples)
#' @param n_top Number of top genes to retain (default: 500)
#' @param method "hvg" (scran), "variance", or "mad"
#' @param norm Whether to apply log2(x + 1) normalization (default: TRUE)
#'
#' @return List with scores, data_log, top_genes, top_mat
get_clustering_matrix <- function(data, n_top = 500,
                                   method = "hvg", norm = TRUE) {

  message(">>> Initial dimensions: ", nrow(data), " genes x ", ncol(data),
      " cells/samples")

  # --- Step 1: Filter low-expression genes ---
  n_cells   <- ncol(data)
  min_cells <- max(3, n_cells * 0.01)
  if (inherits(data, "sparseMatrix")) {
    n_expressed <- Matrix::rowSums(data > 0)
  } else {
    n_expressed <- rowSums(data > 0)
  }
  data <- data[n_expressed >= min_cells, ]
  message(">>> After low-expression filtering: ", nrow(data), " genes")

  # --- Step 2: Log-normalization (optional) ---
  data_log <- if (norm) log2(data + 1) else data
  if (inherits(data_log, "sparseMatrix")) data_log <- as.matrix(data_log)

  # --- Step 3: Feature selection ---
  message(">>> Feature selection method: ", method)

  if (method == "hvg") {
    # ---- scran::modelGeneVar (standard scRNA-seq HVG) ----
    library(SingleCellExperiment)
    library(scuttle)
    library(scran)
    sce <- SingleCellExperiment(list(logcounts = as.matrix(data_log)))
    dec <- modelGeneVar(sce)
    scores <- dec$bio
    names(scores) <- rownames(dec)
    scores[is.na(scores)] <- 0

  } else if (method == "variance") {
    # ---- Simple variance (standard for microarray) ----
    scores <- apply(data_log, 1, var)
    scores[is.na(scores)] <- 0
    names(scores) <- rownames(data_log)

  } else if (method == "mad") {
    # ---- MAD (robust alternative) ----
    scores <- apply(data_log, 1, mad)
    scores[is.na(scores)] <- 0
    names(scores) <- rownames(data_log)

  } else {
    stop("Unknown method '", method, "'. Choose 'hvg', 'variance', or 'mad'.")
  }

  # --- Step 4: Select top N genes ---
  n_select  <- min(n_top, length(scores))
  top_genes <- names(sort(scores, decreasing = TRUE))[1:n_select]
  top_mat   <- data_log[top_genes, ]

  message(">>> Selected top ", n_select, " genes")

  return(list(
    scores    = scores,
    data_log  = data_log,
    top_genes = top_genes,
    top_mat   = top_mat
  ))
}
