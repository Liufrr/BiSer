################################################################################
# BiSer Benchmarking: Feature Selection for scRNA-seq Data
#
# Gene filtering and feature selection for single-cell RNA-seq data.
# Methods: MAD (unsupervised), Variance (unsupervised), ANOVA (supervised)
#
# Dependencies: Matrix (for sparse matrix support)
################################################################################

library(Matrix)

#' Select top genes from an expression matrix for downstream analysis
#'
#' @param data Gene expression matrix (rows = genes, cols = cells).
#'             Accepts matrix or sparseMatrix.
#' @param label Cell type labels (required for ANOVA method; NULL for unsupervised)
#' @param n_top Number of top genes to retain (default: 500)
#' @param method Feature selection method: "mad", "variance", or "anova"
#' @param norm Whether to apply log2(x + 1) normalization (default: TRUE)
#'
#' @return List with: scores (gene-level scores), data_log (normalized matrix)
get_clustering_matrix <- function(data, label = NULL, n_top = 500, norm = TRUE) {
  if (is.null(label)) stop("label is required for ANOVA-based feature selection")
  message(">>> Initial dimensions: ", dim(data)[1], " genes x ", dim(data)[2], " cells")

  # Step 1: Pre-filter low-expression genes (expressed in < 1% of cells)
  n_cells   <- ncol(data)
  min_cells <- max(3, n_cells * 0.01)

  if (inherits(data, "sparseMatrix")) {
    n_expressed <- Matrix::rowSums(data > 0)
  } else {
    n_expressed <- rowSums(data > 0)
  }

  keep_genes <- n_expressed >= min_cells
  data <- data[keep_genes, ]
  message(">>> After filtering: ", dim(data)[1], " genes")

  # Step 2: Log-normalization (optional)
  data_log <- if (norm) log2(data + 1) else data

  # Step 3: Feature scoring
  message(">>> Feature selection")

  pvals <- apply(data_log, 1, function(x) {
    tryCatch(oneway.test(x ~ label, var.equal = FALSE)$p.value,
             error = function(e) 1)
  })
  scores <- -log10(pvals + 1e-300)
  scores[is.na(scores)] <- 0

  return(list(scores = scores, data_log = data_log))
}
