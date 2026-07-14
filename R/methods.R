################################################################################
# BiSer0608: algorithms used in the submitted benchmark
################################################################################

library(TSP)
library(seriation)
library(scran)
library(igraph)
library(NMF)

# The submitted main comparison contains exactly seven methods.
submission_methods <- c(
  "biser", "bs", "tsp_seri", "spec_seri", "Heatmap", "MESBC", "NMF"
)
methname <- submission_methods

# Retained only for optional legacy analyses; these methods are not part of the
# seven-method BiSer0608 benchmark.
Rserimeth <- c("BEA_TSP", "PCA_angle", "Heatmap", "PCA")
bicmeth <- c("bccc", "plaid", "MESBC")

ensure_dimnames <- function(mat) {
  if (is.null(rownames(mat))) rownames(mat) <- paste0("r", seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("c", seq_len(ncol(mat)))
  mat
}

prepare_bipartite_matrix <- function(mat) {
  mat <- ensure_dimnames(as.matrix(mat))
  if (any(!is.finite(mat))) stop("The input matrix contains non-finite values.")

  w <- mat - min(0, min(mat))
  keep_row <- rowSums(w) > 0
  keep_col <- colSums(w) > 0
  if (!any(keep_row) || !any(keep_col)) {
    stop("No positive-degree row or column remains after preprocessing.")
  }

  # Removing zero-degree nodes can expose additional zero-degree nodes on the
  # opposite side. Iterate until the active bipartite graph is stable.
  repeat {
    w_active <- w[keep_row, keep_col, drop = FALSE]
    next_row <- keep_row
    next_col <- keep_col
    next_row[keep_row] <- rowSums(w_active) > 0
    next_col[keep_col] <- colSums(w_active) > 0
    if (identical(next_row, keep_row) && identical(next_col, keep_col)) break
    keep_row <- next_row
    keep_col <- next_col
  }

  list(
    mat = mat[keep_row, keep_col, drop = FALSE],
    w = w[keep_row, keep_col, drop = FALSE],
    removed_rows = rownames(mat)[!keep_row],
    removed_cols = colnames(mat)[!keep_col]
  )
}

safe_cor_similarity <- function(x) {
  sim <- suppressWarnings(cor(t(x), use = "pairwise.complete.obs"))
  if (is.null(dim(sim))) sim <- matrix(sim, nrow = nrow(x), ncol = nrow(x))
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  sim
}

# ==============================================================================
# BiSer: normalized bipartite SVD, shared embedding, and joint TSP path
# ==============================================================================
biser <- function(mat, simmeth = "t", noise = TRUE, pct = 0.3,
                  sparsify = NULL, n_starts = 100L) {
  mat <- ensure_dimnames(as.matrix(mat))
  prep <- prepare_bipartite_matrix(mat)
  mat_active <- prep$mat
  w <- prep$w
  m <- nrow(w)
  p <- ncol(w)

  if (length(prep$removed_rows) + length(prep$removed_cols) > 0L) {
    warning("All-zero rows/columns were removed before BiSer analysis.")
  }

  d1 <- rowSums(w)
  d2 <- colSums(w)
  D1 <- diag(1 / sqrt(d1), nrow = m)
  D2 <- diag(1 / sqrt(d2), nrow = p)
  w_standard <- D1 %*% w %*% D2
  fit <- svd(w_standard)

  # The first normalized singular component is the degree-related trivial
  # component. BiSer0608 uses every available component after that component.
  keep <- if (length(fit$d) > 1L) 2:length(fit$d) else 1L
  U <- fit$u[, keep, drop = FALSE]
  V <- fit$v[, keep, drop = FALSE]
  lambda <- fit$d[keep]
  Y <- rbind(D1 %*% U, D2 %*% V) %*% diag(lambda, nrow = length(lambda))
  joint_names <- c(rownames(mat_active), colnames(mat_active))
  rownames(Y) <- joint_names

  if (identical(simmeth, "t")) {
    # Z has one row per graph node; the dimensionally valid n x n Gram matrix
    # corresponding to the manuscript definition is K = Z Z^T.
    sim <- tcrossprod(Y)
  } else if (identical(simmeth, "cor")) {
    sim <- safe_cor_similarity(Y)
  } else {
    stop("simmeth must be either 'cor' or 't'.")
  }
  rownames(sim) <- colnames(sim) <- joint_names

  # `noise` is retained for backward compatibility with the original scripts.
  # The submitted simulations use noise=TRUE and therefore no sparsification.
  if (is.null(sparsify)) sparsify <- !isTRUE(noise)
  if (isTRUE(sparsify)) {
    nz <- which(sim != 0)
    k <- ceiling(length(nz) * pct)
    if (k > 0L) {
      threshold <- sort(sim[nz], partial = k)[k]
      sim[nz[sim[nz] <= threshold]] <- 0
      diag(sim) <- 1
    }
  }

  tsp_dist <- as.dist(max(sim) - sim)
  tsp_problem <- insert_dummy(TSP(tsp_dist), label = "cut_here")
  tsp_tour <- solve_TSP(
    tsp_problem,
    method = "repetitive_nn",
    control = list(rep = as.integer(n_starts), two_opt = TRUE)
  )
  tour_order <- cut_tour(tsp_tour, cut = "cut_here", exclude_cut = TRUE)
  joint_order <- joint_names[as.integer(tour_order)]

  row_names <- joint_order[joint_order %in% rownames(mat_active)]
  col_names <- joint_order[joint_order %in% colnames(mat_active)]
  reordered_mat <- mat_active[row_names, col_names, drop = FALSE]

  list(
    reordered_mat = reordered_mat,
    Y = Y,
    sim = sim,
    joint_order = joint_order,
    row_order = match(row_names, rownames(mat)),
    col_order = match(col_names, colnames(mat)),
    retained_components = keep,
    removed_rows = prep$removed_rows,
    removed_cols = prep$removed_cols
  )
}

# ==============================================================================
# Independent seriation baselines and ablations
# ==============================================================================
spec_seri <- function(mat) {
  mat <- ensure_dimnames(as.matrix(mat))
  list(
    row_order = get_order(seriate(dist(mat), method = "spectral")),
    col_order = get_order(seriate(dist(t(mat)), method = "spectral"))
  )
}

run_tsp_seriation <- function(mat) {
  mat <- ensure_dimnames(as.matrix(mat))
  row_cor <- suppressWarnings(cor(t(mat), use = "pairwise.complete.obs"))
  col_cor <- suppressWarnings(cor(mat, use = "pairwise.complete.obs"))
  row_cor[!is.finite(row_cor)] <- 0
  col_cor[!is.finite(col_cor)] <- 0
  diag(row_cor) <- diag(col_cor) <- 1

  row_tour <- solve_TSP(TSP(as.dist(1 - row_cor)), method = "nearest_insertion")
  col_tour <- solve_TSP(TSP(as.dist(1 - col_cor)), method = "nearest_insertion")
  list(row_order = as.integer(row_tour), col_order = as.integer(col_tour))
}

Rseriation <- function(mat, method) {
  mat <- ensure_dimnames(as.matrix(mat))

  # "Heatmap" in the manuscript denotes optimal leaf ordering (OLO), applied
  # independently to row and column Euclidean-distance matrices.
  if (method %in% c("Heatmap", "OLO")) {
    return(list(
      row_order = get_order(seriate(dist(mat), method = "OLO")),
      col_order = get_order(seriate(dist(t(mat)), method = "OLO"))
    ))
  }

  out <- seriate(mat - min(mat), method = method)
  reordered <- seriation::permute(mat, out)
  list(
    row_order = match(rownames(reordered), rownames(mat)),
    col_order = match(colnames(reordered), colnames(mat))
  )
}

Yang_bs <- function(mat, method = "bs") {
  if (!identical(method, "bs")) stop("BiSer_SVD is implemented as method='bs'.")
  mat <- ensure_dimnames(as.matrix(mat))
  prep <- prepare_bipartite_matrix(mat)
  W <- prep$w
  d1 <- rowSums(W)
  d2 <- colSums(W)
  W_tilde <- W / sqrt(d1)
  W_tilde <- t(t(W_tilde) / sqrt(d2))
  fit <- svd(W_tilde)

  component <- if (length(fit$d) > 1L) 2L else 1L
  u <- fit$u[, component] / sqrt(d1)
  v <- fit$v[, component] / sqrt(d2)
  row_names <- c(rownames(W)[order(u)], prep$removed_rows)
  col_names <- c(colnames(W)[order(v)], prep$removed_cols)
  list(
    row_order = match(row_names, rownames(mat)),
    col_order = match(col_names, colnames(mat)),
    component = component
  )
}

# ==============================================================================
# Biclustering baselines
# ==============================================================================
bicluster <- function(mat, method, K = 2:10) {
  mat <- ensure_dimnames(as.matrix(mat))

  if (!identical(method, "MESBC")) {
    if (!requireNamespace("biclust", quietly = TRUE)) {
      stop("The optional 'biclust' package is required for legacy method ", method, ".")
    }
    bc_res <- NULL
    if (identical(method, "bcspectral")) {
      bc_res <- biclust::biclust(mat, method = biclust::BCSpectral())
    }
    if (identical(method, "bccc")) {
      bc_res <- biclust::biclust(mat, method = biclust::BCCC())
    }
    if (identical(method, "plaid")) {
      bc_res <- biclust::biclust(
        mat, method = biclust::BCPlaid(), cluster = "b", fit.model = ~m + a + b
      )
    }
    if (is.null(bc_res)) stop("Unknown biclustering method: ", method)

    row_labels <- apply(bc_res@RowxNumber, 1, which.max)
    col_labels <- apply(bc_res@NumberxCol, 2, which.max)
    return(list(
      row_order = order(row_labels), col_order = order(col_labels),
      row_label = row_labels, col_label = col_labels, pred_k = bc_res@Number
    ))
  }

  prep <- prepare_bipartite_matrix(mat)
  w <- prep$w
  if (length(prep$removed_rows) + length(prep$removed_cols) > 0L) {
    stop("MESBC requires positive-degree rows and columns in the benchmark matrix.")
  }
  m <- nrow(w)
  p <- ncol(w)
  D1 <- diag(1 / sqrt(rowSums(w)), nrow = m)
  D2 <- diag(1 / sqrt(colSums(w)), nrow = p)
  fit <- svd(D1 %*% w %*% D2)

  clust <- matrix(NA_integer_, nrow = m + p, ncol = length(K),
                  dimnames = list(NULL, as.character(K)))
  for (k in K) {
    r <- min(k, length(fit$d))
    Y <- rbind(
      D1 %*% fit$u[, seq_len(r), drop = FALSE],
      D2 %*% fit$v[, seq_len(r), drop = FALSE]
    ) %*% diag(fit$d[seq_len(r)], nrow = r)
    clust[, as.character(k)] <- kmeans(
      Y, centers = k, iter.max = 100, nstart = 100
    )$cluster
  }

  Y_full <- rbind(D1 %*% fit$u, D2 %*% fit$v) %*%
    diag(fit$d, nrow = length(fit$d))
  graph <- buildSNNGraph(t(Y_full))
  modularity_by_k <- vapply(K, function(k) {
    modularity(graph, clust[, as.character(k)])
  }, numeric(1))
  best_k <- K[which.max(modularity_by_k)]
  labels <- clust[, as.character(best_k)]
  row_labels <- labels[seq_len(m)]
  col_labels <- labels[m + seq_len(p)]

  list(
    row_order = order(row_labels), col_order = order(col_labels),
    row_label = row_labels, col_label = col_labels,
    pred_k = best_k, modularity = modularity_by_k, Y = Y_full
  )
}

run_nmf <- function(mat, K = 2:10, nrun = 10L, seed = 1L) {
  mat <- ensure_dimnames(as.matrix(mat))
  mat_nn <- mat - min(mat) + .Machine$double.eps

  # Use the same joint SNN graph and modularity selection rule as described for
  # the submitted benchmark, evaluating every K from 2 through 10.
  prep <- prepare_bipartite_matrix(mat_nn)
  w <- prep$w
  D1 <- diag(1 / sqrt(rowSums(w)), nrow = nrow(w))
  D2 <- diag(1 / sqrt(colSums(w)), nrow = ncol(w))
  fit <- svd(D1 %*% w %*% D2)
  Y <- rbind(D1 %*% fit$u, D2 %*% fit$v) %*%
    diag(fit$d, nrow = length(fit$d))
  graph <- buildSNNGraph(t(Y))

  fits <- list()
  modularity_by_k <- setNames(rep(NA_real_, length(K)), as.character(K))
  for (k in K) {
    result <- try(
      nmf(mat_nn, rank = k, method = "brunet", nrun = as.integer(nrun),
          seed = seed + k, .options = ""),
      silent = TRUE
    )
    if (inherits(result, "try-error")) next
    fits[[as.character(k)]] <- result
    row_labels <- predict(result, "rows")
    col_labels <- predict(result, "columns")
    modularity_by_k[as.character(k)] <- modularity(
      graph, c(row_labels, col_labels)
    )
  }

  if (all(is.na(modularity_by_k))) {
    return(list(
      row_order = seq_len(nrow(mat)), col_order = seq_len(ncol(mat)),
      row_label = rep(1L, nrow(mat)), col_label = rep(1L, ncol(mat)),
      row_cluster = rep(1L, nrow(mat)), col_cluster = rep(1L, ncol(mat)),
      success = FALSE
    ))
  }

  best_k <- as.integer(names(which.max(modularity_by_k)))
  result <- fits[[as.character(best_k)]]
  row_labels <- as.numeric(predict(result, "rows"))
  col_labels <- as.numeric(predict(result, "columns"))
  list(
    row_order = order(row_labels), col_order = order(col_labels),
    row_label = row_labels, col_label = col_labels,
    row_cluster = row_labels, col_cluster = col_labels,
    pred_k = best_k, modularity = modularity_by_k, success = TRUE
  )
}
