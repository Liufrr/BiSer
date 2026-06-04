################################################################################
# BiSer Benchmarking: Comparison Methods
#
# All seriation and biclustering methods used in the benchmarking study.
# Methods included:
#   Seriation:    BiSer, Spectral, TSP, BEA_TSP, PCA_angle, Heatmap, PCA
#   Bipartite:    Yang BS
#   Biclustering: MESBC, BCCC, Plaid, NMF
#
# Dependencies:
#   TSP, seriation, scran, igraph, NMF, biclust
################################################################################

library(TSP)
library(seriation)
library(scran)
library(igraph)
library(biclust)
library(NMF)

# --- Global Constants ---
Rserimeth <- c("BEA_TSP", "PCA_angle", "Heatmap", "PCA")
bicmeth   <- c("bccc", "plaid", "MESBC")
methname  <- c("biser", "spec_seri", "tsp_seri", Rserimeth, "bs", bicmeth, "NMF")

# ==============================================================================
# BiSer: Bipartite Seriation via SVD + TSP
# ==============================================================================
biser <- function(mat, simmeth = "cor", noise = FALSE, pct = 0.3) {
  m <- nrow(mat)
  p <- ncol(mat)
  n <- m + p

  # Step 1: Normalized bipartite embedding
  w <- as.matrix(mat + abs(min(0, range(mat)[1])))
  d1 <- apply(w, 1, sum)
  d2 <- apply(w, 2, sum)
  D1 <- diag(1 / d1^0.5)
  D2 <- diag(1 / d2^0.5)
  w.standard <- D1 %*% w %*% D2
  mysvd <- svd(w.standard)
  u <- mysvd$u; v <- mysvd$v; lambda <- mysvd$d
  Y <- rbind(D1 %*% u, D2 %*% v) %*% diag(lambda)

  # Step 2: Similarity matrix
  if (simmeth == "t")   mat2 <- Y %*% t(Y)
  if (simmeth == "cor") mat2 <- cor(t(Y))
  rownames(mat2) <- colnames(mat2) <- c(rownames(mat), colnames(mat))

  # Step 3: Optional sparsification
  if (!noise) {
    sparsify_global <- function(mat, pct = 0.2) {
      res <- mat
      nz_idx  <- which(res != 0)
      nz_vals <- res[nz_idx]
      k <- ceiling(length(nz_vals) * pct)
      if (k > 0) {
        th <- sort(nz_vals)[k]
        to_zero <- nz_idx[res[nz_idx] <= th]
        res[to_zero] <- 0
      }
      return(res)
    }
    mat2 <- sparsify_global(mat2, pct = 0.3)
  }

  # Step 4: TSP-based reordering
  row_sim  <- mat2
  row_dist <- as.dist(max(row_sim) - row_sim)
  row_tsp  <- insert_dummy(TSP(row_dist), label = "cut_here")
  row_tour <- solve_TSP(row_tsp, method = "repetitive_nn",
                        control = list(rep = 100, two_opt = TRUE))

  # Step 5: Extract row/column ordering
  row_order_native <- cut_tour(row_tour, cut = "cut_here", exclude_cut = TRUE)
  tmp  <- mat2[row_order_native, row_order_native]
  mat3 <- mat[rownames(tmp)[rownames(tmp) %in% rownames(mat)],
              rownames(tmp)[rownames(tmp) %in% colnames(mat)]]

  return(list(
    reordered_mat = mat3,
    Y   = Y,
    sim = mat2,
    row_order = match(rownames(mat3), rownames(mat)),
    col_order = match(colnames(mat3), colnames(mat))
  ))
}

# ==============================================================================
# Spectral Seriation
# ==============================================================================
spec_seri <- function(mat) {
  row_dist  <- dist(mat)
  row_ser   <- seriate(row_dist, method = "spectral")
  row_order <- get_order(row_ser)

  col_dist  <- dist(t(mat))
  col_ser   <- seriate(col_dist, method = "spectral")
  col_order <- get_order(col_ser)

  return(list(row_order = row_order, col_order = col_order))
}

# ==============================================================================
# TSP Seriation (nearest insertion)
# ==============================================================================
run_tsp_seriation <- function(mat) {
  dist_mat <- as.dist(1 - cor(t(mat)))
  tsp      <- TSP(dist_mat)
  tour     <- solve_TSP(tsp, method = "nearest_insertion")
  row_order <- as.integer(tour)

  dist_mat <- as.dist(1 - cor(mat))
  tsp      <- TSP(dist_mat)
  tour     <- solve_TSP(tsp, method = "nearest_insertion")
  col_order <- as.integer(tour)

  return(list(row_order = row_order, col_order = col_order))
}

# ==============================================================================
# R seriation package methods (BEA_TSP, PCA_angle, Heatmap, PCA)
# ==============================================================================
Rseriation <- function(mat, method) {
  out <- seriate(mat - min(mat), method = method)
  mat_reordered <- seriation::permute(mat, out)
  return(list(
    row_order = match(rownames(mat_reordered), rownames(mat)),
    col_order = match(colnames(mat_reordered), colnames(mat))
  ))
}

# ==============================================================================
# Yang Bipartite Spectral (BS)
# ==============================================================================
Yang_bs <- function(mat, method) {

  # BS: second singular vector ordering
  bs <- function(W, k.bs = 2) {
    W <- as.matrix(W)
    n <- nrow(W); m <- ncol(W)
    if (is.null(dimnames(W))) {
      rownames(W) <- paste0("r", 1:n)
      colnames(W) <- paste0("c", 1:m)
    }
    rn <- rownames(W); cn <- colnames(W)
    d1 <- apply(W, 1, sum); d2 <- apply(W, 2, sum)
    rn0 <- rn[d1 == 0]; cn0 <- cn[d2 == 0]

    W1 <- W[d1 != 0, d2 != 0]
    d1 <- d1[d1 != 0]; d2 <- d2[d2 != 0]
    W.tilde <- W1 / sqrt(d1)
    W.tilde <- t(t(W.tilde) / sqrt(d2))
    tmp <- svd(W.tilde)
    U <- tmp$u; V <- tmp$v
    U <- U / sqrt(d1); V <- V / sqrt(d2)
    rownames(U) <- rownames(W1); rownames(V) <- colnames(W1)
    U <- as.matrix(U); V <- as.matrix(V)
    if (ncol(U) == 1 | ncol(V) == 1) k.bs <- 1

    u1 <- U[, k.bs]; v1 <- V[, k.bs]
    order.row <- order(u1); order.col <- order(v1)
    tmp2 <- W1[order.row, order.col]
    rn_out <- c(rownames(tmp2), rn0)
    cn_out <- c(colnames(tmp2), cn0)
    return(list(
      row_order = match(rn_out, rownames(W)),
      col_order = match(cn_out, colnames(W))
    ))
  }

  if (method == "bs")  out <- bs(mat)
  return(out)
}

# ==============================================================================
# Biclustering methods (BCCC, Plaid, MESBC)
# ==============================================================================
bicluster <- function(mat, method, K = 2:10) {

  if (method != "MESBC") {

    if (method == "bcspectral") bc_res <- biclust(mat, method = BCSpectral())
    if (method == "bccc")       bc_res <- biclust(mat, method = BCCC())
    if (method == "plaid")      bc_res <- biclust(mat, method = BCPlaid(),
                                                   cluster = "b",
                                                   fit.model = ~m + a + b)

    row_labels <- apply(bc_res@RowxNumber, 1, which.max)
    col_labels <- apply(bc_res@NumberxCol, 2, which.max)
    row_order  <- order(row_labels)
    col_order  <- order(col_labels)

    return(list(
      row_order = row_order, col_order = col_order,
      row_label = row_labels, col_label = col_labels,
      pred_k = bc_res@Number
    ))
  }

  # MESBC: Multiway Embedding Spectral Biclustering
  MESBC_fn <- function(data, K = K) {
    m <- nrow(data); p <- ncol(data); n <- m + p
    w <- as.matrix(data + abs(min(0, range(data)[1])))
    d1 <- apply(w, 1, sum); d2 <- apply(w, 2, sum)
    D1 <- diag(1 / d1^0.5); D2 <- diag(1 / d2^0.5)
    w.standard <- D1 %*% w %*% D2
    mysvd <- svd(w.standard)

    clust <- matrix(nrow = n, ncol = length(K))
    colnames(clust) <- as.character(K)
    for (k in K) {
      if (k <= p & k <= m) {
        u <- mysvd$u[, 1:k]; v <- mysvd$v[, 1:k]; lambda <- mysvd$d[1:k]
      } else {
        u <- mysvd$u; v <- mysvd$v; lambda <- mysvd$d
      }
      u <- as.matrix(u); v <- as.matrix(v)
      Y <- rbind(D1 %*% u, D2 %*% v) %*% diag(lambda)
      clus.out <- kmeans(Y, centers = k, iter.max = 100, nstart = 100)
      clust[, as.character(k)] <- clus.out$cluster
    }

    if (length(K) == 1) clust <- clust[, 1]

    Y_full <- rbind(D1 %*% mysvd$u, D2 %*% mysvd$v) %*% diag(mysvd$d)
    g <- buildSNNGraph(t(Y_full))
    mod <- numeric()
    for (k_idx in 2:(ncol(clust) - 1))
      mod[k_idx] <- modularity(g, clust[, as.character(K[k_idx])])
    kbest <- as.character(which(mod == max(mod, na.rm = TRUE))[1])
    clust <- clust[, kbest]

    rowlabel <- clust[1:m]
    collabel <- clust[(m + 1):n]

    return(list(
      row_order = order(rowlabel), col_order = order(collabel),
      row_label = rowlabel, col_label = collabel,
      Y = Y_full
    ))
  }

  if (method == "MESBC") return(MESBC_fn(mat, K))
}

# ==============================================================================
# NMF Biclustering (with modularity-based K selection)
# ==============================================================================
run_nmf <- function(mat, K = 2:10) {


  # Build bipartite graph for modularity evaluation
  graph_out <- function(data) {
    m <- nrow(data); p <- ncol(data)
    w <- as.matrix(data + abs(min(0, range(data)[1])))
    d1 <- apply(w, 1, sum); d2 <- apply(w, 2, sum)
    D1 <- diag(1 / d1^0.5); D2 <- diag(1 / d2^0.5)
    w.standard <- D1 %*% w %*% D2
    mysvd <- svd(w.standard)
    Y <- rbind(D1 %*% mysvd$u, D2 %*% mysvd$v) %*% diag(mysvd$d)
    g <- buildSNNGraph(t(Y))
    return(g)
  }

  mat_nn <- mat - min(mat) + 1   # ensure non-negative
  nmfout <- list()
  mod <- numeric()
  g <- graph_out(mat)

  for (k in K) {
    try({
      nmfout[[k]] <- nmf(mat_nn, rank = k)
    }, silent = TRUE)
    if (k > length(nmfout)) next
    row_cluster <- predict(nmfout[[k]], "rows")
    col_cluster <- predict(nmfout[[k]], "columns")
    mod[k] <- modularity(g, c(row_cluster, col_cluster))
  }

  # Retry once if all attempts failed
  if (length(mod) == 0) {
    for (k in K) {
      try({
        nmfout[[k]] <- nmf(mat_nn, rank = k)
      }, silent = TRUE)
      if (k > length(nmfout)) next
      row_cluster <- predict(nmfout[[k]], "rows")
      col_cluster <- predict(nmfout[[k]], "columns")
      mod[k] <- modularity(g, c(row_cluster, col_cluster))
    }
  }

  if (length(mod) == 0) {
    return(list(
      row_order = 1:nrow(mat), col_order = 1:ncol(mat),
      row_label = rep(1, nrow(mat)), col_label = rep(1, ncol(mat)),
      success = FALSE
    ))
  }

  res <- nmfout[[which.max(mod)]]
  row_cluster <- predict(res, "rows")
  col_cluster <- predict(res, "columns")

  return(list(
    row_order   = order(row_cluster),
    col_order   = order(col_cluster),
    row_cluster = as.numeric(row_cluster),
    col_cluster = as.numeric(col_cluster),
    success     = TRUE
  ))
}
