################################################################################
# BiSer Benchmarking: Scalability & Computational Complexity
#
# Four experiments:
#   1. BiSer step-wise timing (Panel A)
#   2. All methods cross-comparison (Panel B)
#   3. MESBC / NMF per-K timing (Panel C)
#   4. Non-square matrix scaling & memory analysis
#
# Required data: None (generates synthetic data)
################################################################################

source("R/methods.R")
source("R/data_generation.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(NMF)
library(scran)
library(igraph)
library(biclust)

# ==============================================================================
# Configuration
# ==============================================================================
output_dir <- "output"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

n_rep <- 3   # repeats per size

# Matrix sizes
sizes      <- c(25, 50, 75, 100, 150, 200, 300, 400, 500)
comp_sizes <- c(25, 50, 75, 100, 150, 200, 300, 400, 500)

# ==============================================================================
# BiSer with per-step timing
# ==============================================================================
biser_timed <- function(mat, simmeth = "t", noise = TRUE, pct = 0.2) {
  m <- nrow(mat); p <- ncol(mat); n <- m + p
  timings <- list()

  # Step 1: Normalization
  t0 <- proc.time()
  w <- as.matrix(mat + abs(min(0, range(mat)[1])))
  d1 <- apply(w, 1, sum); d2 <- apply(w, 2, sum)
  D1 <- diag(1 / d1^0.5); D2 <- diag(1 / d2^0.5)
  w.standard <- D1 %*% w %*% D2
  timings$normalization <- (proc.time() - t0)[["elapsed"]]

  # Step 2: SVD
  t0 <- proc.time()
  mysvd <- svd(w.standard)
  keep <- if (length(mysvd$d) > 1L) 2:length(mysvd$d) else 1L
  u <- mysvd$u[, keep, drop = FALSE]
  v <- mysvd$v[, keep, drop = FALSE]
  lambda <- mysvd$d[keep]
  timings$svd <- (proc.time() - t0)[["elapsed"]]

  # Step 3: Joint embedding
  t0 <- proc.time()
  Y <- rbind(D1 %*% u, D2 %*% v) %*% diag(lambda, nrow = length(lambda))
  timings$embedding <- (proc.time() - t0)[["elapsed"]]

  # Step 4: Similarity matrix
  t0 <- proc.time()
  if (simmeth == "t") mat2 <- Y %*% t(Y)
  if (simmeth == "cor") mat2 <- cor(t(Y))
  rownames(mat2) <- colnames(mat2) <- c(rownames(mat), colnames(mat))
  timings$similarity <- (proc.time() - t0)[["elapsed"]]

  # Step 5: Sparsification
  t0 <- proc.time()
  if (!noise) {
    sparsify_global <- function(mat, pct = 0.2) {
      res <- mat
      nz_idx <- which(res != 0); nz_vals <- res[nz_idx]
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
  timings$sparsification <- (proc.time() - t0)[["elapsed"]]

  # Step 6: TSP solving
  t0 <- proc.time()
  row_sim  <- mat2
  row_dist <- as.dist(max(row_sim) - row_sim)
  row_tsp  <- insert_dummy(TSP(row_dist), label = "cut_here")
  row_tour <- solve_TSP(row_tsp, method = "repetitive_nn",
                        control = list(rep = 100, two_opt = TRUE))
  timings$tsp_solving <- (proc.time() - t0)[["elapsed"]]

  # Step 7: Extraction
  t0 <- proc.time()
  row_order_native <- cut_tour(row_tour, cut = "cut_here", exclude_cut = TRUE)
  tmp <- mat2[row_order_native, row_order_native]
  mat3 <- mat[rownames(tmp)[rownames(tmp) %in% rownames(mat)],
              rownames(tmp)[rownames(tmp) %in% colnames(mat)]]
  timings$extraction <- (proc.time() - t0)[["elapsed"]]

  timings$total <- sum(unlist(timings))
  return(list(
    row_order = match(rownames(mat3), rownames(mat)),
    col_order = match(colnames(mat3), colnames(mat)),
    timings = timings
  ))
}

# ==============================================================================
# MESBC with per-K timing
# ==============================================================================
MESBC_timed <- function(mat, K = 2:10) {
  m <- nrow(mat); p <- ncol(mat); n <- m + p
  per_k_times <- numeric(length(K))
  names(per_k_times) <- as.character(K)

  t_total_start <- proc.time()

  w <- as.matrix(mat + abs(min(0, range(mat)[1])))
  d1 <- apply(w, 1, sum); d2 <- apply(w, 2, sum)
  D1 <- diag(1 / d1^0.5); D2 <- diag(1 / d2^0.5)
  w.standard <- D1 %*% w %*% D2
  mysvd <- svd(w.standard)

  clust <- matrix(nrow = n, ncol = length(K))
  colnames(clust) <- as.character(K)

  for (i in seq_along(K)) {
    k <- K[i]
    t_k_start <- proc.time()
    if (k <= p & k <= m) {
      u <- mysvd$u[, 1:k]; v <- mysvd$v[, 1:k]; lambda <- mysvd$d[1:k]
    } else {
      u <- mysvd$u; v <- mysvd$v; lambda <- mysvd$d
    }
    u <- as.matrix(u); v <- as.matrix(v)
    Y <- rbind(D1 %*% u, D2 %*% v) %*% diag(lambda)
    clus.out <- kmeans(Y, centers = k, iter.max = 100, nstart = 100)
    clust[, as.character(k)] <- clus.out$cluster
    per_k_times[i] <- (proc.time() - t_k_start)[["elapsed"]]
  }

  # Model selection via modularity
  Y_full <- rbind(D1 %*% mysvd$u, D2 %*% mysvd$v) %*% diag(mysvd$d)
  g <- buildSNNGraph(t(Y_full))
  mod <- vapply(K, function(k) {
    modularity(g, clust[, as.character(k)])
  }, numeric(1))
  kbest <- as.character(K[which.max(mod)])
  clust_best <- clust[, kbest]

  rowlabel <- clust_best[1:m]
  collabel <- clust_best[(m + 1):n]
  t_total <- (proc.time() - t_total_start)[["elapsed"]]

  return(list(
    row_order  = order(rowlabel),
    col_order  = order(collabel),
    time_total = t_total,
    time_per_k = per_k_times,
    best_k     = as.integer(kbest)
  ))
}

# ==============================================================================
# NMF with per-K timing
# ==============================================================================
NMF_timed <- function(mat, K = 2:10) {
  mat_nn <- mat - min(mat) + 0.01
  per_k_times <- numeric(length(K))
  names(per_k_times) <- as.character(K)
  cophenetic_vals <- numeric(length(K))

  t_total_start <- proc.time()
  best_score  <- -Inf
  best_result <- NULL
  best_k      <- K[1]

  for (i in seq_along(K)) {
    k <- K[i]
    t_k_start <- proc.time()
    tryCatch({
      res  <- nmf(mat_nn, rank = k, method = "brunet", nrun = 10,
                  seed = 123, .options = "v")
      coph <- cophcor(res)
      cophenetic_vals[i] <- coph
      if (coph > best_score) {
        best_score  <- coph
        best_k      <- k
        best_result <- res
      }
    }, error = function(e) {
      cophenetic_vals[i] <<- NA
    })
    per_k_times[i] <- (proc.time() - t_k_start)[["elapsed"]]
  }

  t_total <- (proc.time() - t_total_start)[["elapsed"]]

  if (!is.null(best_result)) {
    row_labels <- predict(best_result, what = "rows")
    col_labels <- predict(best_result, what = "columns")
    row_order  <- order(row_labels)
    col_order  <- order(col_labels)
  } else {
    row_order <- 1:nrow(mat)
    col_order <- 1:ncol(mat)
  }

  return(list(
    row_order  = row_order,
    col_order  = col_order,
    time_total = t_total,
    time_per_k = per_k_times,
    best_k     = best_k,
    cophenetic = cophenetic_vals
  ))
}

# ==============================================================================
# Safe timing wrapper
# ==============================================================================
safe_time <- function(expr, label = "") {
  t0 <- proc.time()
  result <- tryCatch(eval(expr), error = function(e) {
    cat(sprintf("    [%s] Error: %s\n", label, e$message))
    return(NULL)
  })
  elapsed <- (proc.time() - t0)[["elapsed"]]
  return(list(result = result, elapsed = elapsed))
}

# ==============================================================================
# Experiment 1: BiSer step-wise timing
# ==============================================================================
cat("====== Experiment 1: BiSer step-wise timing ======\n")

stepwise_results <- data.frame()

for (sz in sizes) {
  cat(sprintf("  m=p=%d (n=%d) ...\n", sz, 2 * sz))
  for (rep_i in 1:n_rep) {
    mat <- generate_block_matrix(m = sz, p = sz, seed = 42 + rep_i)
    tryCatch({
      res <- biser_timed(mat, simmeth = "t", noise = TRUE)
      for (step_name in names(res$timings)) {
        stepwise_results <- rbind(stepwise_results, data.frame(
          m = sz, p = sz, n = 2 * sz, rep = rep_i,
          step = step_name, time_sec = res$timings[[step_name]],
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      cat(sprintf("    Error at size %d: %s\n", sz, e$message))
    })
  }
  gc()
}

write.csv(stepwise_results, file.path(output_dir, "stepwise_timings.csv"),
          row.names = FALSE)

# ==============================================================================
# Experiment 2: All methods cross-comparison
# ==============================================================================
cat("====== Experiment 2: All methods comparison ======\n")

comp_results <- data.frame()
perk_results <- data.frame()

all_methods <- list(
  list(name = "BiSer",     type = "seriation"),
  list(name = "Spectral",  type = "seriation"),
  list(name = "TSP",       type = "seriation"),
  list(name = "Heatmap",   type = "Rseriation"),
  list(name = "BS",        type = "seriation"),
  list(name = "MESBC",     type = "biclustering"),
  list(name = "NMF",       type = "biclustering")
)

for (sz in comp_sizes) {
  cat(sprintf("  Size m=p=%d ...\n", sz))

  for (rep_i in 1:n_rep) {
    mat    <- generate_block_matrix(m = sz, p = sz, seed = 42 + rep_i)
    mat_nn <- mat - min(mat) + 0.01  # ensure strictly positive

    for (meth_info in all_methods) {
      meth_name <- meth_info$name
      meth_type <- meth_info$type
      cat(sprintf("    %s ...", meth_name))

      elapsed <- NA
      per_k   <- NULL
      best_k  <- NA

      tryCatch({
        if (meth_name == "BiSer") {
          res     <- biser_timed(mat, simmeth = "t", noise = TRUE)
          elapsed <- res$timings$total

        } else if (meth_name == "Spectral") {
          t0      <- proc.time()
          res     <- spec_seri(mat)
          elapsed <- (proc.time() - t0)[["elapsed"]]

        } else if (meth_name == "TSP") {
          t0      <- proc.time()
          res     <- run_tsp_seriation(mat)
          elapsed <- (proc.time() - t0)[["elapsed"]]

        } else if (meth_type == "Rseriation") {
          t0      <- proc.time()
          res     <- Rseriation(mat, method = meth_name)
          elapsed <- (proc.time() - t0)[["elapsed"]]

        } else if (meth_name == "BS") {
          t0      <- proc.time()
          res     <- Yang_bs(mat_nn, method = "bs")
          elapsed <- (proc.time() - t0)[["elapsed"]]

        } else if (meth_name == "MESBC") {
          res     <- MESBC_timed(mat_nn, K = 2:10)
          elapsed <- res$time_total
          per_k   <- res$time_per_k
          best_k  <- res$best_k

        } else if (meth_name == "NMF") {
          res     <- NMF_timed(mat_nn, K = 2:10)
          elapsed <- res$time_total
          per_k   <- res$time_per_k
          best_k  <- res$best_k
        }

        cat(sprintf(" %.3fs\n", elapsed))

      }, error = function(e) {
        cat(sprintf(" ERROR: %s\n", e$message))
        elapsed <<- NA
      })

      # Record total time
      comp_results <- rbind(comp_results, data.frame(
        m = sz, p = sz, n = 2 * sz, rep = rep_i,
        method      = meth_name,
        method_type = ifelse(meth_type == "biclustering", "biclustering", "seriation"),
        time_sec    = elapsed,
        best_k      = best_k,
        stringsAsFactors = FALSE
      ))

      # Record per-K time (biclustering only)
      if (!is.null(per_k)) {
        for (ki in seq_along(per_k)) {
          perk_results <- rbind(perk_results, data.frame(
            m = sz, p = sz, n = 2 * sz, rep = rep_i,
            method   = meth_name,
            K        = as.integer(names(per_k)[ki]),
            time_sec = per_k[ki],
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    gc()
  }
}

write.csv(comp_results, file.path(output_dir, "method_comparison.csv"),
          row.names = FALSE)
write.csv(perk_results, file.path(output_dir, "biclustering_per_k.csv"),
          row.names = FALSE)

# ==============================================================================
# Experiment 3: Non-square scaling (fixed p, varying m)
# ==============================================================================
cat("====== Experiment 3: Non-square scaling ======\n")

fixed_p <- 100
vary_m  <- c(50, 100, 200, 300, 500, 750, 1000)
nonsq_results <- data.frame()

for (m_val in vary_m) {
  cat(sprintf("  m=%d, p=%d ...\n", m_val, fixed_p))
  for (rep_i in 1:n_rep) {
    mat <- generate_block_matrix(m = m_val, p = fixed_p, seed = 42 + rep_i)
    tryCatch({
      res <- biser_timed(mat, simmeth = "t", noise = TRUE)
      for (step_name in names(res$timings)) {
        nonsq_results <- rbind(nonsq_results, data.frame(
          m = m_val, p = fixed_p, n = m_val + fixed_p,
          rep = rep_i, step = step_name,
          time_sec = res$timings[[step_name]],
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      cat(sprintf("    Error: %s\n", e$message))
    })
  }
  gc()
}

write.csv(nonsq_results, file.path(output_dir, "nonsquare_timings.csv"),
          row.names = FALSE)

# ==============================================================================
# Experiment 4: Memory usage analysis
# ==============================================================================
cat("====== Experiment 4: Memory analysis ======\n")

mem_sizes <- c(50, 100, 200, 300, 500, 750, 1000)
memory_results <- data.frame()

for (sz in mem_sizes) {
  cat(sprintf("  m=p=%d ...\n", sz))
  mat <- generate_block_matrix(m = sz, p = sz, seed = 42)
  n <- 2 * sz; r <- sz

  # Theoretical memory (bytes)
  mem_input <- sz * sz * 8
  mem_sim   <- n * n * 8
  mem_embed <- n * r * 8
  mem_svd   <- (sz * r + sz * r + r) * 8
  mem_dist  <- n * (n - 1) / 2 * 8
  mem_total_theory <- mem_input + mem_sim + mem_embed + mem_svd + mem_dist

  gc(reset = TRUE)
  mem_before <- gc(full = TRUE)[2, 2]

  tryCatch({
    res <- biser_timed(mat, simmeth = "t", noise = TRUE)
    mem_after <- gc(full = TRUE)[2, 2]
    memory_results <- rbind(memory_results, data.frame(
      m = sz, p = sz, n = n,
      mem_theory_MB = mem_total_theory / 1024^2,
      mem_actual_MB = max(0, mem_after - mem_before),
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    cat(sprintf("    Error: %s\n", e$message))
  })
  rm(res, mat); gc()
}

write.csv(memory_results, file.path(output_dir, "memory_usage.csv"),
          row.names = FALSE)

# ==============================================================================
# Complexity fitting: T = a * n^b
# ==============================================================================
cat("====== Complexity fitting ======\n")

stepwise_avg <- stepwise_results %>%
  filter(step != "total") %>%
  group_by(n, step) %>%
  summarise(mean_time = mean(time_sec), .groups = "drop") %>%
  filter(mean_time > 0)

complexity_fits <- stepwise_avg %>%
  group_by(step) %>%
  do({
    df_sub <- .
    if (nrow(df_sub) >= 3 && all(df_sub$mean_time > 0)) {
      fit <- lm(log(mean_time) ~ log(n), data = df_sub)
      data.frame(
        empirical_complexity = sprintf("O(n^%.2f)", coef(fit)[2]),
        exponent  = round(coef(fit)[2], 2),
        R_squared = round(summary(fit)$r.squared, 3),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(empirical_complexity = "insufficient data",
                 exponent = NA_real_, R_squared = NA_real_,
                 stringsAsFactors = FALSE)
    }
  }) %>% ungroup()

cat("\n--- BiSer step-wise complexity ---\n")
print(as.data.frame(complexity_fits))
write.csv(complexity_fits, file.path(output_dir, "complexity_fits.csv"),
          row.names = FALSE)

# Cross-method complexity
method_complexity <- comp_results %>%
  filter(!is.na(time_sec)) %>%
  group_by(n, method) %>%
  summarise(mean_time = mean(time_sec), .groups = "drop") %>%
  filter(mean_time > 0) %>%
  group_by(method) %>%
  do({
    df_sub <- .
    if (nrow(df_sub) >= 3) {
      fit <- lm(log(mean_time) ~ log(n), data = df_sub)
      data.frame(
        empirical_complexity = sprintf("O(n^%.2f)", coef(fit)[2]),
        exponent  = round(coef(fit)[2], 2),
        R_squared = round(summary(fit)$r.squared, 3),
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(empirical_complexity = "insufficient data",
                 exponent = NA_real_, R_squared = NA_real_,
                 stringsAsFactors = FALSE)
    }
  }) %>% ungroup()

cat("\n--- Cross-method complexity ---\n")
print(as.data.frame(method_complexity))
write.csv(method_complexity, file.path(output_dir, "method_complexity_fits.csv"),
          row.names = FALSE)

# ==============================================================================
# Summary tables
# ==============================================================================
# Table 1: BiSer step-wise
table1 <- stepwise_results %>%
  group_by(n, step) %>%
  summarise(mean = mean(time_sec), sd = sd(time_sec), .groups = "drop") %>%
  mutate(value = sprintf("%.4f +/- %.4f", mean, sd)) %>%
  select(n, step, value) %>%
  pivot_wider(names_from = step, values_from = value)
write.csv(table1, file.path(output_dir, "table_stepwise.csv"), row.names = FALSE)

# Table 2: Method comparison
table2 <- comp_results %>%
  filter(!is.na(time_sec)) %>%
  group_by(n, method, method_type) %>%
  summarise(mean = mean(time_sec), sd = sd(time_sec), .groups = "drop") %>%
  mutate(value = sprintf("%.3f +/- %.3f", mean, sd)) %>%
  select(n, method, value) %>%
  pivot_wider(names_from = method, values_from = value)
write.csv(table2, file.path(output_dir, "table_methods.csv"), row.names = FALSE)

# Session info
sink(file.path(output_dir, "session_info.txt"))
cat("Benchmark completed at:", format(Sys.time()), "\n\n")
sessionInfo()
sink()

cat("\n========================================\n")
cat("  All results saved to:", output_dir, "\n")
cat("========================================\n")
cat("Done!\n")
