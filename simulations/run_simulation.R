################################################################################
# BiSer Benchmarking: Unified Simulation Runner
#
# Parameterized simulation engine that replaces the duplicated code in
# simu_low.R, simu_mid.R, simu_high.R. Call with different config files
# to run different scenarios.
#
# Usage:
#   source("R/methods.R")
#   source("R/data_generation.R")
#   source("R/metrics.R")
#   source("simulations/config_low_noise.R")   # loads 'configs' list
#   source("simulations/run_simulation.R")      # runs all configs
################################################################################

library(reticulate)
library(abind)
library(gtools)

# --- Load Python boundary detection ---
# Adjust path as needed
source_python("python/auto_boundaries.py")

# ==============================================================================
# Main simulation loop (runs for a single config)
# ==============================================================================
run_single_simulation <- function(config, output_dir = "output") {

  cat("\n========================================\n")
  cat("Running:", config$name, "\n")
  cat("========================================\n")

  # Unpack parameters
  bicrnum    <- config$bicrnum
  biccnum    <- config$biccnum
  overr      <- config$overr
  overc      <- config$overc
  prominence <- config$prominence
  distance   <- config$distance
  n_iter     <- config$n_iter
  gen_fn     <- config$generator     # function: generate_norm, generate_NBP, etc.
  gen_args   <- config$gen_args      # list of extra args (bicmean, bicsd, etc.)

  # Metric names
  metricname <- c("NMI_row", "NMI_col", "NMI", "purity_row", "purity_col", "purity",
                  "pathlen_row", "pathlen_col", "pathlen", "ARI_row", "ARI_col", "ARI",
                  "smooth_row", "smooth_col", "smooth",
                  "bandwidth_row", "bandwidth_col", "bandwidth",
                  "block_contrast_row", "block_contrast_col", "block_contrast")

  allout <- allclus <- allmetric <- trueinfo <- list()

  for (iter in 1:n_iter) {
    cat("  Iteration:", iter, "/", n_iter, "\r")

    # Generate data
    gen_call_args <- c(list(bicrnum = bicrnum, biccnum = biccnum,
                            overr = overr, overc = overc), gen_args)
    mat <- do.call(gen_fn, gen_call_args)
    m <- nrow(mat); p <- ncol(mat)
    rownames(mat) <- paste0("r", 1:m)
    colnames(mat) <- paste0("c", 1:p)

    # Shuffle
    indr <- sample(1:m, m)
    indc <- sample(1:p, p)
    mat1 <- mat[indr, indc]
    rsim <- cor(t(mat1))
    csim <- cor(mat1)

    # Ground truth labels
    truelabel <- generate_true_labels(bicrnum, biccnum, overr, overc, indr, indc)
    trueinfo[[iter]] <- list(mat = mat, label = truelabel)

    # Initialize storage
    allclus[[iter]] <- allout[[iter]] <- list()
    allmetric[[iter]] <- matrix(NA, nrow = length(methname), ncol = length(metricname),
                                dimnames = list(methname, metricname))

    # --- BiSer ---
    out <- biser(mat1, noise = TRUE, pct = 0.3)
    allout[[iter]]$biser <- out
    simmat <- out$sim
    simmat2 <- simmat[c(rownames(out$reordered_mat), colnames(out$reordered_mat)),
                      c(rownames(out$reordered_mat), colnames(out$reordered_mat))]
    py_sim <- r_to_py(simmat2)
    auto_boundaries <- find_auto_boundaries(py_sim, valley = "find_peaks",
                                             prominence = prominence,
                                             distance = distance)
    label <- labelinput(rownames(simmat2),
                        list(row = rownames(mat1), col = colnames(mat1)),
                        auto_boundaries)
    allmetric[[iter]]["biser", ] <- overlapbic_metric(label, truelabel, rsim, csim, m, p)
    allclus[[iter]][["biser"]] <- label

    # --- Spectral seriation ---
    out <- spec_seri(mat1)
    allout[[iter]]$spec_seri <- out
    out2 <- seriout(out, truelabel, rsim, csim, m, p, mat1,
                    r_to_py, find_auto_boundaries, prominence, distance)
    allmetric[[iter]]["spec_seri", ] <- out2$metric
    allclus[[iter]][["spec_seri"]] <- out2$clus

    # --- TSP seriation ---
    out <- run_tsp_seriation(mat1)
    allout[[iter]]$tsp_seri <- out
    out2 <- seriout(out, truelabel, rsim, csim, m, p, mat1,
                    r_to_py, find_auto_boundaries, prominence, distance)
    allmetric[[iter]]["tsp_seri", ] <- out2$metric
    allclus[[iter]][["tsp_seri"]] <- out2$clus

    # --- R seriation methods ---
    for (i in Rserimeth) {
      out <- Rseriation(mat1, i)
      allout[[iter]][[i]] <- out
      out2 <- seriout(out, truelabel, rsim, csim, m, p, mat1,
                      r_to_py, find_auto_boundaries, prominence, distance)
      allmetric[[iter]][i, ] <- out2$metric
      allclus[[iter]][[i]] <- out2$clus
    }

    # --- Yang BS / BS2 ---
    for (i in c("bs", "bs2")) {
      out <- Yang_bs(mat1, i)
      allout[[iter]][[i]] <- out
      out2 <- seriout(out, truelabel, rsim, csim, m, p, mat1,
                      r_to_py, find_auto_boundaries, prominence, distance)
      allmetric[[iter]][i, ] <- out2$metric
      allclus[[iter]][[i]] <- out2$clus
    }

    # --- Biclustering methods ---
    for (i in bicmeth) {
      out <- bicluster(mat1, i)
      if (length(out[["row_order"]]) == 0) {
        allout[[iter]][[i]] <- allclus[[iter]][[i]] <- NA
      } else {
        allout[[iter]][[i]] <- out
        try({
          allclus[[iter]][[i]] <- list(row = out$row_label, col = out$col_label)
          allmetric[[iter]][i, ] <- overlapbic_metric(allclus[[iter]][[i]],
                                                      truelabel, rsim, csim, m, p)
        }, silent = TRUE)
      }
    }

    # --- NMF ---
    nmf_out <- run_nmf(mat1)
    allout[[iter]][["NMF"]] <- nmf_out
    if (nmf_out$success) {
      nmf_label <- list(row = nmf_out$row_cluster, col = nmf_out$col_cluster)
      try({
        allmetric[[iter]]["NMF", ] <- overlapbic_metric(nmf_label, truelabel,
                                                         rsim, csim, m, p)
        allclus[[iter]][["NMF"]] <- nmf_label
      }, silent = TRUE)
    }
  }
  cat("\n")

  # Save results
  ALLout <- list(trueinfo = trueinfo, allout = allout,
                 allclus = allclus, allmetric = allmetric)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  save(ALLout, file = file.path(output_dir, paste0("ALLout_", config$name, ".RData")))
  cat("Saved:", file.path(output_dir, paste0("ALLout_", config$name, ".RData")), "\n")

  return(ALLout)
}

# ==============================================================================
# Batch runner: iterate over all configs
# ==============================================================================
if (exists("configs") && is.list(configs)) {
  for (cfg in configs) {
    run_single_simulation(cfg, output_dir = "output")
  }
  cat("\nAll simulations completed.\n")
}
