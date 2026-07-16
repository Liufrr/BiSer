################################################################################
# BiSer0608 unified simulation runner
################################################################################

library(reticulate)

source_python("python/auto_boundaries.py")

required_config_fields <- c(
  "name", "generator", "n_iter", "seed", "bicrnum", "biccnum",
  "overr", "overc", "gen_args", "window", "smooth", "sigma",
  "prominence", "distance"
)

validate_simulation_config <- function(config) {
  missing <- setdiff(required_config_fields, names(config))
  if (length(missing) > 0L) {
    stop("Configuration ", config$name %||% "<unnamed>",
         " is missing: ", paste(missing, collapse = ", "))
  }
  lengths <- vapply(
    config[c("bicrnum", "biccnum", "overr", "overc")], length, integer(1)
  )
  if (length(unique(lengths)) != 1L) {
    stop("Block-size and overlap vectors must have equal lengths.")
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

safe_matrix_similarity <- function(mat, margin = c("row", "col")) {
  margin <- match.arg(margin)
  x <- if (identical(margin, "row")) mat else t(mat)
  sim <- suppressWarnings(cor(t(x), use = "pairwise.complete.obs"))
  sim[!is.finite(sim)] <- 0
  diag(sim) <- 1
  names <- if (identical(margin, "row")) rownames(mat) else colnames(mat)
  rownames(sim) <- colnames(sim) <- names
  sim
}

evaluate_biser <- function(out, mat1, truelabel, config) {
  joint_order <- out$joint_order
  reordered_sim <- out$sim[joint_order, joint_order, drop = FALSE]
  boundaries <- find_auto_boundaries(
    r_to_py(reordered_sim),
    valley = "find_peaks",
    smooth = config$smooth,
    window = as.integer(config$window),
    sigma = config$sigma,
    prominence = config$prominence,
    distance = as.integer(config$distance)
  )
  matname <- list(row = rownames(mat1), col = colnames(mat1))
  labels <- labelinput(joint_order, matname, boundaries)
  pre_labels <- oracle_joint_partition(joint_order, truelabel, matname)
  metric <- add_discretization_metrics(
    overlapbic_metric(labels, truelabel), pre_labels, truelabel
  )
  list(
    clus = labels, pre_clus = pre_labels, metric = metric,
    boundaries = as.integer(boundaries)
  )
}

run_single_simulation <- function(config, output_dir = "output") {
  validate_simulation_config(config)
  cat("\n========================================\n")
  cat("Running:", config$name, "\n")
  cat("========================================\n")

  n_iter <- as.integer(config$n_iter)
  allout <- allclus <- allpreclus <- allmetric <- trueinfo <- failures <-
    vector("list", n_iter)

  for (iter in seq_len(n_iter)) {
    cat("  Iteration:", iter, "/", n_iter, "\r")
    iteration_seed <- as.integer(config$seed + iter)
    set.seed(iteration_seed)

    generator_args <- c(
      list(
        bicrnum = config$bicrnum, biccnum = config$biccnum,
        overr = config$overr, overc = config$overc
      ),
      config$gen_args
    )
    mat <- do.call(config$generator, generator_args)
    m <- nrow(mat)
    p <- ncol(mat)
    rownames(mat) <- paste0("r", seq_len(m))
    colnames(mat) <- paste0("c", seq_len(p))

    row_permutation <- sample.int(m)
    col_permutation <- sample.int(p)
    mat1 <- mat[row_permutation, col_permutation, drop = FALSE]
    rsim <- safe_matrix_similarity(mat1, "row")
    csim <- safe_matrix_similarity(mat1, "col")
    truelabel <- generate_true_labels(
      config$bicrnum, config$biccnum, config$overr, config$overc,
      row_permutation, col_permutation
    )

    trueinfo[[iter]] <- list(
      mat = mat,
      mat_input = mat1,
      row_permutation = row_permutation,
      col_permutation = col_permutation,
      label = truelabel,
      seed = iteration_seed
    )
    allout[[iter]] <- list()
    allclus[[iter]] <- list()
    allpreclus[[iter]] <- list()
    failures[[iter]] <- list()
    allmetric[[iter]] <- matrix(
      NA_real_, nrow = length(submission_methods),
      ncol = length(submission_metric_names),
      dimnames = list(submission_methods, submission_metric_names)
    )

    record_failure <- function(method, error) {
      failures[[iter]][[method]] <<- conditionMessage(error)
      message("\n[", config$name, ", iteration ", iter, ", ", method,
              "] ", conditionMessage(error))
    }

    # BiSer: one joint order and one joint boundary profile.
    tryCatch({
      out <- biser(
        mat1, simmeth = "t", noise = TRUE,
        n_starts = config$n_starts %||% 100L
      )
      evaluated <- evaluate_biser(out, mat1, truelabel, config)
      out$boundaries <- evaluated$boundaries
      allout[[iter]]$biser <- out
      allclus[[iter]]$biser <- evaluated$clus
      allpreclus[[iter]]$biser <- evaluated$pre_clus
      allmetric[[iter]]["biser", ] <- evaluated$metric
    }, error = function(e) record_failure("biser", e))

    # Independent seriation and the two ablation variants. Row and column
    # profiles are segmented independently, as stated in the manuscript.
    ordering_methods <- list(
      bs = function() Yang_bs(mat1, "bs"),
      tsp_seri = function() run_tsp_seriation(mat1),
      spec_seri = function() spec_seri(mat1),
      Heatmap = function() Rseriation(mat1, "Heatmap")
    )
    for (method in names(ordering_methods)) {
      tryCatch({
        out <- ordering_methods[[method]]()
        evaluated <- seriout(
          out, truelabel, rsim, csim, m, p, mat1,
          r_to_py, find_auto_boundaries,
          prominence = config$prominence, distance = config$distance,
          window = config$window, smooth = config$smooth, sigma = config$sigma
        )
        out$boundaries <- evaluated$boundaries
        allout[[iter]][[method]] <- out
        allclus[[iter]][[method]] <- evaluated$clus
        allpreclus[[iter]][[method]] <- evaluated$pre_clus
        allmetric[[iter]][method, ] <- evaluated$metric
      }, error = function(e) record_failure(method, e))
    }

    # MESBC and NMF return discrete biclusters directly. A boundary-induced
    # discretization penalty is consequently not defined for these methods.
    tryCatch({
      out <- bicluster(mat1, "MESBC", K = config$candidate_k %||% 2:10)
      labels <- list(row = out$row_label, col = out$col_label)
      allout[[iter]]$MESBC <- out
      allclus[[iter]]$MESBC <- labels
      allmetric[[iter]]["MESBC", ] <- overlapbic_metric(labels, truelabel)
    }, error = function(e) record_failure("MESBC", e))

    tryCatch({
      out <- run_nmf(
        mat1, K = config$candidate_k %||% 2:10,
        nrun = config$nmf_nrun %||% 10L, seed = iteration_seed
      )
      allout[[iter]]$NMF <- out
      if (!isTRUE(out$success)) stop("NMF failed for every candidate rank.")
      labels <- list(row = out$row_label, col = out$col_label)
      allclus[[iter]]$NMF <- labels
      allmetric[[iter]]["NMF", ] <- overlapbic_metric(labels, truelabel)
    }, error = function(e) record_failure("NMF", e))
  }
  cat("\n")

  ALLout <- list(
    config = config,
    method_order = submission_methods,
    metric_definitions = list(
      separate_axis_combination = "row/column-size-weighted arithmetic mean",
      ARI_pre = "oracle segmentation using the true number and sizes of groups",
      ARI_post = "data-driven Gaussian-smoothed boundary segmentation",
      discretization_penalty = "ARI_pre - ARI_post"
    ),
    trueinfo = trueinfo,
    allout = allout,
    allclus = allclus,
    allpreclus = allpreclus,
    allmetric = allmetric,
    failures = failures
  )

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  output_file <- file.path(output_dir, paste0("ALLout_", config$name, ".RData"))
  save(ALLout, file = output_file)
  cat("Saved:", output_file, "\n")
  invisible(ALLout)
}

if (exists("configs") && is.list(configs) &&
    !isTRUE(getOption("BiSer.skip_autorun", FALSE))) {
  for (config in configs) run_single_simulation(config, output_dir = "output")
  cat("\nAll simulations completed.\n")
}
