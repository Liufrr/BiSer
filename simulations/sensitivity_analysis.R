################################################################################
# BiSer Boundary Detection Hyperparameter Sensitivity Analysis
#
# Design: Reuse saved BiSer seriation results; only re-run boundary detection
# Usage: setwd("path/to/BiSer"); source("simulations/sensitivity_analysis.R")
# Dependencies: reticulate (+ Python scipy), ggplot2, dplyr, tidyr, viridis, cowplot
################################################################################

library(reticulate)
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(cowplot)
library(mclust)
library(aricode)

# ============================================================
# 0. Load dependencies
# ============================================================
source("R/metrics.R")
source_python("python/auto_boundaries.py")

cat("Python boundary detection functions loaded\n")

# ============================================================
# 1. Parameter grid (3 parameters: window, prominence, distance)
#    BiSer0608 fixes Gaussian smoothing at sigma=3 and sweeps 5 x 7 x 5 = 175
#    combinations of the remaining boundary-detection parameters.
# ============================================================
param_grid <- expand.grid(
  window     = c(5, 8, 10, 15, 20),
  prominence = c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1),  # extended down to 0.001
  distance   = c(3, 5, 8, 10, 15),
  stringsAsFactors = FALSE
)
cat("Total parameter combinations:", nrow(param_grid), "\n")

# ============================================================
# 2. Core: parameter sweep for a single iteration
# ============================================================
run_sensitivity_single <- function(simmat2, mat1, truelabel, rsim, csim,
                                    param_grid) {
  m <- nrow(mat1)
  p <- ncol(mat1)
  sim_names <- rownames(simmat2)
  mat_names <- list(row = rownames(mat1), col = colnames(mat1))

  results <- data.frame(
    window     = param_grid$window,
    prominence = param_grid$prominence,
    distance   = param_grid$distance,
    ARI        = NA_real_,
    NMI        = NA_real_,
    purity     = NA_real_,
    n_clusters = NA_integer_
  )

  for (i in 1:nrow(param_grid)) {
    tryCatch({
      # Use Python find_auto_boundaries (verified correct implementation)
      py_sim <- r_to_py(simmat2)
      boundaries <- find_auto_boundaries(
        py_sim,
        valley     = "find_peaks",
        smooth     = "gaussian",
        sigma      = 3,
        window     = as.integer(param_grid$window[i]),
        prominence = param_grid$prominence[i],
        distance   = as.integer(param_grid$distance[i])
      )

      # Convert to R vector
      boundaries <- as.integer(boundaries)

      # Generate labels from boundaries
      label <- labelinput(sim_names, mat_names, boundaries)

      # Compute metrics
      metric <- overlapbic_metric(label, truelabel, rsim, csim, m, p)

      results$ARI[i]        <- metric["ARI"]
      results$NMI[i]        <- metric["NMI"]
      results$purity[i]     <- metric["purity"]
      results$n_clusters[i] <- max(length(unique(label$row)),
                                    length(unique(label$col)))
    }, error = function(e) {
      # Record as NA
    })
  }
  return(results)
}

# ============================================================
# 3. Batch-load RData files and run parameter sweep
# ============================================================
run_sensitivity_on_saved <- function(rdata_path, n_iter_use = 20,
                                      param_grid, setting_name = "",
                                      orig_prominence = 0.01,
                                      orig_distance = 5) {
  cat("\nLoading:", rdata_path, "\n")
  load(rdata_path)

  n_iter <- min(n_iter_use, length(ALLout$allout))
  cat("Using iterations:", n_iter, "\n")

  all_results <- list()
  for (iter in 1:n_iter) {
    cat("  iter", iter, "/", n_iter, "\r")

    biser_out <- ALLout$allout[[iter]]$biser
    if (is.null(biser_out) || is.null(biser_out$sim)) next

    info <- ALLout$trueinfo[[iter]]
    mat <- info$mat
    m <- nrow(mat); p <- ncol(mat)

    # Reconstruct simmat2 (BiSer-reordered similarity matrix)
    simmat <- biser_out$sim
    joint_order <- biser_out$joint_order
    if (is.null(joint_order)) {
      reordered <- biser_out$reordered_mat
      joint_order <- c(rownames(reordered), colnames(reordered))
    }
    simmat2 <- simmat[joint_order, joint_order, drop = FALSE]

    # New BiSer0608 outputs store the exact disturbed input and its aligned
    # ground truth. The fallback keeps compatibility with older RData files.
    if (!is.null(info$mat_input)) {
      mat1 <- info$mat_input
      truelabel <- info$label
    } else {
      indr <- match(rownames(simmat)[seq_len(m)], rownames(mat))
      indc <- match(rownames(simmat)[m + seq_len(p)], colnames(mat))
      mat1 <- mat[indr, indc, drop = FALSE]
      truelabel <- list(row = info$label$row[indr], col = info$label$col[indc])
    }
    rsim <- cor(t(mat1))
    csim <- cor(mat1)

    # Parameter sweep
    res <- run_sensitivity_single(simmat2, mat1, truelabel, rsim, csim,
                                   param_grid)
    res$iter <- iter
    res$setting <- setting_name

    # Validation: compare ARI using this setting's original parameters
    default_row <- which(abs(res$window - 10) < 0.01 &
                         abs(res$prominence - orig_prominence) < 1e-6 &
                         abs(res$distance - orig_distance) < 0.01)
    if (length(default_row) > 0 && iter == 1) {
      orig_ari <- ALLout$allmetric[[iter]]["biser", "ARI"]
      new_ari  <- res$ARI[default_row[1]]
      cat(sprintf("\n  Validation [prom=%.3f, dist=%d]: original ARI=%.3f, recomputed ARI=%.3f\n",
                  orig_prominence, orig_distance, orig_ari, new_ari))
    }

    all_results[[iter]] <- res
  }
  cat("\n")
  return(bind_rows(all_results))
}

# ============================================================
# 4. Run
# ============================================================

# Original prominence / distance parameters for each setting
settings <- list(
  list(path = file.path("output", "ALLout_norm_low.RData"),
       name = "norm_low",           prom = 0.01, dist = 5),
  list(path = file.path("output", "ALLout_norm_exc_low.RData"),
       name = "norm_exc_low",       prom = 0.01, dist = 5),
  list(path = file.path("output", "ALLout_norm_low_1000.RData"),
       name = "norm_low_1000",      prom = 0.01, dist = 5),
  list(path = file.path("output", "ALLout_norm_exc_low_1000.RData"),
       name = "norm_exc_low_1000",  prom = 0.01, dist = 5)
)

all_sensitivity <- list()
for (s in settings) {
  if (!file.exists(s$path)) { cat("Skipping:", s$path, "\n"); next }
  all_sensitivity[[s$name]] <- run_sensitivity_on_saved(
    s$path, n_iter_use = 20, param_grid = param_grid, setting_name = s$name,
    orig_prominence = s$prom, orig_distance = s$dist
  )
}

df_all <- bind_rows(all_sensitivity)
if (nrow(df_all) == 0) stop("No valid results! Please check RData paths")

# ============================================================
# 5. Summary and visualization
# ============================================================

# Setting name mapping (for display-friendly labels in figures)
setting_labels <- c(
  "norm_low"           = "Overlap, 250×150",
  "norm_exc_low"       = "Exclusive, 250×150",
  "norm_low_1000"      = "Overlap, 1000×650",
  "norm_exc_low_1000"  = "Exclusive, 1000×650"
)

df_all$setting_label <- setting_labels[df_all$setting]
# Set factor order: small matrix first, overlap first
df_all$setting_label <- factor(df_all$setting_label,
  levels = c("Overlap, 250×150", "Exclusive, 250×150",
             "Overlap, 1000×650", "Exclusive, 1000×650"))

# --- 5.1 Global summary ---
cat("\n========== Global Sensitivity Summary ==========\n")
overall <- df_all %>%
  group_by(setting) %>%
  summarise(
    ARI_mean     = mean(ARI, na.rm = TRUE),
    ARI_median   = median(ARI, na.rm = TRUE),
    ARI_sd       = sd(ARI, na.rm = TRUE),
    ARI_min      = min(ARI, na.rm = TRUE),
    ARI_max      = max(ARI, na.rm = TRUE),
    pct_above_09 = mean(ARI >= 0.9, na.rm = TRUE) * 100,
    pct_above_08 = mean(ARI >= 0.8, na.rm = TRUE) * 100,
    .groups = "drop"
  )
print(as.data.frame(overall))

# --- 5.2 Prominence × distance heatmap function ---
make_heatmap_pd <- function(data, setting_lab, w_val = 10) {
  plot_data <- data %>%
    filter(setting_label == setting_lab, window == w_val) %>%
    group_by(prominence, distance) %>%
    summarise(ARI_mean = mean(ARI, na.rm = TRUE), .groups = "drop")
  if (nrow(plot_data) == 0) return(NULL)

  ggplot(plot_data, aes(x = factor(prominence), y = factor(distance),
                         fill = ARI_mean)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", ARI_mean)),
              size = 3.2,
              color = ifelse(plot_data$ARI_mean > 0.5, "black", "white")) +
    scale_fill_viridis(option = "D", limits = c(0, 1), name = "Mean ARI") +
    labs(x = "Prominence threshold (p)",
         y = "Minimum distance (d)",
         title = setting_lab) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
          axis.text = element_text(color = "black", size = 9),
          axis.title = element_text(size = 10),
          legend.key.height = unit(0.4, "cm"),
          panel.grid = element_blank())
}

# --- 5.3 Marginal effect plot function ---
make_marginal_plot <- function(data) {
  data %>%
    pivot_longer(cols = c(window, prominence, distance),
                 names_to = "parameter", values_to = "value") %>%
    mutate(parameter = factor(parameter,
      levels = c("prominence", "distance", "window"),
      labels = c("Prominence threshold (p)",
                  "Minimum distance (d)",
                  "Window width (w)"))) %>%
    group_by(setting_label, parameter, value) %>%
    summarise(ARI_mean = mean(ARI, na.rm = TRUE),
              ARI_sd = sd(ARI, na.rm = TRUE), .groups = "drop") %>%
    ggplot(aes(x = factor(value), y = ARI_mean,
               color = setting_label, group = setting_label)) +
    geom_point(size = 2) +
    geom_line(linewidth = 0.7) +
    geom_errorbar(aes(ymin = pmax(ARI_mean - ARI_sd, 0),
                       ymax = pmin(ARI_mean + ARI_sd, 1)),
                  width = 0.2, linewidth = 0.3) +
    facet_wrap(~ parameter, scales = "free_x", nrow = 1) +
    scale_color_viridis_d(option = "H", name = "Setting") +
    labs(x = "Parameter value", y = "Mean ARI ± SD") +
    theme_classic(base_size = 12) +
    theme(strip.text = element_text(face = "bold", size = 11),
          strip.background = element_blank(),
          legend.position = "bottom",
          legend.text = element_text(size = 9),
          axis.text = element_text(size = 9)) +
    guides(color = guide_legend(nrow = 1))
}

# --- 5.4 Generate Supplementary Figure: heatmaps (A-D) + marginal effects (E) ---
setting_labs <- levels(df_all$setting_label)

pd_plots <- lapply(setting_labs, function(s) make_heatmap_pd(df_all, s))
pd_plots <- pd_plots[!sapply(pd_plots, is.null)]

if (length(pd_plots) > 0) {
  # Upper panel: 4 heatmaps (2×2)
  fig_heatmaps <- plot_grid(plotlist = pd_plots, ncol = 2,
                             labels = LETTERS[1:length(pd_plots)],
                             label_size = 16, label_fontface = "bold")

  # Lower panel: marginal effect plot
  fig_marginal <- make_marginal_plot(df_all)

  # Combined: heatmaps 65%, marginal effects 35%
  fig_combined <- plot_grid(
    fig_heatmaps,
    plot_grid(fig_marginal, labels = "E",
              label_size = 16, label_fontface = "bold"),
    ncol = 1,
    rel_heights = c(0.65, 0.35)
  )

  ggsave(file.path("output", "FigS_sensitivity.png"), fig_combined,
         width = 13, height = 14, dpi = 400, bg = "white")
  cat("Saved: FigS_sensitivity.png\n")
}

# --- 5.5 Save data ---
write.csv(df_all, file.path("output", "sensitivity_raw_results.csv"), row.names = FALSE)

# Paper-ready summary
cat("\n========== Paper-Ready Summary ==========\n")
for (sname in unique(df_all$setting)) {
  sub <- df_all %>% filter(setting == sname)
  cat(sprintf(
    "[%s] ARI: mean=%.3f, median=%.3f, range=[%.3f,%.3f], >=0.8: %.1f%% (%d combos × %d iters)\n",
    sname, mean(sub$ARI, na.rm=T), median(sub$ARI, na.rm=T),
    min(sub$ARI, na.rm=T), max(sub$ARI, na.rm=T),
    mean(sub$ARI >= 0.8, na.rm=T)*100,
    nrow(param_grid), length(unique(sub$iter))
  ))
}
cat("\nDone!\n")
