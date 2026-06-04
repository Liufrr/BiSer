################################################################################
# Metric Comparison Script: Compare the effect of different Spearman
# computation methods on method rankings
#
# Prerequisite: analysis_yan_unsupervised.R has been run in the R environment,
#               so allout, mat1, mat1_labels, group_labels, etc. are available.
#
# This script computes 4 metrics:
#   1. Per-cell Spearman (current scheme)
#   2. Kendall's tau
#   3. Stage-median Spearman
#   4. Stage-mean Spearman
#
# Output:
#   - 4 scatter plots arranged side by side for comparison
#   - Summary table printed to the console
#   - Saved as compare_metrics.png
################################################################################

library(ggplot2)
library(ggrepel)
library(cowplot)
library(dplyr)
library(mclust)
library(combinat)

# ============================================================
# 0. Configuration (consistent with analysis_yan_unsupervised.R)
# ============================================================

method_keys <- c("biser", "bs", "tsp_seri", "spec_seri",
                 "Heatmap", "MESBC", "NMF")

display_name <- c(
  "biser"     = "BiSer",
  "bs"        = "BiSer_SVD",
  "tsp_seri"  = "BiSer_TSP",
  "spec_seri" = "Spectral",
  "Heatmap"   = "Heatmap",
  "MESBC"     = "MESBC",
  "NMF"       = "NMF"
)

stage_levels <- c("Oocyte", "Zygote", "Two-cell", "Four-cell",
                  "Eight-cell", "Morula", "Blastocyst")


# ============================================================
# 1. Metric Computation Functions
# ============================================================

calculate_all_metrics <- function(method_name, out_list,
                                  true_labels, stages) {
  pred_order <- out_list[[method_name]]$col_order
  if (is.null(pred_order)) return(NULL)

  # --- Labels and numeric time codes sorted by predicted order ---
  sorted_labels <- true_labels[pred_order]
  numeric_time  <- as.numeric(factor(sorted_labels, levels = stages))
  pos           <- seq_along(pred_order)

  # --- 1. Per-cell Spearman (current scheme) ---
  spearman_cell <- abs(cor(pos, numeric_time, method = "spearman",
                           use = "complete.obs"))

  # --- 2. Kendall's tau ---
  kendall_tau <- abs(cor(pos, numeric_time, method = "kendall",
                         use = "complete.obs"))

  # --- 3. Stage-median Spearman ---
  # Median position of each stage in the predicted ordering
  stage_median <- tapply(pos, sorted_labels, median)
  stage_median <- stage_median[stages]
  # Compute only for stages with data
  valid <- !is.na(stage_median)
  if (sum(valid) >= 3) {
    spearman_median <- abs(cor(stage_median[valid],
                               seq_along(stages)[valid],
                               method = "spearman"))
  } else {
    spearman_median <- NA
  }

  # --- 4. Stage-mean Spearman ---
  stage_mean <- tapply(pos, sorted_labels, mean)
  stage_mean <- stage_mean[stages]
  valid2 <- !is.na(stage_mean)
  if (sum(valid2) >= 3) {
    spearman_mean <- abs(cor(stage_mean[valid2],
                              seq_along(stages)[valid2],
                              method = "spearman"))
  } else {
    spearman_mean <- NA
  }

  # --- ARI (same as current code) ---
  true_sizes   <- as.numeric(table(true_labels))
  num_clusters <- length(true_sizes)
  all_perms    <- combinat::permn(true_sizes)
  max_ari <- -1
  for (sp in all_perms) {
    pred_clusters <- rep(seq_len(num_clusters), times = sp)
    ari <- adjustedRandIndex(pred_clusters, sorted_labels)
    if (ari > max_ari) max_ari <- ari
  }

  data.frame(
    Method           = method_name,
    Spearman_Cell    = spearman_cell,
    Kendall_Tau      = kendall_tau,
    Spearman_Median  = spearman_median,
    Spearman_Mean    = spearman_mean,
    ARI              = max_ari,
    stringsAsFactors = FALSE
  )
}


# ============================================================
# 2. Compute Metrics for All Methods
# ============================================================

cat("Computing metrics for all methods...\n")

metrics_all <- bind_rows(lapply(method_keys, calculate_all_metrics,
                                out_list    = allout,
                                true_labels = mat1_labels,
                                stages      = stage_levels))

metrics_all$Display <- display_name[metrics_all$Method]
metrics_all$Highlight <- ifelse(metrics_all$Display == "BiSer",
                                "Target", "Others")

# ============================================================
# 3. Print Summary Table
# ============================================================

cat("\n")
cat("========================================================\n")
cat("         Full Comparison Table: Methods x Metrics\n")
cat("========================================================\n\n")

print_df <- metrics_all %>%
  mutate(across(where(is.numeric), ~sprintf("%.4f", .))) %>%
  select(Display, Spearman_Cell, Kendall_Tau,
         Spearman_Median, Spearman_Mean, ARI)
colnames(print_df) <- c("Method", "Per-cell\nSpearman",
                         "Kendall\ntau",
                         "Stage-median\nSpearman",
                         "Stage-mean\nSpearman", "ARI")
print(print_df, row.names = FALSE, right = FALSE)

cat("\n")

# --- BiSer ranking across metrics ---
cat("--- BiSer Ranking ---\n")
for (metric in c("Spearman_Cell", "Kendall_Tau",
                  "Spearman_Median", "Spearman_Mean")) {
  vals <- metrics_all[[metric]]
  biser_val <- vals[metrics_all$Method == "biser"]
  rank_pos <- sum(vals > biser_val, na.rm = TRUE) + 1
  cat(sprintf("  %-20s: %.4f (rank %d / %d)\n",
              metric, biser_val, rank_pos, sum(!is.na(vals))))
}
cat("\n")


# ============================================================
# 4. Draw 4 Scatter Plots for Comparison
# ============================================================

make_scatter <- function(df, x_col, x_label, title) {
  med_x <- median(df[[x_col]], na.rm = TRUE)
  med_y <- median(df$ARI, na.rm = TRUE)

  ggplot(df, aes(x = .data[[x_col]], y = ARI)) +
    geom_hline(yintercept = med_y, linetype = "dashed",
               color = "gray70", linewidth = 0.4) +
    geom_vline(xintercept = med_x, linetype = "dashed",
               color = "gray70", linewidth = 0.4) +
    geom_point(aes(color = Highlight, size = Highlight, shape = Highlight),
               alpha = 0.85) +
    scale_color_manual(values = c("Target" = "#DC0000FF",
                                  "Others" = "#7E6148FF")) +
    scale_size_manual(values  = c("Target" = 5, "Others" = 3)) +
    scale_shape_manual(values = c("Target" = 16, "Others" = 17)) +
    geom_text_repel(aes(label = Display, color = Highlight),
                    size = 5, fontface = "bold",
                    box.padding = 0.5, point.padding = 0.3,
                    max.overlaps = 20, seed = 42) +
    labs(x = x_label, y = "Clustering Purity (ARI)",
         title = title) +
    theme_classic(base_size = 16) +
    theme(legend.position = "none",
          plot.title = element_text(face = "bold", size = 16,
                                   hjust = 0.5),
          axis.title = element_text(face = "bold", size = 14),
          axis.text  = element_text(color = "black", size = 12),
          axis.line  = element_line(linewidth = 0.6))
}

p1 <- make_scatter(metrics_all, "Spearman_Cell",
                   "Per-cell Spearman",
                   "Scheme 1: Per-cell Spearman (current)")

p2 <- make_scatter(metrics_all, "Kendall_Tau",
                   "Kendall's tau",
                   "Scheme 2: Kendall's tau")

p3 <- make_scatter(metrics_all, "Spearman_Median",
                   "Stage-median Spearman",
                   "Scheme 3: Stage-median Spearman")

p4 <- make_scatter(metrics_all, "Spearman_Mean",
                   "Stage-mean Spearman",
                   "Scheme 4: Stage-mean Spearman")


# ============================================================
# 5. Combine and Save
# ============================================================

fig_compare <- plot_grid(p1, p2, p3, p4,
                         ncol = 2, nrow = 2,
                         labels = c("A", "B", "C", "D"),
                         label_size = 20)

ggsave("compare_metrics.png", fig_compare,
       width = 18, height = 14, dpi = 300)

cat("========================================================\n")
cat("Figure saved: compare_metrics.png\n")
cat("========================================================\n")
