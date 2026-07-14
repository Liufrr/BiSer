################################################################################
# BiSer Benchmarking: Yan Embryo Dataset — Unsupervised Feature Selection
#
# Single-cell RNA-seq analysis of the Yan human preimplantation embryo dataset.
# Uses unsupervised feature selection (HVG) instead of supervised methods.
#
# Produces:
#   - Heatmap comparisons (all methods, 2x4 grid with annotation bar)
#   - Spearman vs ARI scatter plot
#   - Marker gene ridge plot
#   - Gene module wave plot
#   - Combined Figure 6
#
# Required data:
#   yan.RData  — containing yan$data and yan$annotation$type
################################################################################

library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(viridis)
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(ggsci)
library(cowplot)
library(mclust)
library(combinat)
library(aricode)
library(TSP)
library(seriation)
library(scran)
library(igraph)
library(NMF)

# ============================================================
# 0. Global Settings
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

highlight_colors <- c("Target" = "#DC0000FF", "Others" = "#7E6148FF")

stage_levels <- c("Oocyte","Zygote", "Two-cell", "Four-cell",
                  "Eight-cell", "Morula", "Blastocyst")

stage_levels_clean <- c("Oocyte", "Zygote", "2-cell", "4-cell",
                        "8-cell", "Morula", "Blastocyst")


# ============================================================
# Source BiSer methods (biser, spec_seri, run_tsp_seriation, etc.)
# ============================================================
source("R/methods.R")

# ============================================================
# Unsupervised feature selection function
# ============================================================
source("R/feature_selection_unsupervised.R")

# ============================================================
# 1. Data Preparation
# ============================================================

# --- 1.1 Load Yan dataset ---
# NOTE: Adjust the path to your yan.RData file
load("data/yan.RData")

# --- 1.2 Unsupervised feature selection (scran HVG, top 200 genes) ---
n_top <- 300
feat <- get_clustering_matrix(yan$data, n_top = n_top, method = "hvg", norm = T)
top_mat <- feat$top_mat

# Z-score standardization (per gene)
final_genes <- as.matrix(top_mat)
class(final_genes) <- "numeric"
data_scaled <- t(scale(t(final_genes)))
data_scaled[is.na(data_scaled)] <- 0

group_labels <- as.factor(yan[["annotation"]][["type"]])
mat0 <- data_scaled

cat("Expression matrix:", nrow(mat0), "genes x", ncol(mat0), "samples\n")
cat("Group labels:", paste(levels(group_labels), collapse = ", "), "\n")

# --- 1.3 Shuffle matrix ---
mat1 <- mat0[sample(1:nrow(mat0), nrow(mat0)),
             sample(1:ncol(mat0), ncol(mat0))]

# ============================================================
# 2. Run all methods
# ============================================================
allout <- list()

# --- BiSer ---
out <- biser(mat1, noise = TRUE, simmeth = "t")
allout$biser <- out

# --- Spectral seriation ---
out <- spec_seri(mat1)
allout$spec_seri <- out

# --- TSP seriation ---
out <- run_tsp_seriation(mat1)
allout$tsp_seri <- out

# --- R seriation (Heatmap) ---
out <- Rseriation(mat1, "Heatmap")
allout$Heatmap <- out

# --- Yang BS ---
out <- Yang_bs(mat1 - min(mat1), "bs")
allout$bs <- out

# --- MESBC ---
out <- bicluster(mat1, "MESBC")
allout$MESBC <- out

# --- NMF ---
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

mat2 <- mat1 - min(mat1) + 1
nmfout <- list()
mod <- numeric()
g <- graph_out(mat1)
for (k in 2:10) {
  try({nmfout[[k]] <- nmf(mat2, rank = k)}, silent = TRUE)
  if (k > length(nmfout)) next
  row_cluster <- predict(nmfout[[k]], "rows")
  col_cluster <- predict(nmfout[[k]], "columns")
  mod[k] <- modularity(g, c(row_cluster, col_cluster))
}
res <- nmfout[[which.max(mod)]]
row_cluster <- predict(res, "rows")
col_cluster <- predict(res, "columns")
allout[["NMF"]] <- list(
  row_order = order(row_cluster), col_order = order(col_cluster),
  row_cluster = row_cluster, col_cluster = col_cluster
)


# ============================================================
# 3. Build mat1 using BiSer sim order (same as original code)
# ============================================================
aa <- allout[["biser"]][["sim"]]
mat1 <- mat0[rownames(aa)[1:nrow(mat0)],
             rownames(aa)[(nrow(mat0) + 1):nrow(aa)]]

# Derive cell labels for mat1 columns
col_shuffle_idx <- match(colnames(mat1), colnames(mat0))
mat1_labels     <- group_labels[col_shuffle_idx]
true_time_num   <- as.numeric(factor(mat1_labels, levels = stage_levels))


# ============================================================
# 4. Panel A: Heatmaps (2x4 grid) with true label annotations
# ============================================================

# --- Stage annotation colors ---
all_stages <- levels(group_labels)
if (is.null(all_stages)) all_stages <- sort(unique(as.character(group_labels)))

stage_palette <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                   "#F39B7FFF", "#8491B4FF", "#B09C85FF", "#7E6148FF")
stage_color_map <- setNames(stage_palette[seq_along(all_stages)], all_stages)

cat("Stage labels detected:", paste(all_stages, collapse = ", "), "\n")

# True labels for mat1 columns
mat1_col_labels <- as.character(mat1_labels)

# --- Heatmap grid function with top annotation bars ---
get_heatmap_grid <- function(mat_list, titles, col_label_list,
                             stage_colors, nrow = 3, ncol = 3,
                             color_palette = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                             legend_title = "Expression") {
  all_vals <- unlist(mat_list)
  max_abs  <- max(abs(all_vals), na.rm = TRUE)
  col_fun  <- colorRamp2(seq(-max_abs, max_abs,
                             length.out = length(color_palette)),
                         color_palette)

  expr_legend <- Legend(col_fun       = col_fun,
                        title         = legend_title,
                        title_gp      = gpar(fontsize = 26, fontface = "bold",
                                             fontfamily = "sans"),
                        labels_gp     = gpar(fontsize = 24, fontfamily = "sans"),
                        legend_height = unit(7, "cm"),
                        grid_width    = unit(1, "cm"),
                        title_gap     = unit(0.6, "cm"))

  stage_legend <- Legend(labels    = names(stage_colors),
                         legend_gp = gpar(fill = stage_colors),
                         title     = "Stage",
                         title_gp  = gpar(fontsize = 26, fontface = "bold",
                                          fontfamily = "sans"),
                         labels_gp = gpar(fontsize = 23, fontfamily = "sans"),
                         grid_width  = unit(0.8, "cm"),
                         grid_height = unit(0.8, "cm"),
                         row_gap     = unit(0.5, "cm"),
                         title_gap   = unit(0.6, "cm"))

  legend_grob <- grid.grabExpr({
    draw(packLegend(expr_legend, stage_legend, direction = "vertical",
                    gap = unit(2, "cm")))
  })

  ht_grobs <- lapply(seq_along(mat_list), function(i) {
    col_labels_i <- col_label_list[[i]]
    top_anno <- HeatmapAnnotation(
      Stage    = col_labels_i,
      col      = list(Stage = stage_colors),
      show_legend      = FALSE,
      show_annotation_name = FALSE,
      annotation_height = unit(0.35, "cm"),
      simple_anno_size  = unit(0.35, "cm")
    )
    ht <- Heatmap(mat_list[[i]],
                  col = col_fun,
                  top_annotation   = top_anno,
                  cluster_rows     = FALSE,
                  cluster_columns  = FALSE,
                  show_row_names   = FALSE,
                  show_column_names = FALSE,
                  column_title     = titles[i],
                  column_title_gp  = gpar(fontsize = 28, fontface = "bold",
                                          fontfamily = "sans"),
                  heatmap_legend_param = list(title = NULL),
                  rect_gp = gpar(col = NA, lwd = 0),
                  name    = paste0("ht_", i))
    grid.grabExpr(draw(ht, show_heatmap_legend = FALSE))
  })

  arrangeGrob(
    arrangeGrob(grobs = ht_grobs, nrow = nrow, ncol = ncol),
    legend_grob,
    nrow = 1, widths = c(4, 1.2)
  )
}

# --- Build ordered matrices and their column labels ---
m <- nrow(mat1)
p <- ncol(mat1)

mat_method_list <- lapply(allout[method_keys], function(x) {
  if (!is.null(x$row_order) &&
      length(x$row_order) == m && length(x$col_order) == p)
    return(mat1[x$row_order, x$col_order])
  return(NA)
})

mat_list <- c(list("Original" = mat1), mat_method_list)
heatmap_titles <- c("Original", display_name[method_keys])

col_label_list <- list(mat1_col_labels)
for (mk in method_keys) {
  x <- allout[[mk]]
  if (!is.null(x$col_order) && length(x$col_order) == p) {
    col_label_list <- c(col_label_list, list(mat1_col_labels[x$col_order]))
  } else {
    col_label_list <- c(col_label_list, list(rep(NA, p)))
  }
}

valid <- !sapply(mat_list, function(x) identical(x, NA))
mat_list       <- mat_list[valid]
heatmap_titles <- heatmap_titles[valid]
col_label_list <- col_label_list[valid]

heatmap_grob <- get_heatmap_grid(mat_list, heatmap_titles, col_label_list,
                                 stage_color_map,
                                 nrow = 2, ncol = 4)


# ============================================================
# 5. Panel B: Spearman vs ARI Scatter Plot
# ============================================================

# --- Pure-R NMI (avoids aricode::sort_pairs compatibility issues) ---
compute_nmi <- function(c1, c2) {
  tab <- table(c1, c2)
  n <- sum(tab)
  p_ij <- tab / n
  p_i <- rowSums(tab) / n
  p_j <- colSums(tab) / n
  mi <- 0
  for (i in seq_len(nrow(tab)))
    for (j in seq_len(ncol(tab)))
      if (tab[i, j] > 0)
        mi <- mi + p_ij[i, j] * log(p_ij[i, j] / (p_i[i] * p_j[j]))
  h1 <- -sum(p_i[p_i > 0] * log(p_i[p_i > 0]))
  h2 <- -sum(p_j[p_j > 0] * log(p_j[p_j > 0]))
  if (h1 + h2 == 0) return(0)
  2 * mi / (h1 + h2)
}

# Exhaustive ARI + Spearman evaluation
evaluate_seriation_rigorous <- function(col_order, true_labels) {
  ordered_labels <- as.character(true_labels[col_order])
  unique_labs    <- unique(as.character(true_labels))
  all_perms      <- combinat::permn(unique_labs)
  label_counts   <- table(true_labels)

  max_acc <- -1
  best_ideal_seq <- NULL

  for (perm in all_perms) {
    ideal_seq   <- rep(perm, times = label_counts[perm])
    current_acc <- sum(ordered_labels == ideal_seq) / length(ordered_labels)
    if (current_acc > max_acc) {
      max_acc        <- current_acc
      best_ideal_seq <- ideal_seq
    }
  }

  best_nmi <- compute_nmi(ordered_labels, best_ideal_seq)
  best_ari <- mclust::adjustedRandIndex(ordered_labels, best_ideal_seq)

  data.frame(NMI = best_nmi, ARI = best_ari)
}

calculate_metrics <- function(method_name, out_list,
                              true_labels, stages) {
  pred_order <- out_list[[method_name]]$col_order
  if (is.null(pred_order)) return(NULL)

  sorted_labels   <- true_labels[pred_order]
  numeric_time    <- as.numeric(factor(sorted_labels, levels = stages))
  spearman_cor    <- abs(cor(seq_along(pred_order), numeric_time,
                            method = "spearman"))

  # Exhaustive max ARI over all label-size permutations
  true_sizes   <- as.numeric(table(true_labels))
  num_clusters <- length(true_sizes)
  all_perms    <- combinat::permn(true_sizes)

  max_ari <- -1
  for (sp in all_perms) {
    pred_clusters <- rep(seq_len(num_clusters), times = sp)
    ari <- adjustedRandIndex(pred_clusters, sorted_labels)
    if (ari > max_ari) max_ari <- ari
  }

  data.frame(Method   = method_name,
             Spearman = spearman_cor,
             ARI      = max_ari)
}

metrics_df <- bind_rows(lapply(method_keys, calculate_metrics,
                               out_list    = allout,
                               true_labels = mat1_labels,
                               stages      = stage_levels))

metrics_clean <- metrics_df %>%
  mutate(Method    = display_name[Method],
         Highlight = ifelse(Method == "BiSer", "Target", "Others"))

med_spearman <- median(metrics_clean$Spearman, na.rm = TRUE)
med_ari      <- median(metrics_clean$ARI,      na.rm = TRUE)

scatter_plot <- ggplot(metrics_clean, aes(x = Spearman, y = ARI)) +
  geom_hline(yintercept = med_ari,      linetype = "dashed",
             color = "gray70", linewidth = 0.5) +
  geom_vline(xintercept = med_spearman, linetype = "dashed",
             color = "gray70", linewidth = 0.5) +
  geom_point(aes(color = Highlight, size = Highlight, shape = Highlight),
             alpha = 0.85) +
  scale_color_manual(values = c("Target" = "#DC0000FF",
                                "Others" = "#7E6148FF")) +
  scale_size_manual(values  = c("Target" = 5, "Others" = 3)) +
  scale_shape_manual(values = c("Target" = 16, "Others" = 17)) +
  geom_text_repel(aes(label = Method, color = Highlight),
                  size = 8, fontface = "bold",
                  box.padding = 0.6, point.padding = 0.4,
                  max.overlaps = 20,
                  segment.color = "gray60", show.legend = FALSE) +
  scale_x_continuous(limits = c(min(metrics_clean$Spearman, na.rm = TRUE) - 0.05,
                                1.05)) +
  scale_y_continuous(limits = c(min(metrics_clean$ARI, na.rm = TRUE) - 0.05,
                                1.05)) +
  theme_classic(base_size = 24) +
  theme(legend.position   = "none",
        axis.title         = element_text(face = "bold", size = 24),
        axis.text          = element_text(color = "black", size = 21),
        panel.grid.major   = element_line(color = "gray95"),
        axis.line          = element_line(linewidth = 0.8)) +
  labs(x = "Trajectory Accuracy (Spearman)",
       y = "Clustering Purity (ARI)")


# ============================================================
# 6. Panel C: Marker Gene Ridge Plot
# ============================================================

biser_ordered_mat <- mat1[allout$biser$row_order, allout$biser$col_order]
num_cells <- ncol(biser_ordered_mat)

# Define three temporal windows
cols_early <- 1:floor(num_cells / 3)
cols_mid   <- (floor(num_cells / 3) + 1):floor(num_cells * 2 / 3)
cols_late  <- (floor(num_cells * 2 / 3) + 1):num_cells

# Compute window means
mean_E <- rowMeans(biser_ordered_mat[, cols_early, drop = FALSE])
mean_M <- rowMeans(biser_ordered_mat[, cols_mid,   drop = FALSE])
mean_L <- rowMeans(biser_ordered_mat[, cols_late,  drop = FALSE])

# Select stage-specific marker genes
gene_early <- rownames(biser_ordered_mat)[which.max(mean_E - pmax(mean_M, mean_L))]
gene_mid   <- rownames(biser_ordered_mat)[which.max(mean_M - pmax(mean_E, mean_L))]
gene_late  <- rownames(biser_ordered_mat)[which.max(mean_L - pmax(mean_E, mean_M))]
target_genes <- c(gene_early, gene_mid, gene_late)
cat("Selected marker genes:", paste(target_genes, collapse = ", "), "\n")

# Map column names to clean stage labels
raw_colnames <- colnames(biser_ordered_mat)
clean_labels <- case_when(
  grepl("Oocyte",     raw_colnames, ignore.case = TRUE) ~ "Oocyte",
  grepl("Zygote",     raw_colnames, ignore.case = TRUE) ~ "Zygote",
  grepl("X2",         raw_colnames, ignore.case = TRUE) ~ "2-cell",
  grepl("X4",         raw_colnames, ignore.case = TRUE) ~ "4-cell",
  grepl("X8",         raw_colnames, ignore.case = TRUE) ~ "8-cell",
  grepl("Morula",     raw_colnames, ignore.case = TRUE) ~ "Morula",
  grepl("blastocyst", raw_colnames, ignore.case = TRUE) ~ "Blastocyst",
  TRUE ~ "Unknown"
)
biser_sorted_labels <- factor(clean_labels, levels = stage_levels_clean)

# Build long-format data for ggridges
marker_df           <- as.data.frame(t(biser_ordered_mat[target_genes, ]))
marker_df$Cell_Stage <- biser_sorted_labels

marker_long <- marker_df %>%
  pivot_longer(cols = all_of(target_genes),
               names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = target_genes))

ridge_plot <- ggplot(marker_long,
                     aes(x = Expression, y = Cell_Stage, fill = Cell_Stage)) +
  geom_density_ridges(alpha = 0.8, scale = 1.5, color = "white",
                      rel_min_height = 0.01) +
  facet_wrap(~ Gene, scales = "free_x") +
  scale_fill_npg() +
  scale_y_discrete(drop = FALSE) +
  theme_classic(base_size = 24) +
  theme(legend.position = "none",
        strip.text      = element_text(face = "bold.italic", size = 24),
        axis.text.y     = element_text(face = "bold", color = "black", size = 21),
        axis.text.x     = element_text(color = "black", size = 21),
        axis.line       = element_line(linewidth = 0.8)) +
  labs(x = "Expression Level", y = "")


# ============================================================
# 7. Panel D: Gene Module Wave Plot
# ============================================================

sorted_mat   <- mat1[allout$biser$row_order, allout$biser$col_order]
num_stages   <- length(stage_levels)
row_splits   <- split(1:nrow(sorted_mat),
                      cut(1:nrow(sorted_mat), num_stages, labels = FALSE))

module_expr <- lapply(row_splits, function(idx) {
  colMeans(sorted_mat[idx, , drop = FALSE])
}) %>%
  bind_rows() %>%
  t() %>%
  as.data.frame()

colnames(module_expr) <- paste0("Module_", seq_len(num_stages))
module_expr$Cell_Order <- seq_len(nrow(module_expr))

module_long <- module_expr %>%
  pivot_longer(starts_with("Module_"),
               names_to = "Gene_Module", values_to = "Mean_Expression")

wave_plot <- ggplot(module_long,
                    aes(x = Cell_Order, y = Mean_Expression,
                        color = Gene_Module)) +
  geom_smooth(method = "loess", span = 0.3, se = FALSE, linewidth = 1.5) +
  geom_vline(xintercept = seq(1, ncol(sorted_mat),
                              length.out = num_stages + 1),
             linetype = "dashed", color = "gray50") +
  scale_color_npg() +
  theme_classic(base_size = 24) +
  theme(legend.position  = "right",
        legend.title      = element_text(face = "bold", size = 22),
        legend.text       = element_text(size = 21),
        axis.title        = element_text(face = "bold", size = 24),
        axis.text         = element_text(color = "black", size = 21),
        axis.line         = element_line(linewidth = 0.8)) +
  labs(x = "Cell Ordering by BiSer", y = "Module Mean Expression",
       color = "Gene Module")


# ============================================================
# 8. Combine All Panels (A / B / C / D)
# ============================================================

panel_A <- ggdraw() +
  draw_grob(heatmap_grob, x = 0, y = 0, width = 1, height = 0.95) +
  draw_label("A", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 35)

panel_B <- ggdraw() +
  draw_plot(scatter_plot, x = 0, y = 0, width = 1, height = 0.95) +
  draw_label("B", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 35)

panel_C <- ggdraw() +
  draw_plot(ridge_plot, x = 0, y = 0, width = 1, height = 0.95) +
  draw_label("C", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 35)

panel_D <- ggdraw() +
  draw_plot(wave_plot, x = 0, y = 0, width = 1, height = 0.95) +
  draw_label("D", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 35)

bottom_row <- plot_grid(panel_B, panel_C, panel_D,
                        nrow = 1, rel_widths = c(0.30, 0.32, 0.38),
                        align = "h")

fig_final <- plot_grid(panel_A, bottom_row,
                       ncol = 1, rel_heights = c(0.60, 0.40))

# ---- Save ----
ggsave("Figure6_Yan_unsupervised.png", fig_final,
       width = 28, height = 24, dpi = 500)
# ---- Print summary ----
cat("\n========================================\n")
cat("Yan Dataset — Unsupervised Feature Selection Summary\n")
cat("========================================\n")
cat("Feature selection: HVG (Highly Variable Genes), top", n_top, "\n")
cat("Matrix dimensions:", nrow(mat1), "genes x", ncol(mat1), "cells\n\n")
cat("--- Metrics ---\n")
print(metrics_clean[, c("Method", "Spearman", "ARI")], row.names = FALSE)
cat("========================================\n")
cat("Figure saved to: Figure6_Yan_unsupervised.png\n")
