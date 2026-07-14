################################################################################
# Li Dataset Analysis — Unsupervised Feature Selection (HVG)
#
# Single-cell RNA-seq data (Li et al.), 7 cell lines:
#   A549, GM12878, H1437, HCT116, IMR90, H1, K562
# Uses unsupervised HVG selection (scran::modelGeneVar), then runs all methods.
#
# Produces:
#   - Heatmap comparisons (2x4 grid with group annotation bar)
#   - ARI vs NMI scatter plot
#   - Combined Figure_Li.png
################################################################################

library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggsci)
library(cowplot)
library(mclust)
library(combinat)
library(scran)
library(igraph)
library(NMF)

# ============================================================
# 0. Source BiSer methods & unsupervised feature selection
# ============================================================
source("R/methods.R")
source("R/feature_selection_unsupervised.R")

# ============================================================
# 1. Global Settings
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

li_colors <- c("A549"    = "#E64B35FF", "GM12878" = "#4DBBD5FF",
               "H1437"   = "#00A087FF", "HCT116"  = "#3C5488FF",
               "IMR90"   = "#F39B7FFF", "H1"      = "#8491B4FF",
               "K562"    = "#B09C85FF")

color_palette <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Pure-R NMI (avoids aricode::sort_pairs compatibility issues)
compute_nmi <- function(c1, c2) {
  tab  <- table(c1, c2)
  n    <- sum(tab)
  p_ij <- tab / n
  p_i  <- rowSums(tab) / n
  p_j  <- colSums(tab) / n
  mi   <- 0
  for (i in seq_len(nrow(tab)))
    for (j in seq_len(ncol(tab)))
      if (tab[i, j] > 0)
        mi <- mi + p_ij[i, j] * log(p_ij[i, j] / (p_i[i] * p_j[j]))
  h1 <- -sum(p_i[p_i > 0] * log(p_i[p_i > 0]))
  h2 <- -sum(p_j[p_j > 0] * log(p_j[p_j > 0]))
  if (h1 + h2 == 0) return(0)
  2 * mi / (h1 + h2)
}


# ============================================================
# 2. Load Data & HVG Feature Selection
# ============================================================

load("data/li.RData")

raw_matrix   <- li$data                # genes x cells (sparse or dense)
group_labels <- as.factor(li[["annotation"]][["type"]])

cat("Raw matrix:", nrow(raw_matrix), "genes x", ncol(raw_matrix), "cells\n")
cat("Group labels:", paste(levels(group_labels), collapse = ", "), "\n")

# Unsupervised HVG selection (scran, appropriate for scRNA-seq)
n_top <- 100
feat  <- get_clustering_matrix(raw_matrix, n_top = n_top, method = "variance", norm = TRUE)
top_mat <- feat$top_mat

# Z-score standardization (per gene)
final_genes <- as.matrix(top_mat)
class(final_genes) <- "numeric"
data_scaled <- t(scale(t(final_genes)))
data_scaled[is.na(data_scaled)] <- 0

mat0 <- data_scaled
cat("After HVG selection:", nrow(mat0), "genes x", ncol(mat0), "cells\n")


# ============================================================
# 3. Shuffle & Run All Methods
# ============================================================

set.seed(42)
mat1 <- mat0[sample(1:nrow(mat0), nrow(mat0)),
             sample(1:ncol(mat0), ncol(mat0))]

allout <- list()

# BiSer
allout$biser <- biser(mat1, noise = TRUE, simmeth = "t")

# Spectral seriation
allout$spec_seri <- spec_seri(mat1)

# TSP seriation
allout$tsp_seri <- run_tsp_seriation(mat1)

# R seriation (Heatmap)
allout$Heatmap <- Rseriation(mat1, "Heatmap")

# Yang BS (BiSer_SVD)
allout$bs <- Yang_bs(mat1 - min(mat1), "bs")

# MESBC
allout$MESBC <- bicluster(mat1, "MESBC")

# NMF
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

cat("All methods completed.\n")


# ============================================================
# 4. Build mat1 using BiSer sim order & derive labels
# ============================================================

aa   <- allout[["biser"]][["sim"]]
mat1 <- mat0[rownames(aa)[1:nrow(mat0)],
             rownames(aa)[(nrow(mat0) + 1):nrow(aa)]]

# Derive cell labels for mat1 columns
col_shuffle_idx <- match(colnames(mat1), colnames(mat0))
mat1_labels     <- as.character(group_labels[col_shuffle_idx])

m <- nrow(mat1)
p <- ncol(mat1)

cat("mat1 dimensions:", m, "genes x", p, "cells\n")
cat("Label distribution:\n")
print(table(mat1_labels))


# ============================================================
# 5. Build Heatmap Matrices (one per method + Original)
# ============================================================

# Ordered matrices
mat_list       <- list("Original" = mat1)
col_label_list <- list(mat1_labels)

for (mk in method_keys) {
  x <- allout[[mk]]
  if (!is.null(x$row_order) && length(x$row_order) == m && length(x$col_order) == p) {
    mat_list[[mk]] <- mat1[x$row_order, x$col_order]
    col_label_list[[length(col_label_list) + 1]] <- mat1_labels[x$col_order]
  } else {
    mat_list[[mk]] <- NA
    col_label_list[[length(col_label_list) + 1]] <- rep(NA, p)
  }
}

heatmap_titles <- c("Original", display_name[method_keys])

# Remove failed methods
valid <- !sapply(mat_list, function(x) identical(x, NA))
mat_list       <- mat_list[valid]
heatmap_titles <- heatmap_titles[valid]
col_label_list <- col_label_list[valid]


# ============================================================
# 6. Draw Heatmap Grid (2x4) with Legends
# ============================================================

# Shared color scale
all_vals <- unlist(mat_list)
max_abs  <- max(abs(all_vals), na.rm = TRUE)
col_fun  <- colorRamp2(seq(-max_abs, max_abs, length.out = length(color_palette)),
                       color_palette)

# Expression legend
expr_legend <- Legend(
  col_fun       = col_fun,
  title         = "Expression",
  title_gp      = gpar(fontsize = 26, fontface = "bold", fontfamily = "sans"),
  labels_gp     = gpar(fontsize = 24, fontfamily = "sans"),
  legend_height = unit(7, "cm"),
  grid_width    = unit(1, "cm"),
  title_gap     = unit(0.6, "cm")
)

# Group legend
group_legend <- Legend(
  labels    = names(li_colors),
  legend_gp = gpar(fill = li_colors),
  title     = "Group",
  title_gp  = gpar(fontsize = 26, fontface = "bold", fontfamily = "sans"),
  labels_gp = gpar(fontsize = 23, fontfamily = "sans"),
  grid_width  = unit(0.8, "cm"),
  grid_height = unit(0.8, "cm"),
  row_gap     = unit(0.5, "cm"),
  title_gap   = unit(0.6, "cm")
)

legend_grob <- grid.grabExpr({
  draw(packLegend(expr_legend, group_legend,
                  direction = "vertical", gap = unit(2, "cm")))
})

# Individual heatmap grobs
ht_grobs <- lapply(seq_along(mat_list), function(i) {
  top_anno <- HeatmapAnnotation(
    Group = col_label_list[[i]],
    col   = list(Group = li_colors),
    show_legend          = FALSE,
    show_annotation_name = FALSE,
    annotation_height    = unit(0.35, "cm"),
    simple_anno_size     = unit(0.35, "cm")
  )
  ht <- Heatmap(
    mat_list[[i]],
    col              = col_fun,
    top_annotation   = top_anno,
    cluster_rows     = FALSE,
    cluster_columns  = FALSE,
    show_row_names   = FALSE,
    show_column_names = FALSE,
    column_title     = heatmap_titles[i],
    column_title_gp  = gpar(fontsize = 28, fontface = "bold", fontfamily = "sans"),
    heatmap_legend_param = list(title = NULL),
    rect_gp = gpar(col = NA, lwd = 0),
    name    = paste0("ht_", i)
  )
  grid.grabExpr(draw(ht, show_heatmap_legend = FALSE))
})

# Assemble 2x4 grid + legends
heatmap_grob <- arrangeGrob(
  arrangeGrob(grobs = ht_grobs, nrow = 2, ncol = 4),
  legend_grob,
  nrow = 1, widths = c(4, 1.2)
)


# ============================================================
# 7. Compute ARI & NMI (Exhaustive Permutation Matching)
# ============================================================

metrics_list <- list()

unique_labs  <- unique(mat1_labels)
all_perms    <- combinat::permn(unique_labs)
label_counts <- table(mat1_labels)

for (mk in method_keys) {
  x <- allout[[mk]]
  if (is.null(x$col_order)) next

  ordered_labels <- mat1_labels[x$col_order]

  max_acc        <- -1
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

  metrics_list[[mk]] <- data.frame(Method = mk, NMI = best_nmi, ARI = best_ari)
}

metrics_df <- bind_rows(metrics_list)

# Print results
metrics_print <- metrics_df
metrics_print$Method <- display_name[metrics_print$Method]
cat("\n========== Li ARI & NMI ==========\n")
print(metrics_print[order(-metrics_print$ARI), c("Method", "ARI", "NMI")],
      row.names = FALSE)
cat("\n")


# ============================================================
# 8. Scatter Plot (ARI vs NMI)
# ============================================================

metrics_clean <- metrics_df %>%
  mutate(Method    = display_name[Method],
         Highlight = ifelse(Method == "BiSer", "Target", "Others"))

scatter_plot <- ggplot(metrics_clean, aes(x = ARI, y = NMI)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray80") +
  geom_point(aes(color = Highlight, size = Highlight, shape = Highlight),
             alpha = 0.85) +
  scale_color_manual(values = c("Target" = "#DC0000FF", "Others" = "#7E6148FF")) +
  scale_size_manual(values  = c("Target" = 6, "Others" = 4)) +
  scale_shape_manual(values = c("Target" = 16, "Others" = 17)) +
  geom_text_repel(aes(label = Method, color = Highlight),
                  size = 8, fontface = "bold",
                  box.padding = 0.6, point.padding = 0.4,
                  max.overlaps = 20,
                  segment.color = "gray60", show.legend = FALSE) +
  scale_x_continuous(limits = c(min(metrics_clean$ARI, na.rm = TRUE) - 0.05, 1.05)) +
  scale_y_continuous(limits = c(min(metrics_clean$NMI, na.rm = TRUE) - 0.05, 1.05)) +
  theme_classic(base_size = 24) +
  theme(legend.position   = "none",
        axis.title         = element_text(face = "bold", size = 24),
        axis.text          = element_text(color = "black", size = 21),
        panel.grid.major   = element_line(color = "gray95"),
        axis.line          = element_line(linewidth = 0.8)) +
  labs(x = "Adjusted Rand Index (ARI)",
       y = "Normalized Mutual Information (NMI)")


# ============================================================
# 9. Combine Panels & Save
# ============================================================

panel_A <- ggdraw() +
  draw_grob(heatmap_grob, x = 0, y = 0, width = 1, height = 0.95) +
  draw_label("A", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 35)

panel_B <- ggdraw() +
  draw_plot(scatter_plot, x = 0, y = 0, width = 1, height = 0.95) +
  draw_label("B", x = 0.01, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 35)

fig_li <- plot_grid(panel_A, panel_B,
                    nrow = 1, rel_widths = c(0.60, 0.40))

ggsave("Figure_Li.png", fig_li,
       width = 28, height = 12, dpi = 500)
save(allout, file = "allout_li.RData")

# ---- Print summary ----
cat("\n========================================\n")
cat("Li Dataset — Unsupervised Feature Selection Summary\n")
cat("========================================\n")
cat("Feature selection: HVG (scran modelGeneVar), top", n_top, "\n")
cat("Matrix dimensions:", nrow(mat1), "genes x", ncol(mat1), "cells\n\n")
cat("--- Metrics ---\n")
print(metrics_clean[, c("Method", "ARI", "NMI")], row.names = FALSE)
cat("========================================\n")
cat("Figure saved to: Figure_Li.png\n")
