################################################################################
# Marker Candidate Gene Selection (based on BiSer ordering)
#
# Original approach retained: use BiSer-ordered tercile splits to define
# early/mid/late windows.
# Improvement: return top-5 candidates for manual selection of genes with
# the best literature support.
#
# Prerequisite: analysis_yan_unsupervised.R has been run; allout, mat1, etc.
#               are available in the environment.
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(ggsci)

stage_levels_clean <- c("Oocyte", "Zygote", "2-cell", "4-cell",
                        "8-cell", "Morula", "Blastocyst")

# BiSer-ordered matrix
biser_ordered_mat <- mat1[allout$biser$row_order, allout$biser$col_order]
num_cells <- ncol(biser_ordered_mat)

# BiSer ordering tercile split
cols_early <- 1:floor(num_cells / 3)
cols_mid   <- (floor(num_cells / 3) + 1):floor(num_cells * 2 / 3)
cols_late  <- (floor(num_cells * 2 / 3) + 1):num_cells

mean_E <- rowMeans(biser_ordered_mat[, cols_early, drop = FALSE])
mean_M <- rowMeans(biser_ordered_mat[, cols_mid,   drop = FALSE])
mean_L <- rowMeans(biser_ordered_mat[, cols_late,  drop = FALSE])

score_early <- mean_E - pmax(mean_M, mean_L)
score_mid   <- mean_M - pmax(mean_E, mean_L)
score_late  <- mean_L - pmax(mean_E, mean_M)

# ============================================================
# Print top-5 candidates
# ============================================================
top_k <- 5

cat("========================================\n")
cat("BiSer Ordering Tercile Marker Candidates\n")
cat("========================================\n\n")

cat("--- Early (BiSer first 1/3) top-5 ---\n")
for (g in names(sort(score_early, decreasing = TRUE))[1:top_k]) {
  cat(sprintf("  %-12s | Score=%+.3f | E=%+.2f  M=%+.2f  L=%+.2f\n",
              g, score_early[g], mean_E[g], mean_M[g], mean_L[g]))
}

cat("\n--- Mid (BiSer middle 1/3) top-5 ---\n")
for (g in names(sort(score_mid, decreasing = TRUE))[1:top_k]) {
  cat(sprintf("  %-12s | Score=%+.3f | E=%+.2f  M=%+.2f  L=%+.2f\n",
              g, score_mid[g], mean_E[g], mean_M[g], mean_L[g]))
}

cat("\n--- Late (BiSer last 1/3) top-5 ---\n")
for (g in names(sort(score_late, decreasing = TRUE))[1:top_k]) {
  cat(sprintf("  %-12s | Score=%+.3f | E=%+.2f  M=%+.2f  L=%+.2f\n",
              g, score_late[g], mean_E[g], mean_M[g], mean_L[g]))
}

cat("\nPlease select the gene with the best literature support from each window.\n")
cat("After selection, modify target_genes below to generate the ridge plot.\n\n")

# ============================================================
# Modify here: replace with your selected genes
# ============================================================
gene_early <- names(sort(score_early, decreasing = TRUE))[1]  # default top-1
gene_mid   <- names(sort(score_mid,   decreasing = TRUE))[1]  # <- replaceable
gene_late  <- names(sort(score_late,  decreasing = TRUE))[1]  # default top-1
target_genes <- c(gene_early, gene_mid, gene_late)

cat("Current selection:", paste(target_genes, collapse = ", "), "\n\n")

# ============================================================
# Ridge Plot
# ============================================================
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
cell_stages <- factor(clean_labels, levels = stage_levels_clean)

marker_df <- as.data.frame(t(biser_ordered_mat[target_genes, ]))
marker_df$Cell_Stage <- cell_stages

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
  theme_classic(base_size = 20) +
  theme(legend.position = "none",
        strip.text  = element_text(face = "bold.italic", size = 20),
        axis.text.y = element_text(face = "bold", color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 14),
        axis.line   = element_line(linewidth = 0.6)) +
  labs(x = "Expression Level", y = "")

ggsave("marker_candidates.png", ridge_plot,
       width = 16, height = 8, dpi = 300)

cat("Ridge plot saved: marker_candidates.png\n")
