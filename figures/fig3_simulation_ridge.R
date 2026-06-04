################################################################################
# Figure 3: Simulation Heatmap Grid + Ridge Plots
#
# Panel A: 3x3 heatmap grid (ground truth, shuffled, and 7 methods)
# Panels B/C/D: Ridge density plots for ARI, NMI, Purity distributions
#
# Required data: output/ALLout_norm_low.RData (or any simulation output)
################################################################################

source("R/visualization.R")

library(ggplot2)
library(ggridges)
library(patchwork)
library(grid)
library(gridExtra)
library(viridis)
library(dplyr)
library(tidyr)
library(tibble)

# ==============================================================================
# Configuration
# ==============================================================================
data_file  <- "output/ALLout_norm_low.RData"
output_dir <- "output"

core_methods <- c("biser", "bs", "tsp_seri", "spec_seri", "Heatmap", "MESBC", "NMF")

# ==============================================================================
# Load data
# ==============================================================================
load(data_file)
allmetric <- ALLout$allmetric

# ==============================================================================
# Panel A: Heatmap grid (3x3)
# ==============================================================================
biser_ari <- sapply(allmetric, function(m) m["biser", "ARI"])
best_iter <- which.max(biser_ari)

true_mat <- ALLout$trueinfo[[best_iter]]$mat
m <- nrow(true_mat); p <- ncol(true_mat)

# Build 9 heatmaps: true, shuffled, then 7 methods
mat_list <- list()
titles <- c("Ground Truth", "Shuffled")
mat_list[[1]] <- true_mat
mat_list[[2]] <- true_mat[sample(m), sample(p)]

for (meth in core_methods) {
  out <- ALLout$allout[[best_iter]][[meth]]
  if (is.null(out) || !is.list(out)) {
    mat_list[[length(mat_list)+1]] <- true_mat
    titles <- c(titles, display_name[meth])
  } else {
    ro <- if (!is.null(out$row_order)) out$row_order else 1:m
    co <- if (!is.null(out$col_order)) out$col_order else 1:p
    mat_list[[length(mat_list)+1]] <- true_mat[ro, co]
    titles <- c(titles, if (meth %in% names(display_name)) display_name[meth] else meth)
  }
}

panel_A <- get_multiple_heatmap_grob(
  mat_list, nrow = 3, ncol = 3,
  title = titles, color = viridis(100),
  legende_title = "Value"
)

# ==============================================================================
# Panels B/C/D: Ridge plots
# ==============================================================================
panel_B <- make_ridge_plot(allmetric, core_methods, method_colors,
                            target_metric = "ARI", x_label = "ARI",
                            show_y_axis = TRUE)

panel_C <- make_ridge_plot(allmetric, core_methods, method_colors,
                            target_metric = "NMI", x_label = "NMI",
                            show_y_axis = FALSE)

panel_D <- make_ridge_plot(allmetric, core_methods, method_colors,
                            target_metric = "purity", x_label = "Purity",
                            show_y_axis = FALSE)

# ==============================================================================
# Combine and save
# ==============================================================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ridge_combined <- panel_B + panel_C + panel_D +
  plot_layout(ncol = 3, widths = c(1.3, 1, 1))

pdf(file.path(output_dir, "fig3_simulation_ridge.pdf"), width = 20, height = 18)
grid.arrange(
  panel_A,
  ggplotGrob(ridge_combined),
  nrow = 2, heights = c(2, 1)
)
dev.off()

cat("Figure 3 saved to", file.path(output_dir, "fig3_simulation_ridge.pdf"), "\n")
