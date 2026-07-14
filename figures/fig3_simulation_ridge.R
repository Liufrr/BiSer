################################################################################
# Figure 3: overlapping, low-noise, large Gaussian simulations
################################################################################

source("R/visualization.R")

library(grid)
library(gridExtra)
library(viridis)

data_file <- "output/ALLout_norm_low_1000.RData"
output_dir <- "output"
core_methods <- c("biser", "bs", "tsp_seri", "spec_seri", "Heatmap", "MESBC", "NMF")

load(data_file)
allmetric <- ALLout$allmetric
biser_ari <- vapply(allmetric, function(x) x["biser", "ARI"], numeric(1))
representative <- which.min(abs(biser_ari - median(biser_ari, na.rm = TRUE)))
info <- ALLout$trueinfo[[representative]]
mat_input <- info$mat_input
if (is.null(mat_input)) stop("Regenerate results with the BiSer0608 runner.")

heatmap_mats <- list(info$mat, mat_input)
heatmap_titles <- c("Ground truth", "Disturbed")
for (method in core_methods) {
  out <- ALLout$allout[[representative]][[method]]
  if (is.null(out) || is.null(out$row_order) || is.null(out$col_order)) {
    stop("Missing ordering output for ", method, " in replicate ", representative)
  }
  heatmap_mats[[length(heatmap_mats) + 1L]] <-
    mat_input[out$row_order, out$col_order, drop = FALSE]
  heatmap_titles <- c(heatmap_titles, display_name[method])
}

panel_A <- get_multiple_heatmap_grob(
  heatmap_mats, nrow = 3, ncol = 3,
  title = heatmap_titles, color = viridis(100), legende_title = "Value"
)
panel_B <- make_ridge_plot(
  allmetric, core_methods, method_colors,
  target_metric = "ARI", x_label = "Adjusted Rand Index (ARI)",
  show_y_axis = TRUE
)

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
pdf(file.path(output_dir, "fig3_simulation_ridge.pdf"), width = 16, height = 15)
grid.arrange(panel_A, ggplotGrob(panel_B), nrow = 2, heights = c(2.2, 1))
dev.off()

cat("Figure 3 saved to", file.path(output_dir, "fig3_simulation_ridge.pdf"), "\n")
