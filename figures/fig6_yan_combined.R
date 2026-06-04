################################################################################
# Figure 6: Yan Embryo Dataset — Panel B
#
# Panel B: Spearman correlation bar chart
#
# Required data: output/yan_results.RData
################################################################################

source("R/visualization.R")

library(ggplot2)
library(ggrepel)
library(patchwork)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(dplyr)

# ==============================================================================
# Configuration
# ==============================================================================
data_file  <- "output/yan_results.RData"
output_dir <- "output"

# ==============================================================================
# Load
# ==============================================================================
load(data_file)

# ==============================================================================
# Panel B: Spearman correlation bar chart
# ==============================================================================
spearman_results$Method_display <- sapply(spearman_results$Method, function(x) {
  if (x %in% names(display_name)) display_name[x] else x
})
spearman_results$Method_display <- factor(spearman_results$Method_display,
                                          levels = method_order)

panel_B <- ggplot(spearman_results, aes(x = Method_display, y = Rho,
                                         fill = Method_display)) +
  geom_col(alpha = 0.85, width = 0.7) +
  scale_fill_manual(values = method_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pub +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 14),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Spearman Rho")

# ==============================================================================
# Save
# ==============================================================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "fig6_yan_combined.pdf"),
       panel_B, width = 10, height = 7, dpi = 300)

cat("Figure 6 saved to", file.path(output_dir, "fig6_yan_combined.pdf"), "\n")
cat("Note: For the full multi-panel Figure 6, run real_data/analysis_yan.R first.\n")
