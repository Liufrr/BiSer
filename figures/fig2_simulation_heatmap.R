################################################################################
# Figure 2: Simulation Heatmap + Point-Plot with Error Bars
#
# Panel A: 3x3 heatmap grid (one simulation replicate)
# Panel B: Mean +/- SD point plot for ARI, NMI, Purity across methods
#
# Required data: output/ALLout_norm_low.RData (or any simulation output)
################################################################################

source("R/visualization.R")

library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)
library(viridis)
library(dplyr)
library(tidyr)

# ==============================================================================
# Configuration
# ==============================================================================
data_file  <- "output/ALLout_norm_low.RData"
output_dir <- "output"

# Methods to display and their display order
core_methods <- c("biser", "bs", "tsp_seri", "spec_seri", "Heatmap", "MESBC", "NMF")

# ==============================================================================
# Load data
# ==============================================================================
load(data_file)
allmetric <- ALLout$allmetric

# ==============================================================================
# Panel A: Heatmap grid
# ==============================================================================
# Select a representative simulation replicate (best BiSer ARI)
biser_ari <- sapply(allmetric, function(m) m["biser", "ARI"])
best_iter <- which.max(biser_ari)

# Get original and reordered matrices
true_mat   <- ALLout$trueinfo[[best_iter]]$mat
biser_out  <- ALLout$allout[[best_iter]]$biser
reord_mat  <- biser_out$reordered_mat

# Build heatmap list: true, original (shuffled), and method-reordered
# (Adjust according to available outputs in your data)
heatmap_titles <- c("Ground Truth", "Shuffled", "BiSer")
heatmap_mats   <- list(true_mat, true_mat[sample(nrow(true_mat)),
                                           sample(ncol(true_mat))],
                        reord_mat)

color_pal <- viridis(100)
panel_A <- get_multiple_heatmap_grob(
  heatmap_mats, nrow = 1, ncol = 3,
  title = heatmap_titles, color = color_pal,
  legende_title = "Value"
)

# ==============================================================================
# Panel B: Mean +/- SD point plot (ARI, NMI, Purity)
# ==============================================================================
metrics_long <- bind_rows(lapply(seq_along(allmetric), function(i) {
  df <- as.data.frame(allmetric[[i]])
  df$Method <- rownames(df)
  df$iter <- i
  df
}))

metrics_sub <- metrics_long %>%
  filter(Method %in% core_methods) %>%
  mutate(Method = display_name[Method]) %>%
  mutate(Method = factor(Method, levels = method_order))

metrics_summary <- metrics_sub %>%
  select(Method, ARI, NMI, purity) %>%
  pivot_longer(cols = c(ARI, NMI, purity), names_to = "Metric", values_to = "Value") %>%
  group_by(Method, Metric) %>%
  summarise(Mean = mean(Value, na.rm = TRUE),
            SD   = sd(Value, na.rm = TRUE), .groups = "drop")

panel_B <- ggplot(metrics_summary, aes(x = Method, y = Mean, color = Method)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = pmin(Mean + SD, 1)),
                width = 0.3, linewidth = 0.8) +
  facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = method_colors) +
  theme_pub +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 13),
    strip.text   = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Score")

# ==============================================================================
# Combine and save
# ==============================================================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

pdf(file.path(output_dir, "fig2_simulation_heatmap.pdf"), width = 16, height = 12)
grid.arrange(
  panel_A, ggplotGrob(panel_B),
  nrow = 2, heights = c(1, 1)
)
dev.off()

cat("Figure 2 saved to", file.path(output_dir, "fig2_simulation_heatmap.pdf"), "\n")
