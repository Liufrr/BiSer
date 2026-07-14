################################################################################
# Figure 2: mutually exclusive, low-noise, large Gaussian simulations
################################################################################

source("R/visualization.R")

library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
library(dplyr)
library(tidyr)

data_file <- "output/ALLout_norm_exc_low_1000.RData"
output_dir <- "output"
core_methods <- c("biser", "bs", "tsp_seri", "spec_seri", "Heatmap", "MESBC", "NMF")

load(data_file)
allmetric <- ALLout$allmetric

# Use the replicate closest to the median BiSer ARI. This is deterministic and
# avoids selecting the best-performing replicate for the representative panel.
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

metric_data <- bind_rows(lapply(seq_along(allmetric), function(i) {
  data.frame(
    Method = rownames(allmetric[[i]]),
    ARI = allmetric[[i]][, "ARI"],
    F1 = allmetric[[i]][, "F1"],
    Accuracy = allmetric[[i]][, "Accuracy"],
    replicate = i,
    check.names = FALSE
  )
})) %>%
  filter(Method %in% core_methods) %>%
  mutate(Method = factor(display_name[Method], levels = method_order)) %>%
  pivot_longer(c(ARI, F1, Accuracy), names_to = "Metric", values_to = "Score") %>%
  group_by(Method, Metric) %>%
  summarise(Mean = mean(Score, na.rm = TRUE), SD = sd(Score, na.rm = TRUE),
            .groups = "drop")

panel_B <- ggplot(metric_data, aes(Method, Mean, color = Method)) +
  geom_point(size = 3.5) +
  geom_errorbar(
    aes(ymin = pmax(Mean - SD, 0), ymax = pmin(Mean + SD, 1)),
    width = 0.25, linewidth = 0.7
  ) +
  facet_wrap(~Metric, nrow = 1) +
  scale_color_manual(values = method_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Mean score ± 1 SD")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
pdf(file.path(output_dir, "fig2_simulation_heatmap.pdf"), width = 16, height = 14)
grid.arrange(panel_A, ggplotGrob(panel_B), nrow = 2, heights = c(2.2, 1))
dev.off()

cat("Figure 2 saved to", file.path(output_dir, "fig2_simulation_heatmap.pdf"), "\n")
