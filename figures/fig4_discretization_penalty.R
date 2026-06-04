################################################################################
# Figure 4: Discretization Penalty (Dumbbell / Shift Plots)
#
# Shows the performance penalty when increasing noise, comparing methods.
# Uses paired (low vs high) simulation results.
#
# Required data:
#   output/ALLout_norm_low.RData, output/ALLout_norm_high.RData
################################################################################

source("R/visualization.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# ==============================================================================
# Configuration
# ==============================================================================
data_low   <- "output/ALLout_norm_low.RData"
data_high  <- "output/ALLout_norm_high.RData"
output_dir <- "output"

core_methods <- c("biser", "bs", "tsp_seri", "spec_seri", "Heatmap", "MESBC", "NMF")
target_metrics <- c("ARI", "NMI", "purity")

# ==============================================================================
# Load and aggregate
# ==============================================================================
load(data_low);  allmetric_low  <- ALLout$allmetric
load(data_high); allmetric_high <- ALLout$allmetric

extract_means <- function(allmetric, methods, metrics) {
  bind_rows(lapply(allmetric, function(m) {
    df <- as.data.frame(m)
    df$Method <- rownames(df)
    df
  })) %>%
    filter(Method %in% methods) %>%
    select(Method, all_of(metrics)) %>%
    group_by(Method) %>%
    summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
}

mean_low  <- extract_means(allmetric_low,  core_methods, target_metrics)
mean_high <- extract_means(allmetric_high, core_methods, target_metrics)

# Combine into paired data for dumbbell plot
combined <- inner_join(
  mean_low  %>% pivot_longer(-Method, names_to = "Metric", values_to = "Low"),
  mean_high %>% pivot_longer(-Method, names_to = "Metric", values_to = "High"),
  by = c("Method", "Metric")
) %>%
  mutate(
    Method = if_else(Method %in% names(display_name), display_name[Method], Method),
    Method = factor(Method, levels = method_order),
    Delta  = High - Low
  )

# ==============================================================================
# Dumbbell plot
# ==============================================================================
dumbbell <- ggplot(combined, aes(y = Method)) +
  geom_segment(aes(x = High, xend = Low, yend = Method, color = Method),
               linewidth = 1.2, alpha = 0.7) +
  geom_point(aes(x = Low, color = Method), size = 4, shape = 16) +
  geom_point(aes(x = High, color = Method), size = 4, shape = 17) +
  facet_wrap(~ Metric, scales = "free_x", nrow = 1) +
  scale_color_manual(values = method_colors) +
  theme_pub +
  theme(
    strip.text      = element_text(size = 16, face = "bold"),
    legend.position = "none",
    axis.text.y     = element_text(size = 14)
  ) +
  labs(x = "Score (circle = Low noise, triangle = High noise)", y = NULL)

# ==============================================================================
# Save
# ==============================================================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "fig4_discretization_penalty.pdf"),
       dumbbell, width = 14, height = 6, dpi = 300)

cat("Figure 4 saved to", file.path(output_dir, "fig4_discretization_penalty.pdf"), "\n")
