################################################################################
# Figure 5: Main Boxplot Comparison (4 Scenarios)
#
# 2x2 panel layout showing ARI distributions for:
#   Exclusive/Low, Exclusive/High, Overlapping/Low, Overlapping/High
#
# Required data:
#   output/ALLout_norm_exc_low.RData, output/ALLout_norm_exc_high.RData
#   output/ALLout_norm_low.RData, output/ALLout_norm_high.RData
################################################################################

source("R/visualization.R")

library(ggplot2)
library(dplyr)
library(patchwork)

# ==============================================================================
# Configuration
# ==============================================================================
output_dir <- "output"

core_methods <- c("biser", "bs", "tsp_seri", "spec_seri", "Heatmap", "MESBC", "NMF")

scenarios <- list(
  list(file = "output/ALLout_norm_exc_low.RData",  title = "Exclusive / Low noise"),
  list(file = "output/ALLout_norm_exc_high.RData", title = "Exclusive / High noise"),
  list(file = "output/ALLout_norm_low.RData",      title = "Overlapping / Low noise"),
  list(file = "output/ALLout_norm_high.RData",     title = "Overlapping / High noise")
)

# ==============================================================================
# Extract ARI from each scenario
# ==============================================================================
process_one_scenario <- function(file_path, scenario_title) {
  load(file_path)
  allmetric <- ALLout$allmetric

  bind_rows(lapply(seq_along(allmetric), function(i) {
    m <- allmetric[[i]]
    data.frame(
      Method   = rownames(m),
      ARI      = m[, "ARI"],
      Scenario = scenario_title,
      stringsAsFactors = FALSE
    )
  })) %>%
    filter(Method %in% core_methods) %>%
    mutate(Method = if_else(Method %in% names(display_name),
                            display_name[Method], Method))
}

all_data <- bind_rows(lapply(scenarios, function(s) {
  process_one_scenario(s$file, s$title)
}))

all_data$Method <- factor(all_data$Method, levels = method_order)
all_data$Scenario <- factor(all_data$Scenario, levels = sapply(scenarios, `[[`, "title"))

# ==============================================================================
# 2x2 boxplot panel
# ==============================================================================
fig5 <- ggplot(all_data, aes(x = Method, y = ARI, fill = Method)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.8, width = 0.7) +
  facet_wrap(~ Scenario, scales = "free_y", nrow = 2) +
  scale_fill_manual(values = method_colors) +
  theme_pub +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 13),
    strip.text      = element_text(size = 15, face = "bold"),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Adjusted Rand Index (ARI)")

# ==============================================================================
# Save
# ==============================================================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "fig5_boxplot_main.pdf"),
       fig5, width = 14, height = 10, dpi = 300)

cat("Figure 5 saved to", file.path(output_dir, "fig5_boxplot_main.pdf"), "\n")
