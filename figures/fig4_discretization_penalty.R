################################################################################
# Figure 4: ARI before and after boundary detection in large Gaussian matrices
################################################################################

source("R/visualization.R")

library(ggplot2)
library(dplyr)

output_dir <- "output"
ordering_methods <- c("biser", "bs", "tsp_seri", "spec_seri", "Heatmap")
scenarios <- data.frame(
  file = file.path("output", c(
    "ALLout_norm_exc_low_1000.RData",
    "ALLout_norm_exc_high_1000.RData",
    "ALLout_norm_low_1000.RData",
    "ALLout_norm_high_1000.RData"
  )),
  Scenario = c(
    "A  Exclusive / low noise",
    "B  Exclusive / high noise",
    "C  Overlapping / low noise",
    "D  Overlapping / high noise"
  ),
  stringsAsFactors = FALSE
)

extract_scenario <- function(path, label) {
  if (!file.exists(path)) stop("Missing simulation output: ", path)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  bind_rows(lapply(env$ALLout$allmetric, function(metric) {
    data.frame(
      Method = rownames(metric),
      ARI_pre = metric[, "ARI_pre"],
      ARI_post = metric[, "ARI_post"],
      stringsAsFactors = FALSE
    )
  })) %>%
    filter(Method %in% ordering_methods) %>%
    group_by(Method) %>%
    summarise(
      ARI_pre = mean(ARI_pre, na.rm = TRUE),
      ARI_post = mean(ARI_post, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Scenario = label, Penalty = ARI_pre - ARI_post)
}

plot_data <- bind_rows(Map(extract_scenario, scenarios$file, scenarios$Scenario)) %>%
  mutate(
    Method = factor(display_name[Method], levels = method_order[1:5]),
    Scenario = factor(Scenario, levels = scenarios$Scenario)
  )

fig4 <- ggplot(plot_data, aes(x = Method, color = Method)) +
  geom_segment(aes(y = ARI_post, yend = ARI_pre, xend = Method),
               linewidth = 1, alpha = 0.75) +
  geom_point(aes(y = ARI_pre), shape = 21, fill = "white", size = 3.4,
             stroke = 1) +
  geom_point(aes(y = ARI_post), shape = 16, size = 3.4) +
  geom_text(aes(y = pmin(ARI_pre, ARI_post) - 0.035,
                label = sprintf("%.2f", Penalty)),
            vjust = 1, size = 3.4, color = "black") +
  facet_wrap(~Scenario, nrow = 2) +
  scale_color_manual(values = method_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pub +
  theme(
    strip.text = element_text(face = "bold", size = 13),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "ARI (open: pre-boundary; filled: post-boundary)")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(file.path(output_dir, "fig4_discretization_penalty.pdf"),
       fig4, width = 13, height = 9, dpi = 300)
cat("Figure 4 saved to", file.path(output_dir, "fig4_discretization_penalty.pdf"), "\n")
