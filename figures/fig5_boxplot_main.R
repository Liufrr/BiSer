################################################################################
# Figure 5: ablation analysis and distributional robustness
################################################################################

source("R/visualization.R")

library(ggplot2)
library(dplyr)
library(patchwork)

output_dir <- "output"
core_methods <- c("biser", "bs", "tsp_seri", "spec_seri", "Heatmap", "MESBC", "NMF")

make_design <- function(prefix, distribution) {
  expand.grid(
    Structure = c("Exclusive", "Overlapping"),
    Noise = c("Low noise", "High noise"),
    Scale = c("250×150", "1000×650"),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Distribution = distribution,
      file = file.path(
        "output",
        paste0(
          "ALLout_", prefix,
          ifelse(Structure == "Exclusive", "_exc", ""),
          ifelse(Noise == "Low noise", "_low", "_high"),
          ifelse(Scale == "1000×650", "_1000", ""),
          ".RData"
        )
      ),
      StructureScale = paste(Structure, Scale, sep = "\n")
    )
}

read_ari <- function(file, distribution, structure, noise, scale, structure_scale) {
  if (!file.exists(file)) stop("Missing simulation output: ", file)
  env <- new.env(parent = emptyenv())
  load(file, envir = env)
  bind_rows(lapply(seq_along(env$ALLout$allmetric), function(i) {
    metric <- env$ALLout$allmetric[[i]]
    data.frame(
      Method = rownames(metric), ARI = metric[, "ARI"], replicate = i,
      stringsAsFactors = FALSE
    )
  })) %>%
    mutate(
      Distribution = distribution, Structure = structure,
      Noise = noise, Scale = scale, StructureScale = structure_scale
    )
}

collect_design <- function(design) {
  bind_rows(Map(
    read_ari, design$file, design$Distribution, design$Structure,
    design$Noise, design$Scale, design$StructureScale
  ))
}

gaussian <- collect_design(make_design("norm", "Gaussian")) %>%
  filter(Method %in% c("biser", "bs", "tsp_seri")) %>%
  mutate(
    Method = factor(display_name[Method],
                    levels = c("BiSer", "BiSer_SVD", "BiSer_TSP")),
    Noise = factor(Noise, levels = c("Low noise", "High noise")),
    StructureScale = factor(
      StructureScale,
      levels = c("Exclusive\n250×150", "Exclusive\n1000×650",
                 "Overlapping\n250×150", "Overlapping\n1000×650")
    )
  )

panel_A <- ggplot(gaussian, aes(StructureScale, ARI, fill = Method)) +
  geom_boxplot(width = 0.72, outlier.size = 0.25, position = position_dodge(.8)) +
  facet_wrap(~Noise, nrow = 1) +
  scale_fill_manual(values = method_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 10),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  ) +
  labs(x = NULL, y = "ARI", fill = NULL, tag = "A")

count_data <- bind_rows(
  collect_design(make_design("poisson", "Poisson")),
  collect_design(make_design("NBP", "Negative binomial"))
) %>%
  filter(Method %in% core_methods) %>%
  mutate(
    Method = factor(display_name[Method], levels = method_order),
    Distribution = factor(Distribution,
                          levels = c("Poisson", "Negative binomial"))
  )

panel_B <- ggplot(count_data, aes(Method, ARI, fill = Method)) +
  geom_boxplot(width = 0.7, outlier.size = 0.2) +
  facet_wrap(~Distribution, nrow = 1) +
  scale_fill_manual(values = method_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_pub +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "ARI", tag = "B")

fig5 <- panel_A / panel_B + plot_layout(heights = c(1, 1))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
ggsave(file.path(output_dir, "fig5_boxplot_main.pdf"),
       fig5, width = 15, height = 11, dpi = 300)
cat("Figure 5 saved to", file.path(output_dir, "fig5_boxplot_main.pdf"), "\n")
