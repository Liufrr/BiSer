################################################################################
# Scalability Figure: Combined A + B + C
#
# Panel A: BiSer step-wise log-log scaling (SVD, Similarity, TSP, Other)
# Panel B: All methods cross-comparison log-log (seriation vs biclustering)
# Panel C: MESBC / NMF per-K timing (faceted, colored by matrix size)
#
# Layout: (A | B) / C
#
# Required data:
#   output/stepwise_timings.csv
#   output/method_comparison.csv
#   output/biclustering_per_k.csv
################################################################################

source("R/visualization.R")

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)
library(forcats)

# ==============================================================================
# Configuration
# ==============================================================================
output_dir <- "output"

# --- Load data ---
stepwise_results <- read.csv(file.path(output_dir, "stepwise_timings.csv"),
                             stringsAsFactors = FALSE)
comp_results     <- read.csv(file.path(output_dir, "method_comparison.csv"),
                             stringsAsFactors = FALSE)
perk_results     <- read.csv(file.path(output_dir, "biclustering_per_k.csv"),
                             stringsAsFactors = FALSE)

# --- Color palettes ---
# Panel A: step colors
step_order <- c("SVD decomposition", "Similarity construction",
                "TSP solving", "Other")
step_group_colors <- c(
  "SVD decomposition"       = "#F28E2B",
  "Similarity construction" = "#76B7B2",
  "TSP solving"             = "#EDC948",
  "Other"                   = "#BAB0AC"
)

# Panel B: method colors (seriation = solid, biclustering = dashed)
method_colors_fig <- c(
  "BiSer"     = "#E15759",
  "BiSer_TSP" = "#4E79A7",
  "BiSer_SVD" = "#8FAACC",
  "Spectral"  = "#A0D8B0",
  "Heatmap"   = "#3E9E6E",
  "MESBC"     = "#F28E2B",
  "NMF"       = "#FFBE7A"
)

method_linetypes <- c(
  "BiSer" = "solid", "BiSer_SVD" = "solid", "BiSer_TSP" = "solid",
  "Spectral" = "solid", "Heatmap" = "solid",
  "MESBC" = "dashed", "NMF" = "dashed"
)

method_order_fig <- c("BiSer", "BiSer_SVD", "BiSer_TSP", "MESBC", "NMF",
                      "Spectral", "Heatmap")

# Panel C: matrix size colors
n_size_colors <- c(
  "n = 300"  = "#4E79A7",
  "n = 800"  = "#E15759",
  "n = 1000" = "#EDC948"
)

# Log10 label formatter
log_labels <- function(x) sprintf("%g", x)

# ==============================================================================
# Panel A: BiSer step-wise log-log scaling
# ==============================================================================
p2_data <- stepwise_results %>%
  filter(step != "total") %>%
  mutate(step_group = case_when(
    step %in% c("tsp_solving")  ~ "TSP solving",
    step %in% c("similarity")   ~ "Similarity construction",
    step %in% c("svd")          ~ "SVD decomposition",
    TRUE                         ~ "Other"
  )) %>%
  group_by(n, rep, step_group) %>%
  summarise(time_sec = sum(time_sec), .groups = "drop") %>%
  group_by(n, step_group) %>%
  summarise(mean_time = mean(time_sec), sd_time = sd(time_sec), .groups = "drop") %>%
  filter(mean_time > 0) %>%
  mutate(step_group = factor(step_group, levels = step_order))

panel_A <- ggplot(p2_data, aes(x = n, y = mean_time, color = step_group)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = pmax(mean_time - sd_time, mean_time * 0.1),
                    ymax = mean_time + sd_time), width = 0.05, alpha = 0.4) +
  scale_x_log10(labels = log_labels) +
  scale_y_log10(labels = log_labels) +
  scale_color_manual(values = step_group_colors, name = "Step",
                     breaks = step_order) +
  labs(x = "n = m + p", y = "Time (s)") +
  theme_pub

# ==============================================================================
# Panel B: All methods cross-comparison (log-log)
# ==============================================================================
p3_data <- comp_results %>%
  filter(!is.na(time_sec)) %>%
  group_by(n, method, method_type) %>%
  summarise(mean_time = mean(time_sec),
            sd_time   = sd(time_sec), .groups = "drop") %>%
  mutate(sd_time = replace_na(sd_time, 0)) %>%
  filter(mean_time > 0)

# Rename methods to display names
p3_data <- p3_data %>%
  mutate(method = fct_recode(method,
    BiSer_SVD = "BS",
    BiSer_TSP = "TSP"
  )) %>%
  mutate(method = factor(method, levels = method_order_fig))

floor_val <- min(p3_data$mean_time, na.rm = TRUE) * 0.1
p3_data <- p3_data %>%
  mutate(ymin_val = pmax(mean_time - sd_time, floor_val),
         ymax_val = mean_time + sd_time) %>%
  na.omit()

panel_B <- ggplot(p3_data, aes(x = n, y = mean_time,
                                color = method, linetype = method)) +
  geom_line(linewidth = 0.9) +
  geom_point(aes(shape = method_type), size = 2.5) +
  geom_ribbon(aes(ymin = ymin_val, ymax = ymax_val, fill = method),
              alpha = 0.06, color = NA) +
  scale_x_log10(labels = log_labels) +
  scale_y_log10(labels = log_labels) +
  scale_color_manual(values = method_colors_fig, name = "Method",
                     breaks = method_order_fig) +
  scale_fill_manual(values = method_colors_fig, guide = "none") +
  scale_linetype_manual(values = method_linetypes, guide = "none") +
  scale_shape_manual(values = c("seriation" = 16, "biclustering" = 17),
                     name = "Type") +
  guides(color = guide_legend(order = 1, nrow = 2, byrow = TRUE),
         shape = guide_legend(order = 2)) +
  labs(x = "n = m + p", y = "Time (s)") +
  theme_pub +
  theme(legend.box       = "vertical",
        legend.box.just  = "left",
        legend.spacing.y = unit(0.1, "cm"))

# ==============================================================================
# Panel C: MESBC / NMF per-K timing (faceted)
# ==============================================================================
if (nrow(perk_results) > 0) {
  # Select representative matrix sizes
  all_n <- sort(unique(perk_results$n))
  selected_n <- c(300, 800, 1000)
  # Fall back to available sizes if not all present
  selected_n <- intersect(selected_n, all_n)
  if (length(selected_n) < 3 && length(all_n) >= 3) {
    selected_n <- all_n[c(1, ceiling(length(all_n)/2), length(all_n))]
  }

  n_labels_ordered <- paste0("n = ", sort(selected_n))
  n_colors <- c("#4E79A7", "#E15759", "#EDC948")[1:length(selected_n)]
  names(n_colors) <- n_labels_ordered

  p7_data <- perk_results %>%
    filter(n %in% selected_n) %>%
    group_by(n, method, K) %>%
    summarise(mean_time = mean(time_sec), .groups = "drop") %>%
    mutate(n_label = factor(paste0("n = ", n), levels = n_labels_ordered))

  panel_C <- ggplot(p7_data, aes(x = K, y = mean_time, color = n_label)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ method, scales = "free_y") +
    scale_color_manual(values = n_colors, name = "Matrix size") +
    scale_x_continuous(breaks = seq(2, 10, 2)) +
    labs(x = "K", y = "Time (s)") +
    theme_pub +
    theme(strip.text       = element_text(size = 14, face = "bold"),
          strip.background = element_blank())
} else {
  panel_C <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "No per-K data available")
}

# ==============================================================================
# Combine: (A | B) / C
# ==============================================================================
tag_theme <- theme(plot.tag = element_text(size = 25, face = "bold",
                                           family = "sans"))

panel_A_combo <- panel_A +
  theme(legend.position = "bottom",
        legend.box      = "horizontal",
        legend.margin   = margin(0, 0, 0, 0),
        plot.margin     = margin(5, 10, 5, 5)) +
  tag_theme

panel_B_combo <- panel_B +
  theme(legend.position  = "bottom",
        legend.box       = "horizontal",
        legend.box.just  = "left",
        legend.margin    = margin(0, 0, 0, 0),
        plot.margin      = margin(5, 5, 5, 10)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE),
         shape = guide_legend(nrow = 2)) +
  tag_theme

panel_C_combo <- panel_C +
  theme(legend.position = "bottom",
        legend.box      = "horizontal",
        legend.margin   = margin(0, 0, 0, 0),
        plot.margin     = margin(5, 10, 5, 5)) +
  tag_theme

fig_combined <- (panel_A_combo | panel_B_combo) / panel_C_combo +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(tag_levels = "A")

# ==============================================================================
# Save
# ==============================================================================
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(file.path(output_dir, "fig_scalability.pdf"),
       fig_combined, width = 18, height = 12, dpi = 300)
ggsave(file.path(output_dir, "fig_scalability.png"),
       fig_combined, width = 18, height = 12, dpi = 500)

# Also save individual panels
ggsave(file.path(output_dir, "fig_scalability_A.pdf"),
       panel_A, width = 5.5, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_scalability_B.pdf"),
       panel_B, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(output_dir, "fig_scalability_C.pdf"),
       panel_C, width = 7, height = 3.5, dpi = 300)

cat("Scalability figure saved to", output_dir, "\n")
cat("  Combined: fig_scalability.pdf / .png\n")
cat("  Individual: fig_scalability_A/B/C.pdf\n")
