################################################################################
# BiSer Benchmarking: Visualization Utilities
#
# Shared plotting functions for publication-quality figures.
# Includes: heatmap grid, ridge plots, scatter plots, themes, color palettes.
#
# Dependencies: ComplexHeatmap, circlize, grid, gridExtra, ggplot2, ggridges,
#               ggrepel, cowplot, viridis, dplyr, tidyr
################################################################################

library(ComplexHeatmap)
library(circlize)
library(grid)
library(gridExtra)
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(cowplot)
library(ggridges)
library(tibble)
library(ggrepel)

# ==============================================================================
# Global color palettes and display name mappings
# ==============================================================================

method_colors <- c(
  "BiSer"     = "#E15759",
  "BiSer_TSP" = "#4E79A7",
  "BiSer_SVD" = "#8FAACC",
  "Spectral"  = "#A0D8B0",
  "Heatmap"   = "#3E9E6E",
  "MESBC"     = "#F28E2B",
  "NMF"       = "#FFBE7A"
)

method_order <- c("BiSer", "BiSer_SVD", "BiSer_TSP", "Spectral",
                  "Heatmap", "MESBC", "NMF")

display_name <- c(
  "biser"     = "BiSer",
  "bs"        = "BiSer_SVD",
  "tsp_seri"  = "BiSer_TSP",
  "spec_seri" = "Spectral",
  "Heatmap"   = "Heatmap",
  "MESBC"     = "MESBC",
  "NMF"       = "NMF"
)

methname_old <- c("biser", "bs", "tsp_seri", "spec_seri",
                  "Heatmap", "MESBC", "NMF")

# --- Scalability figure colors ---
step_group_colors <- c(
  "SVD decomposition"       = "#F28E2B",
  "Similarity construction" = "#76B7B2",
  "TSP solving"             = "#EDC948",
  "Other"                   = "#BAB0AC"
)

# --- Publication theme (Nature / Genome Biology style) ---
theme_pub <- theme_classic(base_size = 16, base_family = "sans") +
  theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.line     = element_blank(),
    axis.ticks    = element_line(colour = "black", linewidth = 0.4),
    axis.ticks.length = unit(0.15, "cm"),
    axis.title    = element_text(size = 16, colour = "black"),
    axis.text     = element_text(size = 14, colour = "black"),
    legend.position   = "right",
    legend.background = element_blank(),
    legend.key        = element_blank(),
    legend.key.size   = unit(0.55, "cm"),
    legend.text       = element_text(size = 14),
    legend.title      = element_text(size = 15, face = "bold"),
    plot.margin       = margin(10, 15, 10, 10)
  )

# ==============================================================================
# Heatmap grid (ComplexHeatmap, no annotation bar)
# ==============================================================================
get_multiple_heatmap_grob <- function(mat_list, nrow = 2, ncol = 3,
                                      title, color, legende_title,
                                      title_fontsize = 18,
                                      legend_fontsize = 16) {
  all_vals <- unlist(mat_list)
  min_val  <- min(all_vals, na.rm = TRUE)
  max_val  <- max(all_vals, na.rm = TRUE)
  col_fun  <- colorRamp2(seq(min_val, max_val, length.out = 100), color)

  legend_obj <- Legend(col_fun = col_fun, title = legende_title,
    title_gp  = gpar(fontsize = legend_fontsize, fontface = "bold", fontfamily = "sans"),
    labels_gp = gpar(fontsize = legend_fontsize - 2, fontfamily = "sans"),
    legend_height = unit(4.5, "cm"), grid_width = unit(0.6, "cm"))

  plot_list <- lapply(seq_along(mat_list), function(i) {
    ht <- Heatmap(mat_list[[i]], col = col_fun,
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = FALSE, show_column_names = FALSE,
      column_title = title[i],
      column_title_gp = gpar(fontsize = title_fontsize, fontface = "bold", fontfamily = "sans"),
      heatmap_legend_param = list(title = NULL),
      rect_gp = gpar(col = NA, lwd = 0), name = NULL)
    grid.grabExpr(draw(ht, show_heatmap_legend = FALSE))
  })

  legend_grob <- grid.grabExpr(draw(legend_obj))
  arrangeGrob(
    arrangeGrob(grobs = plot_list, nrow = nrow, ncol = ncol),
    legend_grob, nrow = 1, widths = c(4, 0.5)
  )
}

# ==============================================================================
# Heatmap grid with annotation bar (for real data figures)
# ==============================================================================
get_heatmap_grid_annotated <- function(mat_list, titles, col_label_list,
                                        group_colors, nrow_grid = 2, ncol_grid = 4,
                                        color_palette = colorRampPalette(
                                          c("navy", "white", "firebrick3"))(100),
                                        legend_title = "Expression",
                                        annotation_name = "Group") {
  all_vals <- unlist(mat_list)
  max_abs  <- max(abs(all_vals), na.rm = TRUE)
  col_fun  <- colorRamp2(seq(-max_abs, max_abs, length.out = length(color_palette)),
                         color_palette)

  expr_legend <- Legend(col_fun = col_fun, title = legend_title,
    title_gp = gpar(fontsize = 26, fontface = "bold", fontfamily = "sans"),
    labels_gp = gpar(fontsize = 24, fontfamily = "sans"),
    legend_height = unit(7, "cm"), grid_width = unit(1, "cm"),
    title_gap = unit(0.6, "cm"))

  group_legend <- Legend(labels = names(group_colors),
    legend_gp = gpar(fill = group_colors), title = annotation_name,
    title_gp = gpar(fontsize = 26, fontface = "bold", fontfamily = "sans"),
    labels_gp = gpar(fontsize = 23, fontfamily = "sans"),
    grid_width = unit(0.8, "cm"), grid_height = unit(0.8, "cm"),
    row_gap = unit(0.5, "cm"), title_gap = unit(0.6, "cm"))

  legend_grob <- grid.grabExpr({
    draw(packLegend(expr_legend, group_legend, direction = "vertical",
                    gap = unit(2, "cm")))
  })

  ht_grobs <- lapply(seq_along(mat_list), function(i) {
    anno_list <- list(col_label_list[[i]])
    names(anno_list) <- annotation_name
    top_anno <- HeatmapAnnotation(
      df = data.frame(Group = col_label_list[[i]]),
      col = setNames(list(group_colors), "Group"),
      show_legend = FALSE, show_annotation_name = FALSE,
      annotation_height = unit(0.35, "cm"), simple_anno_size = unit(0.35, "cm"))
    ht <- Heatmap(mat_list[[i]], col = col_fun,
      top_annotation = top_anno,
      cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = FALSE, show_column_names = FALSE,
      column_title = titles[i],
      column_title_gp = gpar(fontsize = 28, fontface = "bold", fontfamily = "sans"),
      heatmap_legend_param = list(title = NULL),
      rect_gp = gpar(col = NA, lwd = 0), name = paste0("ht_", i))
    grid.grabExpr(draw(ht, show_heatmap_legend = FALSE))
  })

  arrangeGrob(
    arrangeGrob(grobs = ht_grobs, nrow = nrow_grid, ncol = ncol_grid),
    legend_grob, nrow = 1, widths = c(4, 1.2)
  )
}

# ==============================================================================
# Ridge plot generator (for simulation metric distribution)
# ==============================================================================
make_ridge_plot <- function(metric_list, core_methods, method_colors,
                             target_metric, x_label, show_y_axis = TRUE,
                             base_size = 17) {


  df_ridge <- bind_rows(lapply(seq_along(metric_list), function(i) {
    as.data.frame(metric_list[[i]]) %>%
      rownames_to_column("Method") %>%
      filter(Method %in% core_methods) %>%
      select(Method, Score = all_of(target_metric)) %>%
      mutate(Simulation_ID = i)
  }))

  df_ridge <- df_ridge %>%
    mutate(Method = case_when(
      tolower(Method) == "biser" ~ "BiSer",
      Method == "tsp_seri"       ~ "BiSer_TSP",
      Method == "bs"             ~ "BiSer_SVD",
      Method == "spec_seri"      ~ "Spectral",
      TRUE                       ~ Method
    ))

  ridge_order <- df_ridge %>%
    group_by(Method) %>%
    summarize(Med = median(Score, na.rm = TRUE), .groups = "drop") %>%
    arrange(Med) %>% pull(Method)

  df_ridge$Method <- factor(df_ridge$Method, levels = ridge_order)
  df_ridge$Score  <- pmin(df_ridge$Score, 1.0)

  p <- ggplot(df_ridge, aes(x = Score, y = Method, fill = Method, color = Method)) +
    geom_vline(xintercept = 1.0, color = "grey50", linetype = "dashed", linewidth = 0.6) +
    stat_density_ridges(
      jittered_points = TRUE, position = position_nudge(y = -0.1),
      point_shape = "|", point_size = 2.5, point_alpha = 0.5,
      alpha = 0.6, scale = 0.9, rel_min_height = 0.005, bandwidth = 0.03
    ) +
    scale_fill_manual(values = method_colors) +
    scale_color_manual(values = method_colors) +
    coord_cartesian(xlim = c(0.3, 1.0)) +
    theme_classic(base_size = base_size) +
    theme(
      legend.position    = "none",
      axis.text.y        = element_text(size = base_size + 1, color = "black"),
      axis.text.x        = element_text(size = base_size, color = "black", face = "bold"),
      axis.line          = element_line(color = "black", linewidth = 0.8),
      axis.ticks         = element_line(color = "black", linewidth = 0.8),
      axis.title.y       = element_blank(),
      axis.title.x       = element_text(margin = margin(t = 12), face = "bold", size = base_size),
      panel.grid.major.x = element_line(color = "grey90", linetype = "dashed"),
      plot.margin        = margin(10, 15, 10, 10)
    ) +
    labs(x = x_label)

  if (!show_y_axis) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }
  return(p)
}

# ==============================================================================
# Scatter plot for ARI vs NMI (or Spearman vs ARI)
# ==============================================================================
make_scatter_plot <- function(metrics_df, x_col = "ARI", y_col = "NMI",
                               x_label = "Adjusted Rand Index (ARI)",
                               y_label = "Normalized Mutual Information (NMI)") {


  metrics_clean <- metrics_df %>%
    mutate(Method    = display_name[Method],
           Highlight = ifelse(Method == "BiSer", "Target", "Others"))

  ggplot(metrics_clean, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray80") +
    geom_point(aes(color = Highlight, size = Highlight, shape = Highlight), alpha = 0.85) +
    scale_color_manual(values = c("Target" = "#DC0000FF", "Others" = "#7E6148FF")) +
    scale_size_manual(values  = c("Target" = 6, "Others" = 4)) +
    scale_shape_manual(values = c("Target" = 16, "Others" = 17)) +
    geom_text_repel(aes(label = Method, color = Highlight),
      size = 8, fontface = "bold", box.padding = 0.6, point.padding = 0.4,
      max.overlaps = 20, segment.color = "gray60", show.legend = FALSE) +
    theme_classic(base_size = 24) +
    theme(legend.position = "none",
      axis.title = element_text(face = "bold", size = 24),
      axis.text  = element_text(color = "black", size = 21),
      panel.grid.major = element_line(color = "gray95"),
      axis.line = element_line(linewidth = 0.8)) +
    labs(x = x_label, y = y_label)
}
