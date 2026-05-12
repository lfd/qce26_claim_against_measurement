#!/usr/bin/env Rscript
# =============================================================================
# Plot Generator: Desdentado Reanalysis — Authors' Original Data
# =============================================================================
# Reads CSV output from generate_desdentado_original_boxplot.py and produces
# a publication-quality figure.
#
# Output:
#   plots/desdentado_sweep.pdf       — Boxplot of success rate by shot count
#   img-tikz/desdentado_sweep.tex    — TikZ export for paper inclusion
# =============================================================================

source("./reproduction/R/config.R")

# Reference line
P_RANDOM <- 2/32  # 0.0625 success probability


# =============================================================================
# FIGURE 5: Success Rate by Shot Count (Boxplot + Jitter)
# =============================================================================
box_file <- file.path(results_dir, "desdentado_original_boxplot.csv")
if (file.exists(box_file)) {
  cat("Plotting shot-count boxplot...\n")
  df <- read_csv(box_file, show_col_types = FALSE) %>%
    mutate(
      shot_label = factor(
        paste0(shot_count, "\n(", label, ")"),
        levels = c("1024\n(Default)", "3502\n(Mid-Low)",
                    "5980\n(Estimated)", "8458\n(Mid-High)",
                    "10936\n(High)")
      )
    )

  p1 <- ggplot(df, aes(x = shot_label, y = success_rate)) +
    geom_boxplot(
      width = 0.5, outlier.shape = NA,
      fill = NA, colour = LFD$black, linewidth = 0.4
    ) +
    geom_jitter(width = 0.12, size = POINT.SIZE, alpha = 0.7, colour = LFD$black) +
    geom_hline(
      yintercept = P_RANDOM, linetype = "dashed",
      colour = LFD$grey, linewidth = LINE.SIZE
    ) +
    # annotate("text", x = 5.45, y = P_RANDOM + 0.001,
    #          label = "Random (2/32)",
    #          colour = LFD$grey, hjust = 1, vjust = 0, size = 2.5) +
    labs(
      x = "Shot count (configuration)",
      y = "Success rate"
    ) +
    coord_cartesian(ylim = c(min(df$success_rate) - 0.003,
                              max(df$success_rate) + 0.005)) +
    theme_paper()

  save_plot(p1, "desdentado_sweep", width = COLWIDTH, height = 0.4 * COLWIDTH)
  cat("  -> desdentado_sweep\n")
} else {
  cat("  [skip] boxplot CSV not found\n")
}

cat("\nDone.\n")
