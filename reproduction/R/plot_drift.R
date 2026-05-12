#!/usr/bin/env Rscript
# =============================================================================
# Combined Drift Figure — All Three Sessions (Layer 2)
# =============================================================================
# Two-panel figure combining all drift evidence:
#   (a) E(lambda=1) timeseries — all 3 sessions, time rebased to 0
#   (b) ZNE verdict stability (Cohen's d) — all 3 sessions, 0–48h
#
# Input:  data/qexa_drift/raw_data_first_run.csv  (Day 1, 12h)
#         data/qexa_drift/raw_data_day2_full.csv   (Day 2, 12h)
#         data/qexa_drift/raw_data_weekend.csv     (Weekend, 48h)
# Output: plots/drift_combined_all.pdf
# =============================================================================

source("./reproduction/R/config.R")

# All sessions in black (monochrome)
session_colours <- c("Day 1 (12 h)" = LFD$orange,
                     "Day 2 (12 h)" = LFD$grey,
                     "Weekend (48 h)" = LFD$teal)

# =============================================================================
# Load Data
# =============================================================================
data_dir <- file.path(".", "reproduction", "data", "qexa_drift")

d1 <- read.csv(file.path(data_dir, "raw_data_first_run.csv"), stringsAsFactors = FALSE)
d2 <- read.csv(file.path(data_dir, "raw_data_day2_full.csv"), stringsAsFactors = FALSE)
dw <- read.csv(file.path(data_dir, "raw_data_weekend.csv"), stringsAsFactors = FALSE)

d1$session <- "Day 1 (12 h)"
d2$session <- "Day 2 (12 h)"
dw$session <- "Weekend (48 h)"

# Parse timestamps and rebase to hours since each session's start
parse_and_rebase <- function(df) {
  df$timestamp <- as.POSIXct(df$timestamp, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")
  t0 <- min(df$timestamp)
  df$hours <- as.numeric(difftime(df$timestamp, t0, units = "hours"))
  df
}

d1 <- parse_and_rebase(d1)
d2 <- parse_and_rebase(d2)
dw <- parse_and_rebase(dw)

ideal_val <- d1$ideal[1]  # 0.981197

cat(sprintf("Day 1:    %d rows, %d TPs, %.1f h\n", nrow(d1), max(d1$timepoint_idx)+1, max(d1$hours)))
cat(sprintf("Day 2:    %d rows, %d TPs, %.1f h\n", nrow(d2), max(d2$timepoint_idx)+1, max(d2$hours)))
cat(sprintf("Weekend:  %d rows, %d TPs, %.1f h\n", nrow(dw), max(dw$timepoint_idx)+1, max(dw$hours)))

# =============================================================================
# Data Processing
# =============================================================================

compute_tp_means_l1 <- function(df) {
  df %>%
    filter(scale_factor == 1.0) %>%
    group_by(timepoint_idx, session) %>%
    summarise(
      mean_e  = mean(exp_val),
      sd_e    = sd(exp_val),
      se_e    = sd(exp_val) / sqrt(n()),
      hours   = first(hours),
      .groups = "drop"
    )
}

tp1 <- compute_tp_means_l1(d1)
tp2 <- compute_tp_means_l1(d2)
tpw <- compute_tp_means_l1(dw)
tp_all <- bind_rows(tp1, tp2, tpw)

# Richardson ZNE coefficients for lambda = {1, 3, 5}
richardson_coeffs <- c(`1` = 15/8, `3` = -10/8, `5` = 3/8)

compute_verdict_data <- function(df) {
  zne_per_rep <- df %>%
    pivot_wider(
      id_cols = c(timepoint_idx, hours, rep, session),
      names_from = scale_factor,
      values_from = exp_val,
      names_prefix = "lam"
    ) %>%
    mutate(
      zne = richardson_coeffs["1"] * lam1 +
            richardson_coeffs["3"] * lam3 +
            richardson_coeffs["5"] * lam5,
      raw_error = abs(lam1 - ideal_val),
      mit_error = abs(zne - ideal_val),
      delta = raw_error - mit_error
    )

  zne_per_rep %>%
    group_by(timepoint_idx, hours, session) %>%
    summarise(
      mean_raw    = mean(lam1),
      mean_zne    = mean(zne),
      mean_delta  = mean(delta),
      sd_delta    = sd(delta),
      cohens_d    = mean(delta) / sd(delta),
      t_stat      = mean(delta) / (sd(delta) / sqrt(n())),
      p_val       = pt(-abs(t_stat), df = n() - 1) * 2,
      n           = n(),
      .groups     = "drop"
    ) %>%
    mutate(
      verdict = case_when(
        p_val < 0.05 & cohens_d > 0 ~ "Sig. better",
        p_val < 0.05 & cohens_d < 0 ~ "Sig. worse",
        TRUE ~ "Not significant"
      )
    )
}

v1 <- compute_verdict_data(d1)
v2 <- compute_verdict_data(d2)
vw <- compute_verdict_data(dw)
v_all <- bind_rows(v1, v2, vw)

cat("\nVerdict summary across sessions:\n")
print(table(v_all$session, v_all$verdict))

cat(sprintf("\nCohen's d ranges:\n"))
for (s in names(session_colours)) {
  vs <- v_all %>% filter(session == s)
  cat(sprintf("  %s: d = [%.2f, %.2f], mean = %.2f\n",
              s, min(vs$cohens_d), max(vs$cohens_d), mean(vs$cohens_d)))
}

# Shared y-limits across panels for consistent comparison
ts_ylim <- range(c(tp_all$mean_e - 1.96 * tp_all$se_e,
                    tp_all$mean_e + 1.96 * tp_all$se_e)) + c(-0.002, 0.012)
cd_ylim <- range(v_all$cohens_d) + c(-0.2, 0.3)

# Night shading (weekend panels only)
night_rects <- data.frame(
  xmin = c(11.2, 35.2),
  xmax = c(19.2, 43.2),
  ymin = -Inf, ymax = Inf
)


# =============================================================================
# Row 1: E(lambda=1) timeseries — one panel per session
# =============================================================================

make_ts_panel <- function(tp, colour, title, xbreaks, show_ylab,
                          night = NULL, trend_segments = NULL) {
  p <- ggplot(tp, aes(x = hours, y = mean_e))
  if (!is.null(night)) {
    p <- p +
      geom_rect(data = night,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "grey90", alpha = 0.5) +
      annotate("text", x = 15.2, y = ts_ylim[2] - 0.002,
               label = "Night 1", size = 2, colour = "grey50", fontface = "italic") +
      annotate("text", x = 39.2, y = ts_ylim[2] - 0.002,
               label = "Night 2", size = 2, colour = "grey50", fontface = "italic")
  }
  p <- p +
    geom_ribbon(aes(ymin = mean_e - 1.96 * se_e, ymax = mean_e + 1.96 * se_e),
                fill = colour, alpha = 0.15) +
    geom_line(colour = colour, linewidth = 0.6) +
    geom_point(colour = colour, size = 1.0)

  # Linear trend(s)
  if (!is.null(trend_segments)) {
    for (seg in trend_segments) {
      seg_data <- tp %>% filter(hours >= seg[1], hours <= seg[2])
      if (nrow(seg_data) >= 2) {
        p <- p + geom_smooth(data = seg_data, method = "lm", formula = y ~ x,
                             colour = "grey30", linewidth = 0.5,
                             linetype = "dashed", se = FALSE)
      }
    }
  } else {
    p <- p + geom_smooth(method = "lm", formula = y ~ x,
                         colour = "grey30", linewidth = 0.5,
                         linetype = "dashed", se = FALSE)
  }

  p <- p +
    scale_x_continuous(breaks = xbreaks, labels = paste0(xbreaks, "h")) +
    coord_cartesian(ylim = ts_ylim) +
    labs(title = title, x = NULL,
         y = if (show_ylab) "$\\bar{E}(\\lambda_1)$" else NULL) +
    theme_paper() +
    theme(legend.position = "none",
          plot.title = element_text(size = BASE.SIZE, face = "bold"),
          panel.grid.minor.x = element_line(linewidth = 0.15, colour = "grey95"))
  if (!show_ylab) p <- p + theme(axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank())
  p
}

pa1 <- make_ts_panel(tp1, LFD$orange, "(a) Session 1 (12 h)",
                     xbreaks = seq(0, 12, 3), show_ylab = TRUE)
pa2 <- make_ts_panel(tp2, LFD$grey, "(b) Session 2 (12 h)",
                     xbreaks = seq(0, 12, 3), show_ylab = FALSE)
paw <- make_ts_panel(tpw, LFD$teal, "(c) Weekend (48 h)",
                     xbreaks = seq(0, 48, 12), show_ylab = FALSE,
                     night = night_rects,
                     trend_segments = list(c(0, 40), c(43, 48)))


# =============================================================================
# Row 2: ZNE verdict (Cohen's d) — one panel per session
# =============================================================================

make_verdict_panel <- function(vd, colour, title, xbreaks, show_ylab, night = NULL) {
  p <- ggplot(vd, aes(x = hours, y = cohens_d))
  if (!is.null(night)) {
    p <- p +
      geom_rect(data = night,
                aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                inherit.aes = FALSE, fill = "grey90", alpha = 0.5)
  }
  p <- p +
    geom_col(fill = colour, colour = colour, width = 0.35, alpha = 0.7) +
    # geom_point(colour = colour, size = 0.8) +
    geom_hline(yintercept = c(0.8, 2.0), linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    scale_x_continuous(breaks = xbreaks, labels = paste0(xbreaks, "h")) +
    coord_cartesian(ylim = cd_ylim) +
    # labs(title = title, x = "Hours since start",
    #      y = if (show_ylab) "Cohen's d" else NULL) +
      labs(x = NULL,
         y = if (show_ylab) "Cohen's d" else NULL) +
    theme_paper() +
    theme(legend.position = "none",
          plot.title = element_text(size = BASE.SIZE, face = "bold"),
          panel.grid.minor.x = element_line(linewidth = 0.15, colour = "grey95"))
  if (!show_ylab) p <- p + theme(axis.title.y = element_blank(),
                                  axis.text.y = element_blank(),
                                  axis.ticks.y = element_blank())
  p
}

pc1 <- make_verdict_panel(v1, LFD$orange, "(d) Session 1",
                          xbreaks = seq(0, 12, 3), show_ylab = TRUE)
pc2 <- make_verdict_panel(v2, LFD$grey, "(e) Session 2",
                          xbreaks = seq(0, 12, 3), show_ylab = FALSE)
pcw <- make_verdict_panel(vw, LFD$teal, "(f) Weekend",
                          xbreaks = seq(0, 48, 12), show_ylab = FALSE,
                          night = night_rects)

# # Add "large"/"huge" annotations to the weekend verdict panel
# pcw <- pcw +
#   annotate("text", x = 47, y = 0.8, label = "large",
#            colour = "grey50", size = 2, hjust = 1, vjust = -0.3) +
#   annotate("text", x = 47, y = 2.0, label = "huge",
#            colour = "grey50", size = 2, hjust = 1, vjust = -0.3)


# =============================================================================
# Combine: 3 x 2 grid  (1/4 | 1/4 | 1/2)
# =============================================================================

design <- "
AABBCCCC
DDEEFFFF
"

p_grid <- (pa1 + pa2 + paw + pc1 + pc2 + pcw +
    plot_layout(design = design)) +
  plot_annotation(
    caption = "Hours since start",
    theme = theme(
      plot.caption = element_text(hjust = 0.5, size = BASE.SIZE,
                                  margin = margin(t = 2))
    )
  )

save_plot(p_grid, "drift_combined_all", width = TEXTWIDTH, height = 0.42 * TEXTWIDTH)
