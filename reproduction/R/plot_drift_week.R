#!/usr/bin/env Rscript
# =============================================================================
# Week Drift Figure — IQM QExa 7-Day Longitudinal Study
# =============================================================================
# Single-panel timeseries of E(lambda=1) over 163.5 scheduled hours,
# with night shading, red gap band (Easter 2026 QPU outage), and a broken
# line that does NOT connect across the 43-hour gap (h69-h112).
#
# Input:  data/qexa_drift/raw_data_week.csv
# Output: plots/drift_week.pdf + plots/drift_week.tex
# =============================================================================

.script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m))) else "R"
})
source(file.path(.script_dir, "config.R"))

# =============================================================================
# Load & Process Data
# =============================================================================
data_dir <- file.path(.script_dir, "..", "data", "qexa_drift")
raw <- read.csv(file.path(data_dir, "raw_data_week.csv"), stringsAsFactors = FALSE)

cat(sprintf("Raw rows: %d, TPs: %d\n", nrow(raw), length(unique(raw$timepoint_idx))))

# Per-TP mean and SE of E(lambda=1), with scheduled time axis
tp <- raw %>%
  filter(scale_factor == 1.0) %>%
  group_by(timepoint_idx) %>%
  summarise(
    hours_scheduled = first(timepoint_idx) * 0.5,
    mean_e_l1       = mean(exp_val),
    se_e_l1         = sd(exp_val) / sqrt(n()),
    .groups         = "drop"
  ) %>%
  arrange(timepoint_idx)

cat(sprintf("Week drift: %d TPs, %.1f scheduled hours\n",
            nrow(tp), max(tp$hours_scheduled)))

# =============================================================================
# Gap & Night Definitions
# =============================================================================

# QPU outage: Sat 11 April 14:15 — Mon 13 April 15:29 CEST
GAP_START <- 69.0   # h69.0  (last successful TP 138)
GAP_END   <- 112.0  # h112.0 (first recovered TP 224)

# Experiment start: 2026-04-08 14:54 UTC = 16:54 CEST
# Night = 21:00-06:00 CEST => offset from t0:
#   first 21:00 at t0 + (21:00 - 16:54) = t0 + 4.10h
NIGHT_OFFSET_START <- 4.10   # hours from t0 to first 21:00 local
NIGHT_DURATION     <- 9.00   # 21:00 to 06:00 = 9 h

# Build all night windows within and beyond the data span
MAX_NIGHTS <- 8
night_starts <- NIGHT_OFFSET_START + (0:(MAX_NIGHTS - 1)) * 24
night_ends   <- night_starts + NIGHT_DURATION
night_df <- data.frame(
  night   = seq_along(night_starts),
  xmin    = night_starts,
  xmax    = night_ends,
  ymin    = -Inf,
  ymax    =  Inf
)
# Keep only nights that overlap with the data span
MAX_H <- max(tp$hours_scheduled) + 0.5
night_df <- night_df[night_df$xmin <= MAX_H, ]

# =============================================================================
# Split data at gap (so geom_line does NOT span the gap)
# =============================================================================
tp$segment <- ifelse(tp$hours_scheduled <= GAP_START, "before", "after")

# =============================================================================
# Trend-line intervals — EDIT this list only.
# Each entry c(lo, hi) defines one interval for a linear trend line fit.
# Transition jumps (e.g. h19-h24 drop) are intentionally omitted.
# =============================================================================
TREND_INTERVALS <- list(
  c(  0.0,  19.0),         # stable high pre-drop
  c( 24.0,  43.5),         # low stable plateau
  c( 43.5,  GAP_START),    # elevated pre-gap rise
  c( GAP_END, 149.5),      # stable low post-outage
  c(149.5,  MAX_H + 0.1)   # elevated post-gap rise
)

# Assign each TP to its trend interval
tp$tl_id <- NA_integer_
for (.j in seq_along(TREND_INTERVALS)) {
  .iv <- TREND_INTERVALS[[.j]]
  .m  <- tp$hours_scheduled >= .iv[1] & tp$hours_scheduled < .iv[2]
  tp$tl_id[.m] <- .j
}
cat(sprintf("\nTrend intervals assigned: %d TPs\n", sum(!is.na(tp$tl_id))))

# =============================================================================
# Plot
# =============================================================================
ts_ylim <- range(c(tp$mean_e_l1 - 1.96 * tp$se_e_l1,
                   tp$mean_e_l1 + 1.96 * tp$se_e_l1),
                 na.rm = TRUE) + c(-0.005, 0.015)

# Night label y position: just inside top of panel
label_y <- ts_ylim[2] - 0.005

# Classify nights: those overlapping with gap shown in red text
night_df$in_gap <- night_df$xmin >= GAP_START & night_df$xmax <= GAP_END

p <- ggplot(tp, aes(x = hours_scheduled, y = mean_e_l1))

# Layer 1: night shading (grey)
p <- p + geom_rect(
  data = night_df,
  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  inherit.aes = FALSE, fill = "grey90", alpha = 0.5
)

# Layer 2: gap shading (red)
p <- p + annotate(
  "rect",
  xmin = GAP_START, xmax = GAP_END, ymin = -Inf, ymax = Inf,
  fill = LFD$red, alpha = 0.10
)

# Layer 3: night labels
for (i in seq_len(nrow(night_df))) {
  nx  <- night_df$xmin[i] + NIGHT_DURATION / 2
  lbl <- paste0("Night ", night_df$night[i])
  col <- if (night_df$in_gap[i]) LFD$red else "grey50"
  p <- p + annotate(
    "text", x = nx, y = label_y,
    label = lbl, size = 2, colour = col, fontface = "italic"
  )
}

# Layer 4: gap label
gap_center <- (GAP_START + GAP_END) / 2
p <- p + annotate(
  "text", x = gap_center, y = ts_ylim[1] + 0.01,
  label = "QPU outage", size = 2.2, colour = LFD$red, fontface = "bold"
)

# Layer 5: 95% CI ribbon (two segments)
p <- p + geom_ribbon(
  data = tp[tp$segment == "before", ],
  aes(ymin = mean_e_l1 - 1.96 * se_e_l1,
      ymax = mean_e_l1 + 1.96 * se_e_l1),
  fill = "grey60", alpha = 0.20
)
p <- p + geom_ribbon(
  data = tp[tp$segment == "after", ],
  aes(ymin = mean_e_l1 - 1.96 * se_e_l1,
      ymax = mean_e_l1 + 1.96 * se_e_l1),
  fill = "grey60", alpha = 0.20
)

# Layer 6: timeseries line (two segments — no line drawn across gap)
p <- p + geom_line(
  data = tp[tp$segment == "before", ],
  colour = "black", linewidth = 0.5
)
p <- p + geom_line(
  data = tp[tp$segment == "after", ],
  colour = "black", linewidth = 0.5
)

# Layer 7: points
p <- p + geom_point(
  data = tp[tp$segment == "before", ],
  colour = "black", size = POINT.SIZE
)
p <- p + geom_point(
  data = tp[tp$segment == "after", ],
  colour = "black", size = POINT.SIZE
)

# Layer 8: linear trend lines — one per interval, no continuity between them
p <- p + geom_smooth(
  data    = tp[!is.na(tp$tl_id), ],
  aes(x = hours_scheduled, y = mean_e_l1, group = tl_id),
  method  = "lm", formula = y ~ x,
  colour  = "grey40", linewidth = 1, linetype = "dashed",
  se      = FALSE, inherit.aes = FALSE
)

# Axis formatting
x_breaks <- seq(0, ceiling(MAX_H / 24) * 24, by = 24)
x_breaks <- x_breaks[x_breaks <= ceiling(MAX_H)]

p <- p +
  scale_x_continuous(
    breaks = x_breaks,
    labels = paste0(x_breaks, "h"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  coord_cartesian(ylim = ts_ylim) +
  labs(
    x = "Scheduled hours since start",
    y = "$\\bar{E}(\\lambda_1)$"
  ) +
  theme_paper() +
  theme(
    legend.position = "none",
    panel.grid.minor.x = element_line(linewidth = 0.15, colour = "grey95")
  )

# =============================================================================
# Save
# =============================================================================
save_plot(p, "drift_week", width = TEXTWIDTH, height = 0.25 * TEXTWIDTH)
