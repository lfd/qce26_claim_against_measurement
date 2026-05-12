#!/usr/bin/env Rscript
# =============================================================================
# IBM Brussels Drift Figure — 12-Hour Longitudinal Session
# =============================================================================
# Two-panel figure for the IBM ibm_brussels hardware drift evidence:
#   (a) E(lambda=1) timeseries with 95% CI ribbon + linear trend
#   (b) ZNE verdict stability (Cohen's d) per timepoint
#
# Input:  paper_qce/ibm_drift_results/raw_data_ibm_drift.csv
# Output: build/plots/drift_ibm_brussels.tex  (+ PDF preview)
# =============================================================================

.script_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) {
  args <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", args, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m))) else "R"
})
source(file.path(.script_dir, "config.R"))

IDEAL     <- 0.845656    # ideal ⟨ZZZZ⟩ for TC3 with calibrated angles
SCALE_FACTORS <- c(1.0, 3.0, 5.0)
COLOUR    <- LFD$black   # IBM session colour — black to match week drift style — black to match week drift style

# Richardson extrapolation coefficients for lambda = {1, 3, 5}
# Lagrange interpolant through (1,E1),(3,E3),(5,E5) evaluated at lambda=0
RICH_C1 <-  15 / 8
RICH_C3 <- -10 / 8
RICH_C5 <-   3 / 8

# =============================================================================
# Load Data
# =============================================================================
data_file <- file.path("reproduction", "data", "ibm_drift_results", "raw_data_ibm_drift.csv")
if (!file.exists(data_file)) {
  stop("Data not found: ", data_file,
       "\nRun: ./scripts/deploy_ibm_drift.sh --download")
}

raw <- read.csv(data_file, stringsAsFactors = FALSE)

# Remove the broken TP0 (wrong angles, all vals ≈ 0)
raw$timestamp <- as.POSIXct(raw$timestamp, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")

# Keep only the valid run: starts at 2026-04-23T16:37 onward
t_valid_start <- as.POSIXct("2026-04-23T16:37:00", format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
raw <- raw[raw$timestamp >= t_valid_start, ]

# Rebase time to hours since start
t0 <- min(raw$timestamp)
raw$hours <- as.numeric(difftime(raw$timestamp, t0, units = "hours"))

cat(sprintf("Data loaded: %d rows, %d TPs, %.1f h\n",
            nrow(raw),
            length(unique(raw$timepoint_idx)),
            max(raw$hours)))

# =============================================================================
# Panel (a): E(lambda=1) timeseries
# =============================================================================
tp_means_l1 <- raw %>%
  filter(scale_factor == 1.0) %>%
  group_by(timepoint_idx) %>%
  summarise(
    mean_e = mean(exp_val),
    sd_e   = sd(exp_val),
    se_e   = sd(exp_val) / sqrt(n()),
    hours  = first(hours),
    .groups = "drop"
  )

cat(sprintf("\nE(lambda=1) range: [%.4f, %.4f], mean = %.4f\n",
            min(tp_means_l1$mean_e), max(tp_means_l1$mean_e),
            mean(tp_means_l1$mean_e)))

ts_ylim <- range(c(tp_means_l1$mean_e - 1.96 * tp_means_l1$se_e,
                    tp_means_l1$mean_e + 1.96 * tp_means_l1$se_e)) +
           c(-0.01, 0.02)

pa <- ggplot(tp_means_l1, aes(x = hours, y = mean_e)) +
  geom_hline(yintercept = IDEAL, linetype = "dotted",
             colour = "grey50", linewidth = 0.35) +
  annotate("text", x = max(tp_means_l1$hours) * 0.98, y = IDEAL + 0.012,
           label = "ideal", size = SMALL.SIZE / ggplot2::.pt,
           colour = "grey50", hjust = 1, fontface = "italic") +
  geom_ribbon(aes(ymin = mean_e - 1.96 * se_e,
                  ymax = mean_e + 1.96 * se_e),
              fill = COLOUR, alpha = 0.15) +
  geom_line(colour = COLOUR, linewidth = 0.6) +
  geom_point(colour = COLOUR, size = 1.0) +
  geom_smooth(method = "lm", formula = y ~ x,
              colour = "grey30", linewidth = 0.45,
              linetype = "dashed", se = FALSE) +
  scale_x_continuous(breaks = seq(0, 12, 3),
                     labels = paste0(seq(0, 12, 3), "h")) +
  coord_cartesian(ylim = ts_ylim) +
  labs(
       x = "Hours since session start",
       y = "$\\bar{E}(\\lambda_1)$") +
  theme_paper() +
  theme(legend.position = "none",
        plot.title = element_text(size = BASE.SIZE, face = "bold"),
        panel.grid.minor.x = element_line(linewidth = 0.15, colour = "grey95"))

# =============================================================================
# Panel (b): ZNE verdict (Cohen's d)
# =============================================================================
zne_per_rep <- raw %>%
  pivot_wider(
    id_cols     = c(timepoint_idx, hours, rep),
    names_from  = scale_factor,
    values_from = exp_val,
    names_prefix = "lam",
    names_repair = "minimal"
  ) %>%
  rename(lam1 = `lam1`, lam3 = `lam3`, lam5 = `lam5`) %>%
  mutate(
    zne       = RICH_C1 * lam1 + RICH_C3 * lam3 + RICH_C5 * lam5,
    raw_error = abs(lam1 - IDEAL),
    mit_error = abs(zne  - IDEAL),
    delta     = raw_error - mit_error   # positive = ZNE helped
  )

verdict_data <- zne_per_rep %>%
  group_by(timepoint_idx, hours) %>%
  summarise(
    mean_zne   = mean(zne),
    mean_raw   = mean(lam1),
    mean_delta = mean(delta),
    sd_delta   = sd(delta),
    cohens_d   = mean(delta) / sd(delta),
    t_stat     = mean(delta) / (sd(delta) / sqrt(n())),
    p_val      = pt(-abs(mean(delta) / (sd(delta) / sqrt(n()))), df = n() - 1) * 2,
    n          = n(),
    .groups    = "drop"
  ) %>%
  mutate(
    verdict = case_when(
      p_val < 0.05 & cohens_d > 0 ~ "Sig. better",
      p_val < 0.05 & cohens_d < 0 ~ "Sig. worse",
      TRUE ~ "Not significant"
    )
  )

cat(sprintf("\nZNE verdict summary:\n"))
print(table(verdict_data$verdict))
cat(sprintf("\nCohen's d: [%.2f, %.2f], mean = %.2f\n",
            min(verdict_data$cohens_d), max(verdict_data$cohens_d),
            mean(verdict_data$cohens_d)))
cat(sprintf("Mean ZNE: %.4f  (%.1f%% ideal, err = %.4f)\n",
            mean(verdict_data$mean_zne),
            mean(verdict_data$mean_zne) / IDEAL * 100,
            abs(mean(verdict_data$mean_zne) - IDEAL)))
cat(sprintf("Mean raw: %.4f  (%.1f%% ideal, err = %.4f)\n",
            mean(verdict_data$mean_raw),
            mean(verdict_data$mean_raw) / IDEAL * 100,
            abs(mean(verdict_data$mean_raw) - IDEAL)))

cd_ylim <- range(verdict_data$cohens_d) + c(-0.3, 0.5)

pc <- ggplot(verdict_data, aes(x = hours, y = cohens_d)) +
  geom_col(fill = COLOUR, colour = COLOUR, width = 0.35, alpha = 0.7) +
  geom_hline(yintercept = 0.8, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 2.0, linetype = "dashed",
             colour = "grey60", linewidth = 0.3) +
  annotate("text", x = max(verdict_data$hours) * 0.98, y = 0.8 + 0.15,
           label = "large", size = SMALL.SIZE / ggplot2::.pt,
           colour = "grey50", hjust = 1, fontface = "italic") +
  annotate("text", x = max(verdict_data$hours) * 0.98, y = 2.0 + 0.15,
           label = "huge", size = SMALL.SIZE / ggplot2::.pt,
           colour = "grey50", hjust = 1, fontface = "italic") +
  scale_x_continuous(breaks = seq(0, 12, 3),
                     labels = paste0(seq(0, 12, 3), "h")) +
  coord_cartesian(ylim = cd_ylim) +
  labs(title = "(b) ZNE verdict stability — Cohen's $d$",
       x = "Hours since session start",
       y = "Cohen's $d$") +
  theme_paper() +
  theme(legend.position = "none",
        plot.title = element_text(size = BASE.SIZE, face = "bold"),
        panel.grid.minor.x = element_line(linewidth = 0.15, colour = "grey95"))

# =============================================================================
# Save: single panel
# =============================================================================
save_plot(pa, "drift_ibm_brussels",
          width = COLWIDTH, height = 0.467 * COLWIDTH)

cat("\nSaved: build/plots/drift_ibm_brussels.tex\n")
