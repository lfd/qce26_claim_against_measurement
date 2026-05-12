#!/usr/bin/env Rscript
# =============================================================================
# Plot script for Experiment Khan et al. Definitive Parameter-Space Study
# Produces 2 publication-quality figures for QCE 2026 paper
#   Figure 3: Parameter-space heatmap  (khan_heatmap.pdf)
#   Figure 4: Shot-count + power       (khan_sensitivity.pdf)
# =============================================================================

source("./reproduction/R/config.R")
library(ggnewscale)

# ── Backend display names ──
# Heatmap uses multi-line labels to save horizontal space
backend_labels_hm <- c(
  "noiseless"   = "Noiseless\n(neg. control)",
  "depol_kyoto" = "Depol.\nKyoto",
  "fake_osaka"  = "Fake\nOsaka",
  "fake_kyoto"  = "Fake\nKyoto"
)

# Sensitivity uses single-line labels for legend
backend_labels_sl <- c(
  "noiseless"   = "Noiseless (neg. ctrl.)",
  "depol_kyoto" = "Depol. Kyoto",
  "fake_osaka"  = "FakeOsaka",
  "fake_kyoto"  = "FakeKyoto"
)
backend_colours_sl <- c(
  "Noiseless (neg. ctrl.)" = LFD$black,
  "Depol. Kyoto"           = LFD$teal,
  "FakeOsaka"              = LFD$orange,
  "FakeKyoto"              = LFD$grey
)

# ═════════════════════════════════════════════════════════════════════════════
#  FIGURE 3: Parameter-Space Heatmap  (dual tile: Cohen's d | Cliff's δ)
# ═════════════════════════════════════════════════════════════════════════════
cat("─── Figure 3: Heatmap ───\n")

# ── Cliff's delta helper (paired-difference vector) ─────────────────────────
cliffs_delta_vec <- function(d) (sum(d > 0) - sum(d < 0)) / length(d)

# ── Load summary (Cohen's d) ─────────────────────────────────────────────────
df <- read.csv(file.path(results_dir, "khan_summary.csv"),
               stringsAsFactors = FALSE)
df$significant <- df$significant == "True"

CONFIG_LEVELS <- c(
  "Default (Richardson)",
  "Fold: from-right", "Fold: global",
  "Extrap: linear", "Extrap: polynomial", "Extrap: exponential",
  "Scales: 1,2,3", "Scales: 1,1.5,2,2.5,3",
  "Transpiler: 3",
  "Shots: 1024", "Shots: 8192"
)

make_config_label <- function(param_axis, param_value) {
  dplyr::case_when(
    param_axis == "default"                                              ~ "Default (Richardson)",
    param_axis == "folding_strategy" & param_value == "from_right"      ~ "Fold: from-right",
    param_axis == "folding_strategy" & param_value == "global"          ~ "Fold: global",
    param_axis == "extrapolation_method" & param_value == "linear"      ~ "Extrap: linear",
    param_axis == "extrapolation_method" & param_value == "polynomial"  ~ "Extrap: polynomial",
    param_axis == "extrapolation_method" & param_value == "exponential" ~ "Extrap: exponential",
    param_axis == "scale_factors" & grepl("1,2,3", param_value)        ~ "Scales: 1,2,3",
    param_axis == "scale_factors" & grepl("1.5",   param_value)        ~ "Scales: 1,1.5,2,2.5,3",
    param_axis == "transpiler_level" & param_value == "3"              ~ "Transpiler: 3",
    param_axis == "n_shots" & param_value == "1024"                    ~ "Shots: 1024",
    param_axis == "n_shots" & param_value == "8192"                    ~ "Shots: 8192",
    TRUE ~ paste0(param_axis, "=", param_value)
  )
}

df <- df %>%
  mutate(
    config_label = factor(make_config_label(param_axis, param_value),
                          levels = CONFIG_LEVELS)
  )

col_order <- expand.grid(
  tc = c("TC1", "TC3", "TC5"),
  bl = backend_labels_hm,
  stringsAsFactors = FALSE
) %>% mutate(col_label = paste0(bl, "\n", tc))

df <- df %>%
  mutate(
    bl        = factor(backend, levels = names(backend_labels_hm),
                       labels = backend_labels_hm),
    col_label = factor(paste0(bl, "\n", circuit), levels = col_order$col_label)
  )

# ── Load detail data → compute Cliff's δ ────────────────────────────────────
detail <- read.csv(file.path(results_dir, "khan_detail.csv"),
                   stringsAsFactors = FALSE) %>%
  mutate(delta = abs(ideal - raw_exp) - abs(ideal - mit_exp))

cliffs_df <- detail %>%
  group_by(backend, circuit, param_axis, param_value) %>%
  summarise(cliffs = cliffs_delta_vec(delta), .groups = "drop")

# ── Join ─────────────────────────────────────────────────────────────────────
df <- df %>%
  left_join(cliffs_df, by = c("backend", "circuit", "param_axis", "param_value")) %>%
  mutate(
    d_label  = ifelse(abs(cohen_d) >= 10,
                      sprintf("%+.0f", cohen_d),
                      sprintf("%+.1f", cohen_d)),
    cd_label = sprintf("%+.2f", cliffs)
  )

# ── Numeric x positions for split tiles ─────────────────────────────────────
x_levels  <- levels(df$col_label)
n_cols    <- length(x_levels)
df$x_num  <- as.numeric(df$col_label)

df_cohen  <- df %>% mutate(x_tile = x_num - 0.25)
df_cliffs <- df %>% mutate(x_tile = x_num + 0.25)
df_box    <- df %>% select(x_num, config_label) %>% distinct()

# Vertical separator positions (between backend groups) adjusted for half-tiles
sep_x <- c(3.5, 6.5, 9.5)

p_heatmap <- ggplot() +
  # ── Left half: Cohen's d ──
  geom_tile(data = df_cohen,
            aes(x = x_tile, y = config_label, fill = cohen_d),
            width = 0.5, colour = NA) +
  geom_text(data = df_cohen,
            aes(x = x_tile, y = config_label, label = d_label),
            size = 1.9, colour = "black") +
  scale_fill_gradient2(
    low      = LFD$orange,
    mid      = "white",
    high     = LFD$teal,
    midpoint = 0,
    name     = "Cohen's $d$ (left)",
    guide    = guide_colourbar(
      barwidth       = unit(3.5, "cm"),
      barheight      = unit(0.28, "cm"),
      title.position = "top", title.hjust = 0.5, order = 1
    )
  ) +
  # ── Switch to second fill scale ──
  new_scale_fill() +
  # ── Right half: Cliff's δ ──
  geom_tile(data = df_cliffs,
            aes(x = x_tile, y = config_label, fill = cliffs),
            width = 0.5, colour = NA) +
  geom_text(data = df_cliffs,
            aes(x = x_tile, y = config_label, label = cd_label),
            size = 1.9, colour = "black") +
  scale_fill_gradient2(
    low      = LFD$orange,
    mid      = "white",
    high     = LFD$teal,
    midpoint = 0,
    name     = "Cliff's $\\delta$ (right)",
    guide    = guide_colourbar(
      barwidth       = unit(3.5, "cm"),
      barheight      = unit(0.28, "cm"),
      title.position = "top", title.hjust = 0.5, order = 2
    )
  ) +
  # ── Box outlines around each paired (Cohen's d | Cliff's δ) cell ──
  geom_tile(data = df_box,
            aes(x = x_num, y = config_label),
            fill = NA, colour = "grey55", linewidth = 0.35, width = 1.0) +
  # ── Vertical separators between backend groups ──
  geom_vline(xintercept = sep_x, colour = "grey30", linewidth = 0.5) +
  # ── Axes ──
  scale_x_continuous(
    breaks = seq_len(n_cols),
    labels = x_levels,
    expand = c(0.02, 0)
  ) +
  scale_y_discrete(limits = rev) +
  labs(x = NULL, y = NULL) +
  theme_paper() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = BASE.SIZE,
                                    lineheight = 0.85),
    legend.position  = "top",
    legend.direction = "horizontal",
    panel.grid       = element_blank()
  )

save_plot(p_heatmap, "khan_heatmap", width = TEXTWIDTH, height = 0.45 * TEXTWIDTH)
cat("  ✓ khan_heatmap\n")


# ═════════════════════════════════════════════════════════════════════════════
#  FIGURE 4: Combined sensitivity (shot-count + power)
# ═════════════════════════════════════════════════════════════════════════════
cat("─── Figure 4: Sensitivity ───\n")

# ── Panel (a): Shot-count vs Cohen's d ──
ds <- read.csv(file.path(results_dir, "khan_shots.csv"),
               stringsAsFactors = FALSE) %>%
  mutate(backend_label = factor(backend,
                                levels = names(backend_labels_sl),
                                labels = backend_labels_sl))

# Two-row facet: noise backends (large d) vs control/broken (small d)
ds <- ds %>%
  mutate(group = case_when(
    backend %in% c("depol_kyoto", "fake_osaka") ~ "Noise backends (genuine effect)",
    TRUE ~ "Neg. control / broken calibration"
  ),
  group = factor(group, levels = c("Noise backends (genuine effect)",
                                   "Neg. control / broken calibration")))

p_shots <- ggplot(ds, aes(x = n_shots, y = cohen_d,
                           colour = backend_label,
                           linetype = circuit)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.0) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50",
             linewidth = 0.3) +
  facet_wrap(~ group, ncol = 1, scales = "free_y") +
  scale_x_log10(breaks = c(128, 256, 512, 1024, 2048, 4096, 8192),
                labels = c("128", "256", "512", "1K", "2K", "4K", "8K")) +
  scale_colour_manual(values = backend_colours_sl, name = "Backend") +
  scale_linetype_manual(values = c("TC1" = "solid",
                                    "TC3" = "dashed",
                                    "TC5" = "dotted"),
                        name = "Trotter depth") +
  labs(x = "$n_{\\text{shots}}$", y = "Cohen's d",
       title = "(a) Effect size vs. shot count") +
  theme_paper() +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.3, "cm"),
    plot.title       = element_text(size = BASE.SIZE, face = "bold"),
    strip.text       = element_text(size = SMALL.SIZE, face = "italic")
  ) +
  guides(colour   = guide_legend(nrow = 1, order = 1),
         linetype = guide_legend(nrow = 1, order = 2))


# ── Panel (b): Power curves ──
dp <- read.csv(file.path(results_dir, "khan_power.csv"),
               stringsAsFactors = FALSE) %>%
  mutate(backend_label = factor(backend,
                                levels = names(backend_labels_sl),
                                labels = backend_labels_sl))

p_power <- ggplot(dp, aes(x = n_reps, y = power,
                            colour = backend_label,
                            linetype = circuit)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.0) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50",
             linewidth = 0.3) +
  annotate("text", x = 5, y = 0.83, label = "80\\% power", size = 2,
           colour = "grey50", hjust = 0) +
  scale_x_log10(breaks = c(5, 10, 20, 30, 50, 100, 200)) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
  scale_colour_manual(values = backend_colours_sl, name = "Backend") +
  scale_linetype_manual(values = c("TC1" = "solid",
                                    "TC3" = "dashed",
                                    "TC5" = "dotted"),
                        name = "Trotter depth") +
  labs(x = "$n_{\\text{reps}}$", y = "Statistical power",
       title = "(b) Power vs. repetition count") +
  theme_paper() +
  theme(
    legend.position  = "bottom",
    legend.key.size  = unit(0.3, "cm"),
    plot.title       = element_text(size = BASE.SIZE, face = "bold")
  ) +
  guides(colour   = guide_legend(nrow = 1, order = 1),
         linetype = guide_legend(nrow = 1, order = 2))

# Combine with shared legend
p_combined <- p_shots + plot_spacer() + p_power +
  plot_layout(ncol = 3, widths = c(1, 0.05, 1), guides = "collect") &
  theme(legend.position = "bottom")

save_plot(p_combined, "khan_sensitivity", width = TEXTWIDTH, height = 0.4 * TEXTWIDTH)
cat("  ✓ khan_sensitivity\n")

cat("\n═══ Done! ═══\n")
