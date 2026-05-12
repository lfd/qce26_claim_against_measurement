#!/usr/bin/env Rscript
# =============================================================================
# Plot 1: Systematic Review — Criterion Compliance (Stacked Bar)
# =============================================================================
# Input:  data/review_criteria.csv
# Output: plots/review_compliance.pdf
# =============================================================================

source("./reproduction/R/config.R")

# ── Criterion labels ──
criterion_names <- c(
  "sample_size"     = "1. Sample Size",
  "variance"        = "2. Variance",
  "stat_tests"      = "3. Stat. Evidence",
  "drift"           = "4. Drift Control",
  "overhead"        = "5. Overhead",
  "noise_model"     = "6. Noise Model",
  "reproducibility" = "7. Reproducibility",
  "neg_results"     = "8. Neg. Results"
)

criteria_cols <- c("sample_size", "variance", "stat_tests", "drift",
                   "overhead", "noise_model", "reproducibility", "neg_results")

# ── Load data ──
data_path <- file.path("reproduction", "data", "review_criteria.csv")
df <- read_csv(data_path, show_col_types = FALSE)

# Separate core review from foundational
df <- df %>%
  mutate(
    scope = case_when(
      category == "Foundational" ~ "foundational",
      TRUE ~ "core"
    )
  )

core_primary <- df %>%
  filter(scope == "core", category != "Reviews & Surveys")

cat(sprintf("Core primary papers (excl. reviews & foundational): %d\n", nrow(core_primary)))

# ── Compute percentages ──
criterion_stats <- core_primary %>%
  pivot_longer(cols = all_of(criteria_cols), names_to = "criterion", values_to = "rating") %>%
  filter(rating != "na") %>%
  group_by(criterion) %>%
  summarise(
    n_total = n(),
    pct_yes = sum(rating == "yes") / n_total * 100,
    pct_partial = sum(rating == "partial") / n_total * 100,
    pct_no = sum(rating == "no") / n_total * 100,
    .groups = "drop"
  ) %>%
  mutate(criterion = factor(criterion, levels = criteria_cols))

p1_data <- criterion_stats %>%
  pivot_longer(cols = c(pct_yes, pct_partial, pct_no),
               names_to = "rating_type", values_to = "pct") %>%
  mutate(
    rating_type = factor(rating_type,
      levels = c("pct_no", "pct_partial", "pct_yes"),
      labels = c("Missing", "Partial", "Adequate")
    ),
    criterion_label = criterion_names[as.character(criterion)]
  )

# ── Plot ──
p1 <- ggplot(p1_data, aes(x = reorder(criterion_label, -as.numeric(criterion)),
                            y = pct, fill = rating_type)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = c("Adequate" = LFD$teal,
                                "Partial"  = LFD$orange,
                                "Missing"  = LFD$red),
                    name = "Rating") +
  coord_flip() +
  labs(x = NULL, y = "Percentage of applicable papers") +
  scale_y_continuous(labels = latex_percent, limits = c(0, 100)) +
  theme_paper() +
  theme(legend.position = "top",
        legend.direction = "horizontal")

# ── Save ──
save_plot(p1, "review_compliance", width = COLWIDTH, height = 0.75 * COLWIDTH)
