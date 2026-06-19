rm(list = ls())
# ==============================================================
# Script: 11_regression_stringency.R
#
# Purpose: Replace factor(epiweek) (saturated week fixed effects)
#          with a natural spline time trend, ns(week_num, df),
#          to free up identifying variation for stringency_lag1
#          while still controlling for the underlying epidemic
#          time trend.
#
# Background: factor(epiweek) is saturated -- it fully absorbs
#          all variation that occurs at the national-weekly level,
#          including stringency (which only varies by week, not
#          by region, since OxCGRT has no sub-national England
#          data). This causes perfect collinearity. Replacing the
#          saturated week FE with a lower-df spline frees up some
#          of that variation for stringency to be estimated.
#
#          Stringency only varies over time at
#          the national level, and so does epidemic severity --
#          both are driven by the same underlying epidemic
#          dynamics. No amount of spline flexibility identifies a
#          causal policy effect from this data alone; the
#          stringency_lag1 coefficient should be reported as a
#          descriptive association, not a causal estimate.
#
# Input: data/processed/regression_data.rds (script 04)
# ==============================================================

library(MASS)
library(dplyr)
library(broom)
library(splines)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load and prepare data
# ------------------------------------------------------------
reg_df <- readRDS("data/processed/regression_data.rds") %>%
  filter(!is.na(stringency_lag1), population > 0) %>%
  mutate(week_num = as.integer(factor(epiweek)))

cat("Rows after filtering:", nrow(reg_df), "\n")
cat("Number of distinct weeks:", n_distinct(reg_df$week_num), "\n")

# ------------------------------------------------------------
# 2. Compare AIC across a range of spline degrees of freedom
#    Lower df = smoother trend = more room for stringency, but
#    risk of under-fitting real epidemic waves.
#    Higher df = closer to saturated week FE = less room for
#    stringency, but better captures true epidemic shape.
# ------------------------------------------------------------
df_candidates <- c(4, 6, 8, 10, 12, 16)

fit_with_df <- function(k) {
  glm.nb(
    hosp_admissions ~ factor(lad_imd_decile) * factor(itl1_name) +
      ns(week_num, df = k) +
      stringency_lag1 +
      offset(log(population)),
    data = reg_df
  )
}

cat("\nFitting models across spline df candidates (this may take a moment)...\n")
models_by_df <- lapply(df_candidates, function(k) {
  cat("  df =", k, "\n")
  fit_with_df(k)
})
names(models_by_df) <- paste0("df_", df_candidates)

aic_table <- data.frame(
  df  = df_candidates,
  AIC = sapply(models_by_df, AIC)
)
cat("\n--- AIC by spline df ---\n")
print(aic_table)
cat("\nLower AIC = better fit. Check if AIC is still improving at the\n")
cat("highest df tested -- if so, consider testing higher df values too.\n")

best_df <- aic_table$df[which.min(aic_table$AIC)]
cat("\nBest-fitting df by AIC:", best_df, "\n")

# ------------------------------------------------------------
# 3. Stability check: does the stringency coefficient change
#    a lot depending on df? If yes, the estimate is fragile and
#    this should be reported as a limitation.
# ------------------------------------------------------------
stringency_by_df <- lapply(names(models_by_df), function(nm) {
  m <- models_by_df[[nm]]
  tidy(m, exponentiate = TRUE, conf.int = FALSE) %>%
    filter(term == "stringency_lag1") %>%
    mutate(
      df = as.integer(gsub("df_", "", nm)),
      conf.low  = exp(log(estimate) - 1.96 * std.error),
      conf.high = exp(log(estimate) + 1.96 * std.error)
    )
}) %>% bind_rows()

cat("\n--- stringency_lag1 IRR across different spline df ---\n")
print(stringency_by_df %>% dplyr::select(df, estimate, conf.low, conf.high))
cat("\nIf the IRR estimate or its significance flips around as df changes,\n")
cat("treat this coefficient as unstable -- report this instability rather\n")
cat("than picking the df that gives the 'nicest' result.\n")

# ------------------------------------------------------------
# 4. Diagnostic plot: does the chosen spline actually capture
#    the real shape of the epidemic, or does it over-smooth?
#    Compare the spline-implied time trend against raw weekly
#    admissions totals.
# ------------------------------------------------------------
raw_weekly <- reg_df %>%
  group_by(week_num) %>%
  summarise(total_admissions = sum(hosp_admissions), .groups = "drop")

# Predict the spline component only, holding other covariates at
# a reference level, to visualise the shape of the fitted time trend
best_model <- models_by_df[[paste0("df_", best_df)]]

pred_grid <- reg_df %>%
  distinct(week_num) %>%
  arrange(week_num) %>%
  mutate(
    lad_imd_decile = reg_df$lad_imd_decile[1],
    itl1_name      = reg_df$itl1_name[1],
    stringency_lag1 = mean(reg_df$stringency_lag1, na.rm = TRUE),
    population      = mean(reg_df$population, na.rm = TRUE)
  )
pred_grid$spline_fit <- predict(best_model, newdata = pred_grid, type = "response")

p_diag <- ggplot() +
  geom_col(data = raw_weekly, aes(x = week_num, y = total_admissions),
           fill = "grey80", width = 0.8) +
  geom_line(data = pred_grid, aes(x = week_num, y = spline_fit * 50),
            colour = "firebrick", linewidth = 1) +
  labs(
    title    = paste0("Spline time trend (df=", best_df, ") vs raw weekly admissions"),
    subtitle = "Red line: model-implied time trend (rescaled for visibility) | Grey bars: raw weekly totals",
    x = "Week number", y = "Total weekly admissions (raw)"
  ) +
  theme_minimal()

print(p_diag)
dir.create("output", showWarnings = FALSE)
ggsave("output/spline_diagnostic.png", p_diag, width = 10, height = 6, dpi = 200)

# ------------------------------------------------------------
# 5. Final model summary (using best_df, or override manually
#    below if you decide a different df is more defensible)
# ------------------------------------------------------------
cat("\n--- Final model summary (df =", best_df, ") ---\n")
print(summary(best_model)$coefficients["stringency_lag1", ])

cat("\n==============================================================\n")
cat("REMINDER: this coefficient is a descriptive association, not a\n")
cat("causal estimate. Stringency varies only at the national-weekly\n")
cat("level in our data (OxCGRT has no England sub-national policy\n")
cat("granularity), so it is structurally confounded with national\n")
cat("epidemic dynamics. Report accordingly.\n")
cat("==============================================================\n")
