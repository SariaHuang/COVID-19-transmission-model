rm(list = ls())
# ==============================================================
# Script: 05_build_population_strata.R
#
# Purpose: Build the joint population matrix N_{a,i,r} -- population
#          stratified by age band (a), LAD-level IMD decile (i),
#          and NHS region/ITL1 (r) -- needed for the dynamic
#          transmission model's population stratification term:
#            pop_prop_{a,i,r} = N_{a,i,r} / sum_a N_{a,i,r}
# Inputs:
#   data/processed/LAD21_geography_IMD_summary.csv (script 03)
#   data/raw/ukpopestimatesmid2020on2021geography.xls
#
# Output:
#   data/processed/population_age_imd_region.csv
#     (long format: age_band, lad_imd_decile, itl1_name, population)
# ==============================================================

library(readxl)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 1. Read LAD-level geography + IMD decile (from script 03)
# ------------------------------------------------------------
lad_imd <- read.csv("data/processed/LAD21_geography_IMD_summary.csv") %>%
  select(lad_code, lad_name, itl1_name, lad_imd_decile)

# ------------------------------------------------------------
# 2. Read population by single year of age
# ------------------------------------------------------------
pop_raw <- read_excel(
  "data/raw/ukpopestimatesmid2020on2021geography.xls",
  sheet = "MYE2 - Persons",
  skip  = 7
)

# Detect single-year-of-age columns (typically named "0","1",...,"90+")
# CHECK THIS OUTPUT -- if 0 columns are detected, open the Excel file
# and look at the actual column names, then adjust the regex below.
age_cols <- names(pop_raw)[grepl("^[0-9]+\\+?$", names(pop_raw))]
cat("Detected", length(age_cols), "single-year-of-age columns\n")
if (length(age_cols) > 0) {
  cat("First few:", paste(head(age_cols), collapse = ", "), "\n")
  cat("Last few:", paste(tail(age_cols), collapse = ", "), "\n")
} else {
  cat("WARNING: no age columns detected. Run names(pop_raw) and check manually.\n")
}

# ------------------------------------------------------------
# 3. Keep English LADs only, reshape to long format
# ------------------------------------------------------------
pop_long <- pop_raw %>%
  rename(lad_code = Code) %>%
  filter(
    startsWith(lad_code, "E"),
    !Geography %in% c("Country", "Region",
                      "Metropolitan County",
                      "County",
                      "Inner London",
                      "Outer London")
  ) %>%
  select(lad_code, all_of(age_cols)) %>%
  pivot_longer(
    cols = all_of(age_cols),
    names_to = "age_single",
    values_to = "population"
  ) %>%
  mutate(
    age_single = as.integer(gsub("\\+", "", age_single))
  )

# ------------------------------------------------------------
# 4. Bin into 5-year age bands (Knock 2021 / Imperial COVID model
#    convention: 0-4, 5-9, ..., 75-79, 80+)
#    TODO: confirm these cut points match the age bands used in
#    Goodfellow's clinical fraction / IHR / IFR tables exactly --
#    if his bands differ, re-bin to match his, not the other way
#    round, since his pi_a / IHR_a / IFR_a values are fixed inputs.
# ------------------------------------------------------------
age_breaks <- c(seq(0, 80, by = 5), Inf)
age_labels <- c(paste(seq(0, 75, by = 5), seq(4, 79, by = 5), sep = "-"), "80+")

pop_long <- pop_long %>%
  mutate(
    age_band = cut(age_single, breaks = age_breaks, labels = age_labels,
                   right = TRUE, include.lowest = TRUE)
  )

# --- DIAGNOSTIC: isolate where any population inflation comes from ---
# Check 1: how many distinct LADs survived the filter? Should be exactly 309.
n_lads_filtered <- pop_raw %>%
  rename(lad_code = Code) %>%
  filter(
    startsWith(lad_code, "E"),
    !Geography %in% c("Country", "Region",
                      "Metropolitan County",
                      "County",
                      "Inner London",
                      "Outer London")
  ) %>%
  pull(lad_code) %>%
  n_distinct()
cat("DIAGNOSTIC -- distinct LADs after filtering:", n_lads_filtered,
    "(should be 309)\n")

# Check 2: sum of the 'All ages' column directly (bypasses age-band pivot
# entirely) for the same filtered rows, as an independent cross-check.
if ("All ages" %in% names(pop_raw)) {
  total_all_ages <- pop_raw %>%
    rename(lad_code = Code) %>%
    filter(
      startsWith(lad_code, "E"),
      !Geography %in% c("Country", "Region",
                        "Metropolitan County",
                        "County",
                        "Inner London",
                        "Outer London")
    ) %>%
    pull(`All ages`) %>%
    sum(na.rm = TRUE)
  cat("DIAGNOSTIC -- sum of 'All ages' column directly:", total_all_ages, "\n")
  cat("DIAGNOSTIC -- sum from age_cols pivot (should match):",
      sum(pop_long$population, na.rm = TRUE), "\n")
} else {
  cat("DIAGNOSTIC -- no 'All ages' column found; cannot cross-check.\n")
}
# --- END DIAGNOSTIC ---

pop_banded <- pop_long %>%
  group_by(lad_code, age_band) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 5. Join with geography + IMD decile, aggregate to N_{a,i,r}
# ------------------------------------------------------------
pop_strata <- pop_banded %>%
  left_join(lad_imd, by = "lad_code") %>%
  filter(!is.na(lad_imd_decile)) %>%
  group_by(age_band, lad_imd_decile, itl1_name) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 6. Sanity checks
# ------------------------------------------------------------
cat("\n--- Sanity checks ---\n")
cat("Total population (banded):", sum(pop_strata$population), "\n")
cat("Number of age bands:", n_distinct(pop_strata$age_band), "\n")
cat("Number of IMD deciles:", n_distinct(pop_strata$lad_imd_decile), "\n")
cat("Number of regions:", n_distinct(pop_strata$itl1_name), "\n")

expected_cells <- length(age_labels) * 10 * n_distinct(pop_strata$itl1_name)
cat("Expected cells (age x decile x region):", expected_cells, "\n")
cat("Actual cells with data:", nrow(pop_strata), "\n")
cat("Missing cells:", expected_cells - nrow(pop_strata), "\n")

# ------------------------------------------------------------
# 7. Save
# ------------------------------------------------------------
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write.csv(pop_strata, "data/processed/population_age_imd_region.csv", row.names = FALSE)

