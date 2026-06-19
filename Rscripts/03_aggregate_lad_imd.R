# ==============================================================
# Script: 03_aggregate_lad_imd.R
#
# Purpose: Aggregate LSOA-level IMD 2019 deprivation scores up
#          to LAD level, then attach geography (STP/ITL) info
#          to produce a LAD-level table with both deprivation
#          summary stats and geography.
#
# Input:   data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv
#          (output of script 02)
#
# Output:  data/processed/LAD21_IMD_summary.csv
#          data/processed/LAD21_geography_IMD_summary.csv
# ==============================================================

library(readr)
library(dplyr)

# 1. Read the LSOA-level geography + IMD lookup
lookup_lsoa_imd <- read_csv(
  "data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv",
  show_col_types = FALSE
)

# 2. Aggregate LSOA-level IMD to LAD-level summary stats,
#    then rank and decile LADs by mean IMD score
lad21_imd_summary <- lookup_lsoa_imd %>%
  group_by(lad_code, lad_name) %>%
  summarise(
    n_lsoa = n(),
    mean_imd_score = mean(imd_score, na.rm = TRUE),
    median_imd_score = median(imd_score, na.rm = TRUE),
    mean_lsoa_imd_decile = mean(imd_decile, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_imd_score)) %>%
  mutate(
    lad_imd_rank = row_number(),
    lad_imd_decile = ntile(desc(mean_imd_score), 10)
  )

# 3. Check output
cat("Number of LAD21 areas:", nrow(lad21_imd_summary), "\n")
cat("LAD-level IMD deciles:\n")
print(table(lad21_imd_summary$lad_imd_decile))
print(head(lad21_imd_summary))

# 4. Save LAD-level IMD summary
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(
  lad21_imd_summary,
  "data/processed/LAD21_IMD_summary.csv"
)

# 5. Attach geography (STP/ITL) info to the LAD-level IMD summary
#    NOTE: a small number of LADs (mostly those created by the 2021
#    boundary reorganisation, e.g. Buckinghamshire, North/West
#    Northamptonshire) span more than one STP, because STP boundaries
#    were not redrawn to match the new LAD boundaries. To guarantee
#    exactly one row per LAD, we keep the STP/ITL combination
#    supported by the most LSOAs within that LAD (majority rule).
lad21_geography_imd <- lookup_lsoa_imd %>%
  count(
    lad_code,
    lad_name,
    stp_code,
    stp_name,
    itl3_code,
    itl3_name,
    itl2_code,
    itl2_name,
    itl1_code,
    itl1_name
  ) %>%
  group_by(lad_code, lad_name) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-n) %>%
  left_join(
    lad21_imd_summary,
    by = c("lad_code", "lad_name")
  )

# Confirm the fix: should be exactly 309 rows, 309 distinct LADs
cat("Final geography+IMD table rows:", nrow(lad21_geography_imd), "\n")
cat("Distinct LADs:", n_distinct(lad21_geography_imd$lad_code), "\n")

write_csv(
  lad21_geography_imd,
  "data/processed/LAD21_geography_IMD_summary.csv"
)
