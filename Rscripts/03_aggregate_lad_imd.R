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

lookup <- read_csv(
  "data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv",
  show_col_types = FALSE
)

# Aggregate LSOA IMD to LAD level
lad_imd <- lookup %>%
  group_by(lad_code, lad_name) %>%
  summarise(
    n_lsoa            = n(),
    mean_imd_score    = mean(imd_score,   na.rm = TRUE),
    median_imd_score  = median(imd_score, na.rm = TRUE),
    mean_lsoa_imd_decile = mean(imd_decile, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_imd_score)) %>%
  mutate(
    lad_imd_rank   = row_number(),
    lad_imd_decile = ntile(desc(mean_imd_score), 10)
  )

cat("LADs:", nrow(lad_imd), "\n")
print(table(lad_imd$lad_imd_decile))

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write_csv(lad_imd, "data/processed/LAD21_IMD_summary.csv")

# Attach STP/ITL geography using majority rule
# (some LADs span multiple STPs after 2021 boundary changes;
#  keep the STP/ITL supported by the most LSOAs within each LAD)
lad_geo_imd <- lookup %>%
  count(lad_code, lad_name,
        stp_code, stp_name,
        itl2_code, itl2_name,
        itl1_code, itl1_name) %>%
  group_by(lad_code, lad_name) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-n) %>%
  left_join(lad_imd, by = c("lad_code", "lad_name"))

cat("Rows:", nrow(lad_geo_imd), "\n")
cat("Distinct LADs:", n_distinct(lad_geo_imd$lad_code), "\n")

write_csv(lad_geo_imd, "data/processed/LAD21_geography_IMD_summary.csv")
