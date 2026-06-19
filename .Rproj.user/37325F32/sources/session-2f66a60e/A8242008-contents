# ==============================================================
# Script: 02_add_imd_clean.R
#
# Purpose: Join the geography lookup (from 01) with IMD 2019
#          deprivation data to produce a clean LSOA-level table
#          with geography + deprivation score/rank/decile.
#
# Input:   data/processed/lookup_LSOA11_STP_ITL2.csv (from 01)
#          data/raw/File_7_-_All_IoD2019_Scores_..._3.csv
#
# Output:  data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv
# ==============================================================

library(readr)
library(dplyr)

# 1. Read geography lookup, drop duplicate LSOA rows
lookup_itl2 <- read_csv(
  "data/processed/lookup_LSOA11_STP_ITL2.csv",
  show_col_types = FALSE
)

lookup_itl2_clean <- lookup_itl2 %>%
  distinct(lsoa_code, .keep_all = TRUE)

cat("Rows before dedup:", nrow(lookup_itl2), "\n")
cat("Rows after dedup:", nrow(lookup_itl2_clean), "\n")
cat("Unique LSOAs:", n_distinct(lookup_itl2_clean$lsoa_code), "\n")

# 2. Read IMD 2019 data, keep relevant columns
imd <- read_csv(
  "data/raw/File_7_-_All_IoD2019_Scores__Ranks__Deciles_and_Population_Denominators_3.csv",
  show_col_types = FALSE
)

imd_clean <- imd %>%
  transmute(
    lsoa_code  = `LSOA code (2011)`,
    imd_score  = `Index of Multiple Deprivation (IMD) Score`,
    imd_rank   = `Index of Multiple Deprivation (IMD) Rank (where 1 is most deprived)`,
    imd_decile = `Index of Multiple Deprivation (IMD) Decile (where 1 is most deprived 10% of LSOAs)`
  )

# 3. Join geography lookup with IMD
lookup_lsoa_stp_itl_imd_clean <- lookup_itl2_clean %>%
  left_join(imd_clean, by = "lsoa_code")

# 4. Sanity checks
cat("Rows:", nrow(lookup_lsoa_stp_itl_imd_clean), "\n")
cat("Unique LSOAs:", n_distinct(lookup_lsoa_stp_itl_imd_clean$lsoa_code), "\n")
cat("Missing IMD:", sum(is.na(lookup_lsoa_stp_itl_imd_clean$imd_decile)), "\n")
cat("Unique STPs:", n_distinct(lookup_lsoa_stp_itl_imd_clean$stp_name), "\n")
cat("Unique LAD/LTLA:", n_distinct(lookup_lsoa_stp_itl_imd_clean$lad_name), "\n")
cat("Unique ITL3:", n_distinct(lookup_lsoa_stp_itl_imd_clean$itl3_name), "\n")
cat("Unique ITL2:", n_distinct(lookup_lsoa_stp_itl_imd_clean$itl2_name), "\n")
cat("Unique Region/ITL1:", n_distinct(lookup_lsoa_stp_itl_imd_clean$itl1_name), "\n")
cat("Unique IMD deciles:", n_distinct(lookup_lsoa_stp_itl_imd_clean$imd_decile), "\n")

# Confirm no duplicate LSOAs remain
lookup_lsoa_stp_itl_imd_clean %>%
  count(lsoa_code) %>%
  filter(n > 1)

# 5. Save
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(
  lookup_lsoa_stp_itl_imd_clean,
  "data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv"
)
