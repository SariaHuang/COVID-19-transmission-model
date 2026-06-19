# ==============================================================
# Script: 01_build_geography_lookup.R
#
# Purpose:
#   Build a lookup table linking LSOA -> STP -> LAD -> ITL
#   geography codes, using ONS/NHS geography reference files.
#
# Inputs:
#   data/raw/LSOA11_CCG21_STP21_LAD21_EN_LU_ae1a442ee397483cab1a31f2e7b24029_997658846778534204.csv
#   data/raw/Local_Authority_District_(April_2021)_to_LAU1_to_ITL3_to_ITL2_to_ITL1_(January_2021)_Lookup_in_United_Kingdom.csv
#
# Output:
#   data/processed/lookup_LSOA11_STP_ITL2.csv
# ==============================================================

library(readr)
library(dplyr)

# ------------------------------------------------------------
# Read geography lookup files
# ------------------------------------------------------------
lsoa_stp_lad <- read_csv(
  "data/raw/LSOA11_CCG21_STP21_LAD21_EN_LU_ae1a442ee397483cab1a31f2e7b24029_997658846778534204.csv",
  show_col_types = FALSE
)

lad_itl <- read_csv(
  "data/raw/Local_Authority_District_(April_2021)_to_LAU1_to_ITL3_to_ITL2_to_ITL1_(January_2021)_Lookup_in_United_Kingdom.csv",
  show_col_types = FALSE
)

# ------------------------------------------------------------
# Clean LSOA -> STP/LAD lookup
# ------------------------------------------------------------
lookup_lsoa_stp_lad <- lsoa_stp_lad %>%
  transmute(
    lsoa_code = LSOA11CD,
    lsoa_name = LSOA11NM,
    stp_code  = STP21CD,
    stp_name  = STP21NM,
    lad_code  = LAD21CD,
    lad_name  = LAD21NM
  )

# ------------------------------------------------------------
# Clean LAD -> ITL lookup
# ------------------------------------------------------------
lookup_lad_itl <- lad_itl %>%
  transmute(
    lad_code     = LAD21CD,
    lad_name_itl = LAD21NM,
    itl3_code    = ITL321CD,
    itl3_name    = ITL321NM,
    itl2_code    = ITL221CD,
    itl2_name    = ITL221NM,
    itl1_code    = ITL121CD,
    itl1_name    = ITL121NM
  )

# ------------------------------------------------------------
# Join to create LSOA -> STP -> ITL lookup
# ------------------------------------------------------------
lookup_lsoa_stp_itl <- lookup_lsoa_stp_lad %>%
  left_join(lookup_lad_itl, by = "lad_code")

# ------------------------------------------------------------
# Check result
# ------------------------------------------------------------
cat("Unique LSOAs:", n_distinct(lookup_lsoa_stp_itl$lsoa_code), "\n")
cat("Unique STPs:", n_distinct(lookup_lsoa_stp_itl$stp_name), "\n")
cat("Unique ITL2:", n_distinct(lookup_lsoa_stp_itl$itl2_name), "\n")
cat("Missing ITL2 rows:", sum(is.na(lookup_lsoa_stp_itl$itl2_name)), "\n")

lookup_lsoa_stp_itl %>%
  select(lsoa_code, stp_name, lad_name, itl2_name, itl1_name) %>%
  head()

# ------------------------------------------------------------
# Save result
# ------------------------------------------------------------
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

write_csv(
  lookup_lsoa_stp_itl,
  "data/processed/lookup_LSOA11_STP_ITL2.csv"
)
