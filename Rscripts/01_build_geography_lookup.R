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
#   data/processed/lookup_LSOA11_STP_ITL2.csv (LSOA -> STP -> LAD -> ITL1/2/3)
# ==============================================================

library(readr)
library(dplyr)

# LSOA to STP/LAD (England)
lsoa_stp_lad <- read_csv(
  "data/raw/LSOA11_CCG21_STP21_LAD21_EN_LU_ae1a442ee397483cab1a31f2e7b24029_997658846778534204.csv",
  show_col_types = FALSE
) %>%
  transmute(
    lsoa_code = LSOA11CD,
    lsoa_name = LSOA11NM,
    stp_code  = STP21CD,
    stp_name  = STP21NM,
    lad_code  = LAD21CD,
    lad_name  = LAD21NM
  )

# LAD to ITL regions (UK)
lad_itl <- read_csv(
  "data/raw/Local_Authority_District_(April_2021)_to_LAU1_to_ITL3_to_ITL2_to_ITL1_(January_2021)_Lookup_in_United_Kingdom.csv",
  show_col_types = FALSE
) %>%
  transmute(
    lad_code  = LAD21CD,
    itl2_code = ITL221CD,
    itl2_name = ITL221NM,
    itl1_code = ITL121CD,
    itl1_name = ITL121NM
  ) %>%
  distinct(lad_code, .keep_all = TRUE)  

# Join: LSOA -> STP -> LAD -> ITL
lookup <- lsoa_stp_lad %>%
  left_join(lad_itl, by = "lad_code")

# check
cat("LSOAs:", n_distinct(lookup$lsoa_code), "\n")
cat("STPs:", n_distinct(lookup$stp_name), "\n")
cat("ITL2 regions:", n_distinct(lookup$itl2_name), "\n")
cat("Missing ITL2:", sum(is.na(lookup$itl2_name)), "\n")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write_csv(lookup, "data/processed/lookup_LSOA11_STP_ITL2.csv")
