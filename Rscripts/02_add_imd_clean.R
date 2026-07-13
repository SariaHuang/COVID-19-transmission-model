# ==============================================================
# Script: 02_add_imd_clean.R
#
# Purpose: Join the geography lookup (from script 01) with IMD 2019
#          deprivation data at LSOA level.
#
# Input:   data/processed/lookup_LSOA11_STP_ITL2.csv (from 01)
#          data/raw/File_7_-_All_IoD2019_Scores_..._3.csv (IMD 2019)
#
# Output:  data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv
# ==============================================================

library(readr)
library(dplyr)

# Geography lookup (from script 01)
lookup <- read_csv(
  "data/processed/lookup_LSOA11_STP_ITL2.csv",
  show_col_types = FALSE
)

# IMD 2019
imd <- read_csv(
  "data/raw/File_7_-_All_IoD2019_Scores__Ranks__Deciles_and_Population_Denominators_3.csv",
  show_col_types = FALSE
) %>%
  transmute(
    lsoa_code  = `LSOA code (2011)`,
    imd_score  = `Index of Multiple Deprivation (IMD) Score`,
    imd_rank   = `Index of Multiple Deprivation (IMD) Rank (where 1 is most deprived)`,
    imd_decile = `Index of Multiple Deprivation (IMD) Decile (where 1 is most deprived 10% of LSOAs)`
  )

# Join
lookup <- lookup %>%
  left_join(imd, by = "lsoa_code")

# checks
cat("LSOAs:", n_distinct(lookup$lsoa_code), "\n")
cat("Missing IMD decile:", sum(is.na(lookup$imd_decile)), "\n")
cat("STPs:", n_distinct(lookup$stp_name), "\n")
cat("LADs:", n_distinct(lookup$lad_name), "\n")
cat("ITL2 regions:", n_distinct(lookup$itl2_name), "\n")
cat("ITL1 regions:", n_distinct(lookup$itl1_name), "\n")
cat("IMD deciles:", n_distinct(lookup$imd_decile), "\n")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write_csv(lookup, "data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv")
