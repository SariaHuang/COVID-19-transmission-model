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
library(readr)

# LAD geography + IMD decile (from script 03)
lad_imd <- read.csv("data/processed/LAD21_geography_IMD_summary.csv") %>%
  select(lad_code, lad_name, itl1_name, lad_imd_decile)

# ONS mid-year 2020 population by single year of age
pop_raw <- read_excel(
  "data/raw/ukpopestimatesmid2020on2021geography.xls",
  sheet = "MYE2 - Persons",
  skip  = 7
)

# Single-year-of-age columns
age_cols <- names(pop_raw)[grepl("^[0-9]+\\+?$", names(pop_raw))]
cat("Age columns detected:", length(age_cols), "\n")

# Keep English LADs, reshape to long
pop_long <- pop_raw %>%
  rename(lad_code = Code) %>%
  filter(
    startsWith(lad_code, "E"),
    !Geography %in% c("Country", "Region", "Metropolitan County",
                      "County", "Inner London", "Outer London")
  ) %>%
  select(lad_code, all_of(age_cols)) %>%
  pivot_longer(cols = all_of(age_cols),
               names_to = "age_single",
               values_to = "population") %>%
  mutate(age_single = as.integer(gsub("\\+", "", age_single)))

cat("Distinct LADs:", n_distinct(pop_long$lad_code), "(should be 309)\n")

# Bin into 5-year age bands matching Goodfellow et al. (2024) / Knock et al. (2021)
age_breaks <- c(seq(0, 80, by = 5), Inf)
age_labels <- c(paste(seq(0, 75, by = 5), seq(4, 79, by = 5), sep = "-"), "80+")

pop_banded <- pop_long %>%
  mutate(age_band = cut(age_single, breaks = age_breaks, labels = age_labels,
                        right = FALSE, include.lowest = TRUE)) %>%
  group_by(lad_code, age_band) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop")

# Cross-check: banded total should match 'All ages' column
total_banded   <- sum(pop_banded$population)
total_all_ages <- pop_raw %>%
  rename(lad_code = Code) %>%
  filter(startsWith(lad_code, "E"),
         !Geography %in% c("Country", "Region", "Metropolitan County",
                           "County", "Inner London", "Outer London")) %>%
  pull(`All ages`) %>% sum(na.rm = TRUE)

cat("Banded total:  ", total_banded, "\n")
cat("All ages total:", total_all_ages, "\n")

# Aggregate to age x IMD decile x ITL1 region
pop_strata <- pop_banded %>%
  left_join(lad_imd, by = "lad_code") %>%
  filter(!is.na(lad_imd_decile)) %>%
  group_by(age_band, lad_imd_decile, itl1_name) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop")

# checks
cat("Total population:", sum(pop_strata$population), "\n")
cat("Age bands:", n_distinct(pop_strata$age_band), "\n")
cat("IMD deciles:", n_distinct(pop_strata$lad_imd_decile), "\n")
cat("ITL1 regions:", n_distinct(pop_strata$itl1_name), "\n")

n_expected <- length(age_labels) * 10 * n_distinct(pop_strata$itl1_name)
cat("Expected cells:", n_expected, "| Actual:", nrow(pop_strata),
    "| Missing:", n_expected - nrow(pop_strata), "\n")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
write_csv(pop_strata, "data/processed/population_age_imd_region.csv")
