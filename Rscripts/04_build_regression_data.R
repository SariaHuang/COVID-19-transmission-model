# ==============================================================
# Script: 04_build_regression_data.R
#
# Purpose: Build the final LAD-by-week regression dataset by
#          combining NHS hospital admissions, LAD-level IMD
#          deprivation/geography, OxCGRT government policy
#          stringency data, and ONS mid-year population
#          estimates.
#
# Inputs:
#   - NHS admissions data (pulled via covid19.nhs.data::get_admissions())
#   - data/processed/LAD21_geography_IMD_summary.csv (output of script 03)
#   - OxCGRT policy data (pulled directly from GitHub)
#   - data/raw/ukpopestimatesmid2020on2021geography.xls (ONS population estimates)
#
# Output:
#   - data/processed/regression_data.rds
#
# Note: `release_date` below is still TBD — update once the final
#       NHS data release date for this analysis is decided.
# ==============================================================

library(covid19.nhs.data)
library(dplyr)
library(lubridate)
library(readxl)

release_date <- as.Date("2022-03-03")

# NHS hospital admissions (LAD level)
adm_full <- get_admissions(
  keep_vars    = "new_adm",
  level        = "ltla",
  release_date = release_date
)

# LAD geography + IMD (from script 03)
lad_imd <- read.csv("data/processed/LAD21_geography_IMD_summary.csv")
cat("LADs:", nrow(lad_imd), "\n")

# Recode pre-2021 LAD codes to LAD21
# Buckinghamshire and Northamptonshire were reorganised in April 2021
lad_recode <- tibble::tribble(
  ~old_code,    ~new_code,
  "E07000004",  "E06000060",
  "E07000005",  "E06000060",
  "E07000006",  "E06000060",
  "E07000007",  "E06000060",
  "E07000150",  "E06000061",
  "E07000152",  "E06000061",
  "E07000153",  "E06000061",
  "E07000156",  "E06000061",
  "E07000151",  "E06000062",
  "E07000154",  "E06000062",
  "E07000155",  "E06000062"
)

adm_imd <- adm_full %>%
  filter(!is.na(geo_code)) %>%
  left_join(lad_recode, by = c("geo_code" = "old_code")) %>%
  mutate(geo_code = coalesce(new_code, geo_code)) %>%
  select(-new_code) %>%
  left_join(
    lad_imd %>% select(lad_code, itl1_code, itl1_name, lad_imd_decile),
    by = c("geo_code" = "lad_code")
  )

unmatched <- adm_imd %>%
  filter(is.na(lad_imd_decile)) %>%
  distinct(geo_code, geo_name)
cat("Unmatched areas:", nrow(unmatched), "\n")
if (nrow(unmatched) > 0) print(unmatched)

# OxCGRT policy stringency (England, weekly with 1-week lag)
oxcgrt <- read.csv(
  "https://raw.githubusercontent.com/OxCGRT/covid-policy-dataset/main/data/OxCGRT_compact_subnational_v1.csv"
)

policy_weekly <- oxcgrt %>%
  filter(CountryCode == "GBR", RegionCode == "UK_ENG") %>%
  mutate(date = as.Date(as.character(Date), format = "%Y%m%d")) %>%
  select(date, StringencyIndex_Average,
         C1M_School.closing, C2M_Workplace.closing,
         C6M_Stay.at.home.requirements) %>%
  filter(date >= as.Date("2020-08-01"), date <= release_date) %>%
  mutate(epiweek = floor_date(date, "week", week_start = 1)) %>%
  group_by(epiweek) %>%
  summarise(
    stringency_mean   = mean(StringencyIndex_Average, na.rm = TRUE),
    school_closing    = max(C1M_School.closing,           na.rm = TRUE),
    workplace_closing = max(C2M_Workplace.closing,        na.rm = TRUE),
    stay_at_home      = max(C6M_Stay.at.home.requirements, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(epiweek) %>%
  mutate(stringency_lag1 = lag(stringency_mean, 1))

# Aggregate admissions to weekly LAD level, join policy and population
pop_lad <- read_excel(
  "data/raw/ukpopestimatesmid2020on2021geography.xls",
  sheet = "MYE2 - Persons",
  skip  = 7
) %>%
  rename(geo_code = Code, population = `All ages`) %>%
  filter(
    startsWith(geo_code, "E"),
    !Geography %in% c("Country", "Region", "Metropolitan County",
                      "County", "Inner London", "Outer London")
  ) %>%
  select(geo_code, population)

cat("LADs in population data:", nrow(pop_lad), "\n")

regression_data <- adm_imd %>%
  filter(!is.na(lad_imd_decile)) %>%
  mutate(epiweek = floor_date(date, "week", week_start = 1)) %>%
  group_by(geo_code, geo_name, itl1_code, itl1_name,
           lad_imd_decile, epiweek) %>%
  summarise(hosp_admissions = sum(admissions, na.rm = TRUE),
            .groups = "drop") %>%
  left_join(policy_weekly, by = "epiweek") %>%
  left_join(pop_lad,       by = "geo_code")

cat("Missing population:", sum(is.na(regression_data$population)), "\n")
cat("Total rows:", nrow(regression_data), "\n")
cat("Weeks:", n_distinct(regression_data$epiweek), "\n")
cat("LADs:", n_distinct(regression_data$geo_code), "\n")
cat("IMD deciles:", sort(unique(regression_data$lad_imd_decile)), "\n")
cat("Missing stringency:", sum(is.na(regression_data$stringency_mean)), "\n")

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
saveRDS(regression_data, "data/processed/regression_data.rds")
