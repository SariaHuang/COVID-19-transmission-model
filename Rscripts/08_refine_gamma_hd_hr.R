# ==============================================================
# Script: 08_refine_gamma_hd_hr.R
#
# Purpose: Replace the provisional 80/20-weighted gamma_hd, gamma_hr
#          placeholders with age-varying values, properly aggregated
#          from Knock et al. (2021)'s hospital pathway structure
#          (Table S2: durations; Table S6: region-level absolute
#          branch probabilities; Table S8: age-scaling factors).
#
# METHOD:
#   Knock et al. (2021) do not report a single England-wide absolute
#   probability for each branch point -- Table S6 gives region-specific
#   posterior means (7 NHS regions), and Table S8 gives age-scaling
#   factors (relative to the age group with the highest value, not
#   absolute). There is no England-wide number to look up directly.
#
#   Knock's own method for aggregating region-level estimates into an
#   England-wide figure is to weight by each region's share of total
#   attack rate (their model's own output). We don't have access to
#   Knock's fitted attack rates, so instead we weight by each region's
#   REAL share of total observed COVID hospital admissions (from
#   data/processed/regression_data.rds, i.e. actual NHS England data,
#   not another model's estimate). This follows the same logic --
#   weight by share of burden -- using directly observed data instead
#   of a model-derived proxy.
#
# Pathway structure (Knock Figure S4 / S11):
#   Hospitalised patient either:
#     (a) NOT triaged to ICU (prob 1 - p_ICU(a)):
#           -> general ward -> dies (prob p_HD(a), duration 10.3d)
#                            or recovers (prob 1-p_HD(a), duration 10.7d)
#     (b) triaged to ICU (prob p_ICU(a)):
#           -> ICU_pre (2.5d) -> then:
#              - dies directly in ICU (prob p_ICUD(a), duration 11.8d)
#              - OR survives ICU, reaches stepdown (prob 1-p_ICUD(a)):
#                   -> dies in stepdown (prob p_WD(a), duration 7.0+8.1=15.1d)
#                   -> OR recovers in stepdown (prob 1-p_WD(a), duration 15.6+12.2=27.8d)
#
# Output:
#   data/parameters/gamma_hd_hr_by_age.csv
#   -- age-specific gamma_hd, gamma_hr, and P(death|hosp) columns
# ==============================================================

library(dplyr)
library(readr)

age_levels_gf <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                   "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                   "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                   "70 to 74","75+")

# ==============================================================
# SECTION 1 -- Table S8: age-scaling factors (relative, max = 1)
#   Exact values transcribed from Knock et al. (2021) supplementary
#   Table S8. These scale each branch probability RELATIVE to the
#   age group with the highest probability -- they are not absolute.
#   Knock's 17 age bands are: 0-5,5-10,...,75-80,80+ (5-year bands),
#   which we map onto Goodfellow's 17 bands below.
# ==============================================================
table_s8 <- tribble(
  ~age_band_knock, ~scale_hosp, ~scale_icu, ~scale_gen_death, ~scale_icu_death, ~scale_stepdown_death,
  "0-5",   0.039, 0.243, 0.039, 0.282, 0.091,
  "5-10",  0.001, 0.289, 0.037, 0.286, 0.083,
  "10-15", 0.006, 0.338, 0.035, 0.291, 0.077,
  "15-20", 0.009, 0.389, 0.035, 0.299, 0.074,
  "20-25", 0.026, 0.443, 0.036, 0.310, 0.074,
  "25-30", 0.040, 0.503, 0.039, 0.328, 0.076,
  "30-35", 0.042, 0.570, 0.045, 0.353, 0.080,
  "35-40", 0.045, 0.653, 0.055, 0.390, 0.086,
  "40-45", 0.050, 0.756, 0.074, 0.446, 0.093,
  "45-50", 0.074, 0.866, 0.107, 0.520, 0.102,
  "50-55", 0.138, 0.954, 0.157, 0.604, 0.117,
  "55-60", 0.198, 1.000, 0.238, 0.705, 0.148,
  "60-65", 0.247, 0.972, 0.353, 0.806, 0.211,
  "65-70", 0.414, 0.854, 0.502, 0.899, 0.332,
  "70-75", 0.638, 0.645, 0.675, 0.969, 0.526,
  "75-80", 1.000, 0.402, 0.832, 1.000, 0.753,
  "80+",   0.873, 0.107, 1.000, 0.918, 1.000
)

# Map Knock's 5-year bands (0-5,...,80+) onto Goodfellow's 17 bands
# (Under 1, 1-4, 5-9, ..., 75+). Knock's "0-5" is duplicated for
# Goodfellow's "Under 1" and "1 to 4" (same approximation used in
# script 07 for IHR/IFR). Knock's "75-80" and "80+" are averaged
# into Goodfellow's "75+".
row_0_5   <- table_s8 %>% filter(age_band_knock == "0-5")
row_75up  <- table_s8 %>% filter(age_band_knock %in% c("75-80","80+")) %>%
  summarise(across(starts_with("scale_"), mean))
rows_mid  <- table_s8 %>% filter(!age_band_knock %in% c("0-5","75-80","80+"))

table_s8_gf <- bind_rows(
  row_0_5 %>% dplyr::select(-age_band_knock),  # Under 1
  row_0_5 %>% dplyr::select(-age_band_knock),  # 1 to 4
  rows_mid %>% dplyr::select(-age_band_knock), # 5 to 9 ... 70 to 74
  row_75up                                     # 75+
) %>%
  mutate(age_band_gf = age_levels_gf) %>%
  dplyr::select(age_band_gf, everything())

cat("Table S8 age-scaling factors mapped to Goodfellow's 17 age bands:\n")
print(table_s8_gf)

# ==============================================================
# SECTION 2 -- Table S6: England-representative absolute
#   probabilities (max-scaling-group value), weighted by each
#   region's REAL observed share of total COVID hospital admissions
#   (from data/processed/regression_data.rds, script 04).
#
#   This follows Knock's logic (weight by each region's share of
#   total burden) but uses real observed NHS admissions data as the
#   weight, rather than Knock's own model-derived attack rate
#   estimates (which we don't have access to). Arguably more
#   defensible: it weights by directly observed burden rather than
#   another model's output.
#
#   Values transcribed from Knock et al. (2021) supplementary Table
#   S6, posterior means, one column per region:
#   NW, NEY, MID, EE, LON, SW, SE
# ==============================================================
table_s6_regional <- tribble(
  ~parameter,  ~NW,   ~NEY,  ~MID,  ~EE,   ~LON,  ~SW,   ~SE,
  "p_ICU_max", 0.14,  0.13,  0.18,  0.28,  0.27,  0.17,  0.23,
  "p_HD_max",  0.41,  0.44,  0.42,  0.58,  0.33,  0.36,  0.38,
  "p_ICUD_max",0.68,  0.72,  0.71,  0.70,  0.63,  0.68,  0.72,
  "p_WD_max",  0.39,  0.40,  0.37,  0.41,  0.32,  0.40,  0.38
)

# --- Compute real regional weights from NHS admissions data ---
regression_data <- readRDS("data/processed/regression_data.rds")

# Map our itl1_name values to Knock's 7-region abbreviations.
# Adjust the left-hand names below if your itl1_name values differ --
# check with unique(regression_data$itl1_name) first.
region_name_map <- c(
  "North West (England)"          = "NW",
  "North East (England)"          = "NEY",   # Knock combines North East + Yorkshire
  "Yorkshire and The Humber"      = "NEY",   # into a single "NEY" region
  "West Midlands (England)"       = "MID",
  "East Midlands (England)"       = "MID",   # Knock combines both Midlands into "MID"
  "East"                          = "EE",
  "London"                        = "LON",
  "South West (England)"          = "SW",
  "South East (England)"          = "SE"
)

region_weights <- regression_data %>%
  filter(!is.na(itl1_name), !is.na(hosp_admissions)) %>%
  mutate(knock_region = region_name_map[itl1_name]) %>%
  filter(!is.na(knock_region)) %>%
  group_by(knock_region) %>%
  summarise(total_admissions = sum(hosp_admissions, na.rm = TRUE), .groups = "drop") %>%
  mutate(weight = total_admissions / sum(total_admissions))

cat("\n--- Region weights from real NHS admissions data ---\n")
print(region_weights)
cat("(Weights should sum to 1; check that all 7 Knock regions are present)\n")

# Reshape weights into a named vector matching table_s6_regional's columns
w_vec <- setNames(region_weights$weight, region_weights$knock_region)
w_vec <- w_vec[c("NW","NEY","MID","EE","LON","SW","SE")]  # enforce column order
if (any(is.na(w_vec))) {
  warning("Missing weight for one or more regions -- check region_name_map against ",
          "unique(regression_data$itl1_name). Falling back to equal weights for ",
          "any missing region.")
  w_vec[is.na(w_vec)] <- 1/7
  w_vec <- w_vec / sum(w_vec)
}

table_s6_england <- table_s6_regional %>%
  rowwise() %>%
  mutate(england_mean = sum(c(NW, NEY, MID, EE, LON, SW, SE) * w_vec)) %>%
  ungroup() %>%
  dplyr::select(parameter, england_mean)

cat("\nEngland-representative absolute probabilities (real-admissions-weighted):\n")
print(table_s6_england)

p_ICU_max  <- table_s6_england$england_mean[table_s6_england$parameter == "p_ICU_max"]
p_HD_max   <- table_s6_england$england_mean[table_s6_england$parameter == "p_HD_max"]
p_ICUD_max <- table_s6_england$england_mean[table_s6_england$parameter == "p_ICUD_max"]
p_WD_max   <- table_s6_england$england_mean[table_s6_england$parameter == "p_WD_max"]

# ==============================================================
# SECTION 3 -- Combine Table S6 (absolute max) x Table S8 (relative
#   age-scaling) to get age-specific absolute probabilities
# ==============================================================
branch_probs <- table_s8_gf %>%
  mutate(
    p_ICU  = p_ICU_max  * scale_icu,
    p_HD   = p_HD_max   * scale_gen_death,
    p_ICUD = p_ICUD_max * scale_icu_death,
    p_WD   = p_WD_max   * scale_stepdown_death
  ) %>%
  dplyr::select(age_band_gf, p_ICU, p_HD, p_ICUD, p_WD)

cat("\nAge-specific absolute branch probabilities:\n")
print(branch_probs)

# ==============================================================
# SECTION 4 -- Durations from Table S2 (exact values, days)
# ==============================================================
dur_general_death     <- 10.3  # H_D
dur_general_recovery  <- 10.7  # H_R
dur_icu_direct_death  <- 11.8  # ICU_D
dur_icu_stepdown_death    <- 7.0 + 8.1   # ICU_WD + W_D = 15.1
dur_icu_stepdown_recovery <- 15.6 + 12.2 # ICU_WR + W_R = 27.8

# ==============================================================
# SECTION 5 -- Aggregate into age-specific gamma_hd, gamma_hr
#
#   Our model's HD/HR compartments are a 2-state simplification of
#   Knock's full hospital pathway (general ward + ICU + stepdown
#   combined). We aggregate Knock's pathway probabilities/durations
#   into a single "died in hospital" and "recovered in hospital"
#   rate for each age band, weighting each pathway's duration by its
#   probability of occurring.
# ==============================================================
gamma_table <- branch_probs %>%
  mutate(
    p_death_general = (1 - p_ICU) * p_HD,
    p_death_icu_direct = p_ICU * p_ICUD,
    p_death_icu_stepdown = p_ICU * (1 - p_ICUD) * p_WD,
    p_death_total = p_death_general + p_death_icu_direct + p_death_icu_stepdown,
    
    p_recovery_general = (1 - p_ICU) * (1 - p_HD),
    p_recovery_icu_stepdown = p_ICU * (1 - p_ICUD) * (1 - p_WD),
    p_recovery_total = p_recovery_general + p_recovery_icu_stepdown,
    
    mean_death_duration = (
      p_death_general    * dur_general_death +
        p_death_icu_direct * dur_icu_direct_death +
        p_death_icu_stepdown * dur_icu_stepdown_death
    ) / p_death_total,
    
    mean_recovery_duration = (
      p_recovery_general * dur_general_recovery +
        p_recovery_icu_stepdown * dur_icu_stepdown_recovery
    ) / p_recovery_total,
    
    gamma_hd = 1 / mean_death_duration,
    gamma_hr = 1 / mean_recovery_duration
  )

cat("\n--- Final age-specific gamma_hd, gamma_hr ---\n")
print(gamma_table %>% dplyr::select(age_band_gf, p_death_total, p_recovery_total,
                                    mean_death_duration, mean_recovery_duration,
                                    gamma_hd, gamma_hr))

cat("\nSanity check -- compare to previous flat placeholders:\n")
cat("  Old gamma_hd (flat, all ages): 0.0885  (mean duration 11.3 days)\n")
cat("  Old gamma_hr (flat, all ages): 0.0709  (mean duration 14.1 days)\n")
cat("  New gamma_hd range across ages:", round(range(gamma_table$gamma_hd), 4), "\n")
cat("  New gamma_hr range across ages:", round(range(gamma_table$gamma_hr), 4), "\n")


# ==============================================================
dir.create("data/parameters", recursive = TRUE, showWarnings = FALSE)
write_csv(
  gamma_table %>% dplyr::select(age_band_gf, p_death_total, p_recovery_total,
                                gamma_hd, gamma_hr),
  "data/parameters/gamma_hd_hr_by_age.csv"
)
cat("\nSaved to data/parameters/gamma_hd_hr_by_age.csv\n")
cat("\nNEXT STEP: update scripts 08/09 to read age-specific gamma_hd[a],\n")
cat("gamma_hr[a] from this file, instead of the flat placeholder values\n")
cat("gamma_hd=0.0885, gamma_hr=0.0709 currently hardcoded in both scripts.\n")
