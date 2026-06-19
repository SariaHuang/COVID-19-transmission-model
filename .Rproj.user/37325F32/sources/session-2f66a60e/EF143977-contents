# ==============================================================
# Script: 07_prepare_model_inputs.R
#
# Purpose: Process Goodfellow et al. (2024) data files and
#          IHR/IFR table (Knock et al. 2021) into the arrays
#          and matrices needed by the age-stratified odin model.
#
# Inputs:
#   data/parameters/clin_frac.csv       (Goodfellow repo /data/)
#   data/parameters/G.csv               (Goodfellow repo /data/)
#   data/parameters/rural_age.csv       (Goodfellow repo /data/)
#   data/parameters/ihr_ifr_by_age_knock2021.csv
#
# Outputs (all saved to data/parameters/):
#   pi_matrix.csv         -- clinical fraction, 17 ages x 10 IMD deciles
#   contact_matrix_imd1.csv ... contact_matrix_imd10.csv
#                          -- density-corrected 17x17 contact matrix
#                             for each IMD decile (urban, matching our
#                             LAD-level urban-weighted population)
#   h_mu_by_age.csv       -- h_a and mu_ca_h by age band
# ==============================================================

library(dplyr)
library(tidyr)
library(readr)

# ------------------------------------------------------------
# 1. Age band order 
# ------------------------------------------------------------
age_levels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                "70 to 74","75+")

# ------------------------------------------------------------
# 2. Clinical fraction matrix: pi[age, imd_decile]
# ------------------------------------------------------------
clin_frac <- read_csv("data/parameters/clin_frac.csv", show_col_types = FALSE) %>%
  mutate(Ageband = factor(Ageband, levels = age_levels))

pi_matrix <- clin_frac %>%
  pivot_wider(names_from = IMD, values_from = clin_frac,
              names_prefix = "imd_") %>%
  arrange(Ageband) %>%
  select(Ageband, imd_1:imd_10)

cat("pi_matrix dimensions:", nrow(pi_matrix), "ages x", ncol(pi_matrix)-1, "IMD deciles\n")
cat("pi range:", round(range(pi_matrix[,-1], na.rm=TRUE), 3), "\n")
cat("Sample: IMD decile 1 (most deprived), ages 60-64 to 75+:\n")
print(pi_matrix %>% filter(Ageband %in% c("60 to 64","65 to 69","70 to 74","75+")) %>%
        select(Ageband, imd_1, imd_10))

write_csv(pi_matrix, "data/parameters/pi_matrix.csv")

# ------------------------------------------------------------
# 3. Contact matrices: one per IMD decile
#    M[i,j] = G[i,j] * N[j] / sum(N)
#    where N = age-specific urban population for that IMD decile
#    (using urban population as the representative case, consistent
#    with Goodfellow's main analysis which focused on urban areas)
# ------------------------------------------------------------
G <- as.matrix(read_csv("data/parameters/G.csv", show_col_types = FALSE))
dim(G)  # should be 17 x 17

rural_age <- read_csv("data/parameters/rural_age.csv", show_col_types = FALSE) %>%
  mutate(Age = factor(Age, levels = age_levels))

compute_contact_matrix <- function(imd_decile, urban = TRUE) {
  setting <- if (urban) "Urban" else "Rural"
  pop_vec <- rural_age %>%
    filter(IMD == imd_decile, rural == setting) %>%
    arrange(Age) %>%
    pull(Population)
  M <- matrix(nrow = 17, ncol = 17)
  for (i in 1:17) {
    for (j in 1:17) {
      M[i, j] <- G[i, j] * pop_vec[j] / sum(pop_vec)
    }
  }
  return(M)
}

# Compute and save for all 10 IMD deciles (urban)
dir.create("data/parameters", recursive = TRUE, showWarnings = FALSE)
for (d in 1:10) {
  M <- compute_contact_matrix(d, urban = TRUE)
  fname <- paste0("data/parameters/contact_matrix_imd", d, ".csv")
  write.table(M, fname, sep = ",", row.names = FALSE, col.names = FALSE)
}
cat("Contact matrices saved for IMD deciles 1-10 (urban)\n")

# Quick sanity check: sum of row 1 of IMD decile 1 contact matrix
M1 <- compute_contact_matrix(1, urban = TRUE)
cat("Row 1 sum (Under 1, IMD decile 1) -- avg daily contacts:", round(sum(M1[1,]), 2), "\n")
cat("Row 17 sum (75+, IMD decile 1):", round(sum(M1[17,]), 2), "\n")

# ------------------------------------------------------------
# 4. h_a and mu_ca_h by age band (Knock 2021 Table S9)
#    h_a     = IHR_a / pi_a
#    mu_ca_h = IFR_a / IHR_a
#
#    NOTE: pi_a here is the general population (non-IMD-specific)
#    average clinical fraction, as used in Goodfellow's original
#    mu_ca formula. We use the IMD-decile-averaged pi (mean across
#    deciles) as the denominator. In the full stratified model,
#    h_a will vary by IMD decile because pi_a varies.
# ------------------------------------------------------------
ihr_ifr <- read_csv("data/parameters/ihr_ifr_by_age_knock2021.csv",
                    show_col_types = FALSE) %>%
  filter(age_band != "care_home")  # exclude care home row

# Match Knock's 17 age bands to Goodfellow's 17 age bands
# Knock uses 0-4, 5-9... 80+; Goodfellow uses Under 1, 1-4, 5-9... 75+
# The mapping is approximate for the youngest bands
knock_to_goodfellow <- c(
  "0-4"   = "Under 1",   # split across Under 1 and 1 to 4 in Goodfellow
  "0-4"   = "1 to 4",
  "5-9"   = "5 to 9",
  "10-14" = "10 to 14",
  "15-19" = "15 to 19",
  "20-24" = "20 to 24",
  "25-29" = "25 to 29",
  "30-34" = "30 to 34",
  "35-39" = "35 to 39",
  "40-44" = "40 to 44",
  "45-49" = "45 to 49",
  "50-54" = "50 to 54",
  "55-59" = "55 to 59",
  "60-64" = "60 to 64",
  "65-69" = "65 to 69",
  "70-74" = "70 to 74",
  "75-79" = "75+"        # 75+ in Goodfellow absorbs Knock's 75-79 and 80+
)

# Build a 17-row table matching Goodfellow's age bands
# Use 0-4 IHR/IFR for both "Under 1" and "1 to 4"
# Use average of 75-79 and 80+ for "75+"
ihr_0_4  <- ihr_ifr %>% filter(age_band == "0-4")
ihr_75up <- ihr_ifr %>% filter(age_band %in% c("75-79","80+")) %>%
  summarise(age_band = "75+", ihr = mean(ihr), ifr = mean(ifr))

ihr_middle <- ihr_ifr %>%
  filter(!age_band %in% c("0-4","75-79","80+")) %>%
  mutate(age_band_goodfellow = age_band)

# Mean pi across IMD deciles for each age band (for denominator)
pi_by_age_mean <- clin_frac %>%
  group_by(Ageband) %>%
  summarise(pi_mean = mean(clin_frac), .groups = "drop") %>%
  arrange(factor(Ageband, levels = age_levels))

# Assemble h and mu_ca_h table
h_mu_table <- data.frame(
  age_band_goodfellow = age_levels,
  ihr = c(rep(ihr_0_4$ihr, 2),
          ihr_middle$ihr,
          ihr_75up$ihr),
  ifr = c(rep(ihr_0_4$ifr, 2),
          ihr_middle$ifr,
          ihr_75up$ifr),
  pi_mean = pi_by_age_mean$pi_mean
) %>%
  mutate(
    h_a     = ihr / pi_mean,   # P(hosp | clinical case)
    mu_ca_h = ifr / ihr        # P(death | hospitalised)
  )

cat("\nh_a and mu_ca_h by age band:\n")
print(h_mu_table %>% select(age_band_goodfellow, ihr, pi_mean, h_a, mu_ca_h))

write_csv(h_mu_table, "data/parameters/h_mu_by_age.csv")
cat("\nAll model inputs saved to data/parameters/\n")
