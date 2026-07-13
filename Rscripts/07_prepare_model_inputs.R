# ==============================================================
# Script: 07_build_model_inputs.R
#
# Purpose: Build all static model inputs for the age x IMD stratified
#          transmission model:
#          (1) pi_matrix: clinical fraction by age x IMD decile
#          (2) contact matrices: one per IMD decile (urban population)
#          (3) h_mu_by_age: h_a and mu_ca_h derived from Knock et al. (2021)
#
# Inputs:
#   data/parameters/clin_frac.csv         (Goodfellow et al. 2024)
#   data/parameters/G.csv                 (POLYMOD intrinsic connectivity)
#   data/parameters/rural_age.csv         (ONS urban/rural age structure)
#   data/parameters/ihr_ifr_by_age_knock2021.csv
#
# Outputs:
#   data/parameters/pi_matrix.csv
#   data/parameters/contact_matrix_imd{1-10}.csv
#   data/parameters/h_mu_by_age.csv
# ==============================================================
library(dplyr)
library(tidyr)
library(readr)

age_levels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                "70 to 74","75+")

# 1. Clinical fraction matrix pi[age, IMD decile]
clin_frac <- read_csv("data/parameters/clin_frac.csv",
                      show_col_types = FALSE) %>%
  mutate(Ageband = factor(Ageband, levels = age_levels))

pi_matrix <- clin_frac %>%
  pivot_wider(names_from = IMD, values_from = clin_frac,
              names_prefix = "imd_") %>%
  arrange(Ageband) %>%
  select(Ageband, imd_1:imd_10)

cat("pi_matrix:", nrow(pi_matrix), "ages x", ncol(pi_matrix)-1, "deciles\n")
cat("pi range:", round(range(pi_matrix[,-1], na.rm=TRUE), 3), "\n")

write_csv(pi_matrix, "data/parameters/pi_matrix.csv")

# 2. Contact matrices: M[i,j] = G[i,j] * N[j] / sum(N)
#    N = urban age-specific population for each IMD decile
#    Urban population used throughout, consistent with Goodfellow et al. (2024)
G <- as.matrix(read_csv("data/parameters/G.csv", show_col_types = FALSE))
cat("G dimensions:", dim(G), "\n")

rural_age <- read_csv("data/parameters/rural_age.csv",
                      show_col_types = FALSE) %>%
  mutate(Age = factor(Age, levels = age_levels))

compute_contact_matrix <- function(imd_decile, urban = TRUE) {
  pop_vec <- rural_age %>%
    filter(IMD == imd_decile,
           rural == if (urban) "Urban" else "Rural") %>%
    arrange(Age) %>%
    pull(Population)
  M <- matrix(nrow = 17, ncol = 17)
  for (i in 1:17) for (j in 1:17) M[i,j] <- G[i,j] * pop_vec[j] / sum(pop_vec)
  M
}

dir.create("data/parameters", recursive = TRUE, showWarnings = FALSE)
for (d in 1:10) {
  M <- compute_contact_matrix(d, urban = TRUE)
  write_csv(as.data.frame(M),
            paste0("data/parameters/contact_matrix_imd", d, ".csv"),
            col_names = FALSE)
}
cat("Contact matrices saved for IMD deciles 1-10\n")

M1 <- compute_contact_matrix(1)
cat("Row sums (avg daily contacts): Under 1 =", round(sum(M1[1,]), 2),
    "| 75+ =", round(sum(M1[17,]), 2), "\n")

# 3. h_a and mu_ca_h by age (Knock et al. 2021, Table S9)
#    h_a     = IHR_a / pi_a   (P(hospitalised | clinical case))
#    mu_ca_h = IFR_a / IHR_a  (P(death | hospitalised))
#
#    Age band mapping from Knock to Goodfellow:
#    - Knock "0-4"        -> Goodfellow "Under 1" and "1 to 4" (same IHR/IFR)
#    - Knock "75-79","80+" -> Goodfellow "75+" (simple average)
#    - All other bands: direct 1-to-1 match
#
#    pi_a denominator: mean across IMD deciles (general population average)
ihr_ifr <- read_csv("data/parameters/ihr_ifr_by_age_knock2021.csv",
                    show_col_types = FALSE) %>%
  filter(age_band != "care_home")

ihr_0_4  <- ihr_ifr %>% filter(age_band == "0-4")
ihr_75up <- ihr_ifr %>% filter(age_band %in% c("75-79","80+")) %>%
  summarise(age_band = "75+", ihr = mean(ihr), ifr = mean(ifr))
ihr_mid  <- ihr_ifr %>% filter(!age_band %in% c("0-4","75-79","80+"))

pi_mean_by_age <- clin_frac %>%
  group_by(Ageband) %>%
  summarise(pi_mean = mean(clin_frac), .groups = "drop") %>%
  arrange(factor(Ageband, levels = age_levels))

h_mu <- data.frame(
  age_band = age_levels,
  ihr = c(rep(ihr_0_4$ihr, 2), ihr_mid$ihr, ihr_75up$ihr),
  ifr = c(rep(ihr_0_4$ifr, 2), ihr_mid$ifr, ihr_75up$ifr),
  pi_mean = pi_mean_by_age$pi_mean
) %>%
  mutate(
    h_a     = ihr / pi_mean,
    mu_ca_h = ifr / ihr
  )

cat("\nh_a and mu_ca_h by age:\n")
print(h_mu %>% select(age_band, ihr, pi_mean, h_a, mu_ca_h))

write_csv(h_mu, "data/parameters/h_mu_by_age.csv")
cat("\nAll model inputs saved to data/parameters/\n")
