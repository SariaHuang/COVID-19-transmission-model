# ==============================================================
# Script: 10_region_stratified_model.R
#
# Purpose: Add the region dimension to the age x IMD model (script 09).
#          Region enters only through population composition -- contact
#          matrices and proportion vectors are recomputed per region x IMD
#          combination using region-specific age structures and urban/rural
#          blending weights.
#
# Urban/rural blending:
#   M_blended = w_urban * M_urban + (1 - w_urban) * M_rural
#   w_urban is region x IMD specific, computed from LSOA-level Rural Urban
#   Classification (2011). Falls back to national IMD-level share for cells
#   with fewer than 5 LSOAs.
#
# Inputs:
#   data/processed/population_age_imd_region.csv       (script 05)
#   data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv (script 02)
#   data/raw/ruc_lsoa_2fold.csv
#   data/parameters/G.csv
#   data/parameters/rural_age.csv
#   data/parameters/urban_share_by_decile.csv
#   data/parameters/pi_matrix.csv, h_mu_by_age.csv     (script 07)
#   data/parameters/gamma_hd_hr_by_age.csv             (script 08)
#
# Output:
#   output/region_imd_hospital_gradient.png
# ==============================================================

library(odin2)
library(dust2)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

age_levels_gf <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                   "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                   "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                   "70 to 74","75+")

# Load model if not already compiled (from script 09)
if (!exists("age_seird_hosp")) {
  SKIP09_RUN <- TRUE
  source("Rscripts/09_age_imd_stratified_odin.R")
  rm(SKIP09_RUN)
}

# Under-1 / 1-to-4 split ratio from Goodfellow's rural_age data
rural_age_gf <- read_csv("data/parameters/rural_age.csv",
                         show_col_types = FALSE)

split_ratio_under1 <- rural_age_gf %>%
  filter(Age %in% c("Under 1","1 to 4")) %>%
  group_by(Age) %>%
  summarise(pop = sum(Population), .groups = "drop") %>%
  { .$pop[.$Age == "Under 1"] / sum(.$pop) }

cat("Under-1 share of (Under 1 + 1 to 4):", round(split_ratio_under1, 3), "\n")

# Re-bin region population to Goodfellow's 17 age bands
pop_region <- read_csv("data/processed/population_age_imd_region.csv",
                       show_col_types = FALSE) %>%
  mutate(age_band = as.character(age_band))

pop_0_4 <- pop_region %>% filter(age_band == "0-4")
pop_rebinned <- bind_rows(
  pop_0_4 %>% mutate(age_band_gf = "Under 1",
                     population  = population * split_ratio_under1),
  pop_0_4 %>% mutate(age_band_gf = "1 to 4",
                     population  = population * (1 - split_ratio_under1)),
  pop_region %>%
    filter(age_band != "0-4") %>%
    mutate(age_band_gf = case_when(
      age_band == "5-9"              ~ "5 to 9",
      age_band == "10-14"            ~ "10 to 14",
      age_band == "15-19"            ~ "15 to 19",
      age_band == "20-24"            ~ "20 to 24",
      age_band == "25-29"            ~ "25 to 29",
      age_band == "30-34"            ~ "30 to 34",
      age_band == "35-39"            ~ "35 to 39",
      age_band == "40-44"            ~ "40 to 44",
      age_band == "45-49"            ~ "45 to 49",
      age_band == "50-54"            ~ "50 to 54",
      age_band == "55-59"            ~ "55 to 59",
      age_band == "60-64"            ~ "60 to 64",
      age_band == "65-69"            ~ "65 to 69",
      age_band == "70-74"            ~ "70 to 74",
      age_band %in% c("75-79","80+") ~ "75+",
      TRUE                            ~ NA_character_
    ))
) %>%
  group_by(itl1_name, lad_imd_decile, age_band_gf) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
  mutate(age_band_gf = factor(age_band_gf, levels = age_levels_gf))

cat("Rebinned population rows:", nrow(pop_rebinned),
    "| Expected:", 9*10*17, "\n")

# Region x IMD urban share from LSOA-level RUC (2011)
ruc_lsoa   <- read_csv("data/raw/ruc_lsoa_2fold.csv", show_col_types = FALSE)
lsoa_lookup <- read_csv(
  "data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv",
  show_col_types = FALSE
)

lsoa_ruc <- lsoa_lookup %>%
  select(lsoa_code, itl1_name, imd_decile) %>%
  inner_join(ruc_lsoa, by = "lsoa_code")

urban_share_region_imd <- lsoa_ruc %>%
  group_by(itl1_name, imd_decile) %>%   
  summarise(
    n_urban = sum(urban_rural == "Urban"),
    n_total = n(),
    w_urban = n_urban / n_total,
    .groups = "drop"
  )

# National fallback urban share per IMD decile
urban_share_national <- read_csv("data/parameters/urban_share_by_decile.csv",
                                 show_col_types = FALSE)

# Contact matrix from population vector
G <- as.matrix(read_csv("data/parameters/G.csv", show_col_types = FALSE))

rural_age_full <- read_csv("data/parameters/rural_age.csv",
                           show_col_types = FALSE) %>%
  mutate(Age = factor(Age, levels = age_levels_gf))

make_contact_matrix <- function(pop_vec) {
  M <- matrix(nrow = 17, ncol = 17)
  for (i in 1:17) for (j in 1:17) M[i,j] <- G[i,j] * pop_vec[j] / sum(pop_vec)
  M
}

get_region_decile_inputs <- function(region, decile) {
  pop_vec <- pop_rebinned %>%
    filter(itl1_name == region, lad_imd_decile == decile) %>%
    arrange(age_band_gf) %>%
    pull(population)
  
  if (length(pop_vec) != 17 || sum(pop_vec, na.rm = TRUE) == 0) return(NULL)
  
  # Urban weight: region x IMD specific, with national fallback
  w_row <- urban_share_region_imd %>%
    filter(itl1_name == region, imd_decile == decile)  
  
  w_urban <- if (nrow(w_row) == 0 || w_row$n_total < 5) {
    urban_share_national$w_urban[urban_share_national$IMD == decile]
  } else {
    w_row$w_urban
  }
  
  pop_urban <- rural_age_full %>%
    filter(IMD == decile, rural == "Urban") %>%
    arrange(Age) %>% pull(Population)
  pop_rural <- rural_age_full %>%
    filter(IMD == decile, rural == "Rural") %>%
    arrange(Age) %>% pull(Population)
  
  M_blended <- w_urban * make_contact_matrix(pop_urban) +
    (1 - w_urban) * make_contact_matrix(pop_rural)
  
  list(contact    = M_blended,
       proportion = pop_vec / sum(pop_vec))
}

# Model parameters
pi_matrix    <- read_csv("data/parameters/pi_matrix.csv", show_col_types = FALSE)
h_mu         <- read_csv("data/parameters/h_mu_by_age.csv", show_col_types = FALSE)
gamma_hd_hr  <- read_csv("data/parameters/gamma_hd_hr_by_age.csv", show_col_types = FALSE)
gamma_hd_vec <- gamma_hd_hr$gamma_hd
gamma_hr_vec <- gamma_hd_hr$gamma_hr

# State index for Adm[1:17] in dust2 output (170 states total)
adm_idx <- 154:170

run_region_decile <- function(region, decile) {
  inputs <- get_region_decile_inputs(region, decile)
  if (is.null(inputs)) return(NULL)
  
  sys <- dust2::dust_system_create(age_seird_hosp, list(
    S0         = inputs$proportion - c(rep(0,7), 1e-3, rep(0,9)),
    Ip0        = c(rep(0,7), 1e-3, rep(0,9)),
    proportion = inputs$proportion,
    pi_a       = pi_matrix[[paste0("imd_", decile)]],
    h_a        = h_mu$h_a,
    mu_ca_h    = h_mu$mu_ca_h,
    contact    = inputs$contact,
    gam_hd     = gamma_hd_vec,
    gam_hr     = gamma_hr_vec
  ))
  dust2::dust_system_set_state_initial(sys)
  out <- dust2::dust_system_simulate(sys, seq(0, 365, by = 1))
  
  data.frame(
    day            = 0:365,
    region         = region,
    imd_decile     = decile,
    cum_admissions = colSums(out[adm_idx, ])
  ) %>%
    mutate(new_adm_per1000 = c(0, diff(cum_admissions)) * 1000)
}

regions <- unique(pop_rebinned$itl1_name)
cat("\nRunning model for", length(regions), "regions x 10 IMD deciles...\n")

all_combos <- expand.grid(region = regions, decile = 1:10,
                          stringsAsFactors = FALSE)

results_region <- bind_rows(Map(function(r, d) {
  cat("  ", r, "decile", d, "\n")
  run_region_decile(r, d)
}, all_combos$region, all_combos$decile))

cat("Completed:", n_distinct(paste(results_region$region,
                                   results_region$imd_decile)),
    "of", nrow(all_combos), "combinations\n")

# Summary
summary_region <- results_region %>%
  group_by(region, imd_decile) %>%
  summarise(peak_adm_per1000 = max(new_adm_per1000, na.rm = TRUE),
            .groups = "drop")

cat("\nPeak admissions by region and IMD decile:\n")
print(summary_region %>% arrange(region, imd_decile), n = 50)

# Plot
decile_colours <- c(
  "#67001f","#b2182b","#d6604d","#f4a582","#fddbc7",
  "#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"
)

p <- ggplot(results_region,
            aes(x = day, y = new_adm_per1000,
                colour = factor(imd_decile, levels = 10:1),
                group  = imd_decile)) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  facet_wrap(~ region, scales = "free_y") +
  scale_colour_manual(values = rev(decile_colours),
                      name = "IMD decile\n(1 = most\ndeprived)") +
  labs(
    title    = "Modelled hospital admissions by region and IMD decile",
    subtitle = "Region enters via population composition only (age x IMD structure per region)",
    x        = "Day of epidemic",
    y        = "New hospital admissions per 1,000 population",
    caption  = "Parameters: Goodfellow et al. (2024) + Knock et al. (2021)."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "#555555"),
    plot.caption     = element_text(size = 7, colour = "#888888", hjust = 0),
    strip.text       = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

print(p)
dir.create("output", showWarnings = FALSE)
ggsave("output/region_imd_hospital_gradient.png", p,
       width = 13, height = 9, dpi = 200)
cat("Plot saved: output/region_imd_hospital_gradient.png\n")
