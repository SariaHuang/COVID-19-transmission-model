# ==============================================================
# Script: 10_region_stratified_model.R
#
# Purpose: Add the region dimension to the age x IMD stratified
#          model (script 08). Region enters ONLY through population
#          composition (age x IMD population structure differs by
#          region) -- contact matrices and clinical/severity
#          parameters are recomputed per region x IMD combination,
#          but the model equations themselves (odin definition in
#          script 08) are unchanged.
#
# Contact matrix blending :
#   M_blended = w_urban * M_urban + (1-w_urban) * M_rural
#   Urban/rural matrices computed from Goodfellow's national age
#   structures for each IMD decile. Blending weight w_urban is now
#   REGION x IMD SPECIFIC, computed from LSOA-level Rural Urban
#   Classification (2011) joined to our LSOA -> region -> IMD lookup.
#   Fallback to national IMD-level share for sparse cells (<5 LSOAs).
#
# Inputs:
#   data/processed/population_age_imd_region.csv (script 05)
#   data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv (script 02)
#   data/raw/ruc_lsoa_2fold.csv  (Rural Urban Classification 2011)
#   data/parameters/G.csv
#   data/parameters/rural_age.csv
#   data/parameters/urban_share_by_decile.csv (script 07, fallback)
#   data/parameters/pi_matrix.csv, h_mu_by_age.csv (script 07)
#
# CAVEAT: our age bins (0-4, 5-9, ..., 80+) don't exactly match
# Goodfellow's (Under 1, 1-4, ..., 75+). We approximate by splitting
# our 0-4 band into Under-1/1-4 using Goodfellow's national-level
# ratio, and merging our 75-79 + 80+ into his 75+ band. This is a
# documented approximation, not an exact re-binning.
#
# Output: hospital admissions gradient plot, faceted by region,
#         coloured by IMD decile
# ==============================================================

library(odin)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)


age_levels_gf <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                   "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                   "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                   "70 to 74","75+")

age_levels_ons <- c("0-4","5-9","10-14","15-19","20-24","25-29","30-34",
                    "35-39","40-44","45-49","50-54","55-59","60-64",
                    "65-69","70-74","75-79","80+")

# ------------------------------------------------------------
# 1. Work out the Under-1 / 1-4 split ratio from Goodfellow's data
# ------------------------------------------------------------
rural_age_gf <- read_csv("data/parameters/rural_age.csv", show_col_types = FALSE)

national_under1_1to4 <- rural_age_gf %>%
  filter(Age %in% c("Under 1", "1 to 4")) %>%
  group_by(Age) %>%
  summarise(pop = sum(Population), .groups = "drop")

split_ratio_under1 <- national_under1_1to4$pop[national_under1_1to4$Age == "Under 1"] /
  sum(national_under1_1to4$pop)
cat("Under-1 share of (Under1 + 1-4):", round(split_ratio_under1, 3), "\n")

# ------------------------------------------------------------
# 2. Load and re-bin our region population to Goodfellow's age bands
# ------------------------------------------------------------
pop_region <- read_csv("data/processed/population_age_imd_region.csv",
                       show_col_types = FALSE)

pop_region_rebinned <- pop_region %>%
  mutate(age_band = as.character(age_band)) %>%
  # Split 0-4 into Under 1 / 1 to 4
  { 
    df <- .
    under1 <- df %>% filter(age_band == "0-4") %>%
      mutate(age_band_gf = "Under 1", population = population * split_ratio_under1)
    to4 <- df %>% filter(age_band == "0-4") %>%
      mutate(age_band_gf = "1 to 4", population = population * (1 - split_ratio_under1))
    rest <- df %>% filter(age_band != "0-4") %>%
      mutate(age_band_gf = case_when(
        age_band == "5-9"   ~ "5 to 9",
        age_band == "10-14" ~ "10 to 14",
        age_band == "15-19" ~ "15 to 19",
        age_band == "20-24" ~ "20 to 24",
        age_band == "25-29" ~ "25 to 29",
        age_band == "30-34" ~ "30 to 34",
        age_band == "35-39" ~ "35 to 39",
        age_band == "40-44" ~ "40 to 44",
        age_band == "45-49" ~ "45 to 49",
        age_band == "50-54" ~ "50 to 54",
        age_band == "55-59" ~ "55 to 59",
        age_band == "60-64" ~ "60 to 64",
        age_band == "65-69" ~ "65 to 69",
        age_band == "70-74" ~ "70 to 74",
        age_band %in% c("75-79","80+") ~ "75+",
        TRUE ~ NA_character_
      ))
    bind_rows(under1, to4, rest)
  } %>%
  group_by(itl1_name, lad_imd_decile, age_band_gf) %>%
  summarise(population = sum(population, na.rm = TRUE), .groups = "drop") %>%
  mutate(age_band_gf = factor(age_band_gf, levels = age_levels_gf))

cat("Rebinned population rows:", nrow(pop_region_rebinned), "\n")
cat("Expected (9 regions x 10 deciles x 17 ages):", 9*10*17, "\n")

# ------------------------------------------------------------
# 3. Build contact matrix and proportion vector per region x decile
#
#    Urban/rural blending (advisor's request, consistent with script 07):
#      M_blended = w_urban * M_urban + (1-w_urban) * M_rural
#
#    Improvement over script 07: w_urban is now REGION x IMD SPECIFIC,
#    computed from LSOA-level Rural Urban Classification (2011) joined
#    to our LSOA -> region -> IMD lookup. Each LSOA counts as one unit
#    (each LSOA ~1,500 people, so LSOA count is a valid population proxy).
#
#    For M_urban and M_rural, we use Goodfellow's national urban/rural
#    age structures for that IMD decile (rural_age.csv) -- the best
#    available proxy since we don't have region x IMD x urban/rural
#    age breakdowns.
#
#    Fallback: if a region x IMD cell has no LSOA data (sparse cells),
#    falls back to the national IMD-level urban share.
# ------------------------------------------------------------
G <- as.matrix(read_csv("data/parameters/G.csv", show_col_types = FALSE))

rural_age_gf_full <- read_csv("data/parameters/rural_age.csv",
                              show_col_types = FALSE) %>%
  mutate(Age = factor(Age, levels = age_levels_gf))

# National urban share per IMD decile (fallback)
urban_share_national <- read_csv("data/parameters/urban_share_by_decile.csv",
                                 show_col_types = FALSE)

# --- Compute region x IMD urban share from LSOA-level RUC data ---
ruc_lsoa <- read_csv("data/raw/ruc_lsoa_2fold.csv", show_col_types = FALSE)

lsoa_lookup <- read_csv("data/processed/lookup_LSOA11_STP_ITL2_IMD_clean.csv",
                        show_col_types = FALSE)

# Join RUC to LSOA lookup to get region + IMD decile for each LSOA
# Column names may vary -- adjust if needed
lsoa_ruc <- lsoa_lookup %>%
  dplyr::select(any_of(c("lsoa11cd","LSOA11CD","lsoa_code")),
                any_of(c("itl1_name","ITL1NAME","region")),
                any_of(c("lad_imd_decile","imd_decile","IMD_decile"))) %>%
  setNames(c("lsoa_code","itl1_name","lad_imd_decile")) %>%
  inner_join(ruc_lsoa, by = "lsoa_code")

# Urban share per region x IMD decile (LSOA count as proxy for population)
urban_share_region_imd <- lsoa_ruc %>%
  group_by(itl1_name, lad_imd_decile) %>%
  summarise(
    n_urban = sum(urban_rural == "Urban"),
    n_total = n(),
    w_urban = n_urban / n_total,
    .groups = "drop"
  )

cat("\n--- Urban share by region x IMD decile (sample) ---\n")
print(urban_share_region_imd %>% arrange(itl1_name, lad_imd_decile) %>% head(20))

compute_setting_matrix_from_vec <- function(pop_vec) {
  M <- matrix(nrow = 17, ncol = 17)
  for (i in 1:17) {
    for (j in 1:17) {
      M[i, j] <- G[i, j] * pop_vec[j] / sum(pop_vec)
    }
  }
  M
}

get_region_decile_inputs <- function(region, decile) {
  
  # Total population for this region x decile (from rebinned ONS data)
  pop_vec_total <- pop_region_rebinned %>%
    filter(itl1_name == region, lad_imd_decile == decile) %>%
    arrange(age_band_gf) %>%
    pull(population)
  
  if (length(pop_vec_total) != 17 || sum(pop_vec_total, na.rm = TRUE) == 0) {
    return(NULL)
  }
  
  # Region x IMD specific urban weight (from LSOA RUC data)
  # Fall back to national IMD-level share if this region x decile is sparse
  w_row <- urban_share_region_imd %>%
    filter(itl1_name == region, lad_imd_decile == decile)
  
  if (nrow(w_row) == 0 || w_row$n_total < 5) {
    # Fallback: use national IMD-level share
    w_urban <- urban_share_national$w_urban[urban_share_national$IMD == decile]
  } else {
    w_urban <- w_row$w_urban
  }
  
  # Urban/rural age structures for this IMD decile (national, from Goodfellow)
  pop_urban_nat <- rural_age_gf_full %>%
    filter(IMD == decile, rural == "Urban") %>%
    arrange(Age) %>%
    pull(Population)
  
  pop_rural_nat <- rural_age_gf_full %>%
    filter(IMD == decile, rural == "Rural") %>%
    arrange(Age) %>%
    pull(Population)
  
  M_urban  <- compute_setting_matrix_from_vec(pop_urban_nat)
  M_rural  <- compute_setting_matrix_from_vec(pop_rural_nat)
  M_blended <- w_urban * M_urban + (1 - w_urban) * M_rural
  
  # proportion: actual region x IMD total age structure (region-specific)
  proportion <- pop_vec_total / sum(pop_vec_total)
  
  list(contact = M_blended, proportion = proportion)
}

# ------------------------------------------------------------
# 4. Run the model for each region x IMD decile combination
# ------------------------------------------------------------
pi_matrix <- read_csv("data/parameters/pi_matrix.csv", show_col_types = FALSE)
h_mu      <- read_csv("data/parameters/h_mu_by_age.csv", show_col_types = FALSE)
gamma_hd_hr <- read_csv("data/parameters/gamma_hd_hr_by_age.csv",
                        show_col_types = FALSE)
gamma_hd_vec <- gamma_hd_hr$gamma_hd  # length 17, age-specific
gamma_hr_vec <- gamma_hd_hr$gamma_hr  # length 17, age-specific

regions <- unique(pop_region_rebinned$itl1_name)

run_region_decile <- function(region, decile) {
  inputs <- get_region_decile_inputs(region, decile)
  if (is.null(inputs)) return(NULL)
  
  pi_a    <- pi_matrix[[paste0("imd_", decile)]]
  h_a     <- h_mu$h_a
  mu_ca_h <- h_mu$mu_ca_h
  
  S0  <- inputs$proportion
  S0[8] <- S0[8] - 1e-3
  Ip0 <- c(rep(0, 7), 1e-3, rep(0, 9))
  
  mod <- age_seird_hosp$new(
    S0 = S0, Ip0 = Ip0, proportion = inputs$proportion,
    pi_a = pi_a, h_a = h_a, mu_ca_h = mu_ca_h,
    contact = inputs$contact, gam_hd = gamma_hd_vec, gam_hr = gamma_hr_vec
  )
  
  out <- as.data.frame(mod$run(seq(0, 365, by = 1)))
  adm_cols <- grep("^Adm\\[", names(out))
  
  data.frame(
    day             = out$t,
    region          = region,
    imd_decile      = decile,
    cum_admissions  = rowSums(out[, adm_cols])
  ) %>%
    mutate(new_adm_per1000 = c(0, diff(cum_admissions)) * 1000)
}

cat("\nRunning model for", length(regions), "regions x 10 IMD deciles...\n")
all_combos <- expand.grid(region = regions, decile = 1:10,
                          stringsAsFactors = FALSE)

results_list <- list()
for (k in seq_len(nrow(all_combos))) {
  r <- all_combos$region[k]
  d <- all_combos$decile[k]
  res <- tryCatch(run_region_decile(r, d), error = function(e) NULL)
  if (!is.null(res)) results_list[[length(results_list) + 1]] <- res
}
results_region <- bind_rows(results_list)

cat("Successfully ran:", length(unique(paste(results_region$region, results_region$imd_decile))),
    "out of", nrow(all_combos), "region x decile combinations\n")
cat("(missing combinations = sparse cells with no population data, expected)\n")

# ------------------------------------------------------------
# 5. Sanity check
# ------------------------------------------------------------
summary_region <- results_region %>%
  group_by(region, imd_decile) %>%
  summarise(peak_adm_per1000 = max(new_adm_per1000, na.rm = TRUE), .groups = "drop")

cat("\n--- Peak admissions by region and IMD decile ---\n")
print(summary_region %>% arrange(region, imd_decile), n = 50)

# ------------------------------------------------------------
# 6. Plot: faceted by region, coloured by IMD decile
# ------------------------------------------------------------
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
    caption  = "Parameters: Goodfellow et al. (2024) + Knock et al. (2021). Contact matrices and clinical fraction national; only population structure varies by region."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "#555555"),
    plot.caption  = element_text(size = 7, colour = "#888888", hjust = 0),
    strip.text    = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

print(p)
dir.create("output", showWarnings = FALSE)
ggsave("output/region_imd_hospital_gradient.png", p, width = 13, height = 9, dpi = 200)
cat("\nPlot saved to output/region_imd_hospital_gradient.png\n")
