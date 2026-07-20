# ==============================================================
# Script: 16_counterfactual.R
#
# Purpose: Counterfactual analysis -- quantify how much of the
#          COVID-19 hospitalisation and death burden is attributable
#          to health inequality in clinical severity (pi_a).
#
# Two counterfactual scenarios:
#   CF1: all deciles assigned decile 10 (least deprived) pi_a
#        -- eliminates health inequality entirely
#   CF2: all deciles assigned decile 5 (median) pi_a
#        -- partial equalisation; more realistic policy target
#
# h_a is NOT replaced in either scenario: h_a = IHR_a / pi_a,
#   and IHR_a is deprivation-invariant (Knock et al. 2021).
#   Replacing h_a would mechanically inflate it for deprived
#   deciles, producing artefactual admission increases.
#
# Beta: fixed at decile 1 posterior median throughout, consistent
#   with sensitivity analysis in script 15.
#
# Population: blended urban/rural (consistent with scripts 09/13).
#   adm_avoided_abs uses blended pop_size per decile.
#
# Outputs (output/plots/counterfactual/):
#   16_fig1_attributable_fraction.png
#   16_fig2_burden_comparison_region.png
#   16_fig3_national_summary.png
#   16_counterfactual_full.csv
#   16_counterfactual_decile_summary.csv
# ==============================================================

library(odin2)
library(dust2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(purrr)
library(stringr)
library(scales)
library(forcats)

if (!exists("age_seird_hosp")) {
  SKIP09_RUN <- TRUE
  source("Rscripts/08_refine_gamma_hd_hr.R")
  source("Rscripts/09_age_imd_stratified_odin.R")
  rm(SKIP09_RUN)
  cat("odin2 model loaded.\n")
}
stopifnot(exists("get_blended_inputs"))

if (!exists("age_labels")) {
  age_labels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                  "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                  "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                  "70 to 74","75+")
}

if (!exists("beta_posterior")) {
  burnin <- 500
  beta_posterior <- sapply(1:10, function(d) {
    fit <- readRDS(paste0("output/fitting/fitted_samples_imd", d, ".rds"))
    median(sapply(1:3, function(ch)
      median(exp(fit$pars[1, (burnin+1):2000, ch]))))
  })
}

beta_fixed <- beta_posterior[1]
cat(sprintf("Fixed beta: %.5f (decile 1 posterior median)\n", beta_fixed))

dir.create("output/plots/counterfactual", recursive = TRUE,
           showWarnings = FALSE)

pi_a_d10 <- pi_matrix[["imd_10"]]
pi_a_d5  <- pi_matrix[["imd_5"]]

cat("CF1 pi_a (decile 10 mean):", round(mean(pi_a_d10), 3), "\n")
cat("CF2 pi_a (decile 5 mean): ", round(mean(pi_a_d5),  3), "\n")
cat("Baseline pi_a (decile 1 mean):",
    round(mean(pi_matrix[["imd_1"]]), 3), "\n\n")

adm_idx <- 154:170
d_idx   <- 120:136

# Run odin2 model with optional pi_a override using blended population
run_epidemic_fit_age_cf <- function(imd_decile, beta,
                                    pi_a_override = NULL) {
  contact <- as.matrix(read.csv(
    paste0("data/parameters/contact_matrix_imd", imd_decile, ".csv"),
    header = FALSE
  ))
  
  pi_a    <- if (!is.null(pi_a_override)) pi_a_override else
    pi_matrix[[paste0("imd_", imd_decile)]]
  h_a     <- h_mu$h_a
  mu_ca_h <- h_mu$mu_ca_h
  
  blended    <- get_blended_inputs(imd_decile)
  proportion <- blended$proportion
  pop_by_age <- proportion * blended$pop_size
  
  S0    <- proportion
  S0[8] <- S0[8] - 1e-4
  Ip0   <- c(rep(0, 7), 1e-4, rep(0, 9))
  
  sys <- dust2::dust_system_create(age_seird_hosp, list(
    S0         = S0,
    Ip0        = Ip0,
    proportion = proportion,
    pi_a       = pi_a,
    h_a        = h_a,
    mu_ca_h    = mu_ca_h,
    contact    = contact,
    gam_hd     = gamma_hd_vec,
    gam_hr     = gamma_hr_vec,
    susc       = beta
  ))
  dust2::dust_system_set_state_initial(sys)
  out <- dust2::dust_system_simulate(sys, seq(0, 365, by = 1))
  
  map_dfr(seq_along(age_labels), function(a) {
    data.frame(
      day               = 0:365,
      imd_decile        = imd_decile,
      age_group         = age_labels[a],
      age_idx           = a,
      pop_blended       = pop_by_age[a],
      cum_adm_per1000   = out[adm_idx[a], ] * 1000,
      cum_death_per1000 = out[d_idx[a],   ] * 1000
    )
  })
}

run_all_deciles <- function(pi_a_override = NULL, label = "") {
  cat("Running:", label, "\n")
  map_dfr(1:10, function(d) {
    cat("  Decile", d, "\n")
    run_epidemic_fit_age_cf(d, beta_fixed, pi_a_override = pi_a_override)
  }) %>%
    group_by(imd_decile, age_group, age_idx, pop_blended) %>%
    slice_max(day, n = 1) %>%
    ungroup()
}

base_fd <- run_all_deciles(pi_a_override = NULL,     label = "Baseline")
cf1_fd  <- run_all_deciles(pi_a_override = pi_a_d10, label = "CF1 (decile 10 pi_a)")
cf2_fd  <- run_all_deciles(pi_a_override = pi_a_d5,  label = "CF2 (decile 5 pi_a)")

make_results <- function(base, cf, cf_label) {
  base %>%
    rename(adm_base = cum_adm_per1000, death_base = cum_death_per1000) %>%
    left_join(
      cf %>% select(imd_decile, age_group, age_idx,
                    adm_cf = cum_adm_per1000, death_cf = cum_death_per1000),
      by = c("imd_decile","age_group","age_idx")
    ) %>%
    mutate(
      scenario          = cf_label,
      adm_diff          = adm_base - adm_cf,
      death_diff        = death_base - death_cf,
      adm_attr_frac     = pmax(adm_diff   / adm_base   * 100, 0),
      death_attr_frac   = pmax(death_diff / death_base * 100, 0),
      adm_avoided_abs   = adm_diff   / 1000 * pop_blended,
      death_avoided_abs = death_diff / 1000 * pop_blended
    )
}

results_cf1 <- make_results(base_fd, cf1_fd, "CF1: decile 10 \u03c0\u2090")
results_cf2 <- make_results(base_fd, cf2_fd, "CF2: decile 5 \u03c0\u2090")
results_all <- bind_rows(results_cf1, results_cf2)
results     <- results_cf1

decile_summary <- results_all %>%
  group_by(scenario, imd_decile) %>%
  summarise(
    adm_base_total       = sum(adm_base),
    adm_cf_total         = sum(adm_cf),
    adm_attr_frac_mean   = mean(adm_attr_frac),
    death_attr_frac_mean = mean(death_attr_frac),
    adm_avoided_abs      = sum(adm_avoided_abs),
    death_avoided_abs    = sum(death_avoided_abs),
    .groups = "drop"
  )

cat("\n--- National totals by scenario ---\n")
results_all %>%
  group_by(scenario) %>%
  summarise(
    avoided_adm        = round(sum(adm_avoided_abs)),
    avoided_death      = round(sum(death_avoided_abs)),
    mean_attr_frac_adm = round(mean(adm_attr_frac), 1),
    .groups = "drop"
  ) %>%
  print()

theme_pub <- theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12,
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 9, colour = "#444444",
                                    margin = margin(b = 8)),
    plot.caption     = element_text(size = 7.5, colour = "#888888",
                                    hjust = 0, margin = margin(t = 8)),
    axis.title       = element_text(size = 9.5),
    axis.text        = element_text(size = 8.5, colour = "#333333"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#eeeeee"),
    legend.title     = element_text(size = 8.5, face = "bold"),
    legend.text      = element_text(size = 8),
    plot.margin      = margin(12, 12, 8, 12)
  )

cf_caption <- paste0(
  "CF1: \u03c0\u2090 replaced with decile 10 (least deprived) values. ",
  "CF2: \u03c0\u2090 replaced with decile 5 (median) values. ",
  "h\u2090 unchanged (age-specific IHR, Knock et al. 2021). ",
  "\u03b2 fixed at decile 1 posterior median (0.031). ",
  "Population: blended urban/rural per decile. ",
  "Model: age \u00d7 IMD SEIRD + hospital (odin2). ",
  "Parameters: Goodfellow et al. (2024) + Knock et al. (2021)."
)

# Figure 1: Attributable fraction
cat("\nPlot 1: Attributable fraction (CF1 and CF2)...\n")

fig1_data <- decile_summary %>%
  select(scenario, imd_decile, adm_attr_frac_mean) %>%
  mutate(scenario = factor(scenario,
                           levels = c("CF1: decile 10 \u03c0\u2090",
                                      "CF2: decile 5 \u03c0\u2090")))

p1 <- ggplot(fig1_data,
             aes(x = imd_decile, y = adm_attr_frac_mean, fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
  geom_text(aes(label = paste0(round(adm_attr_frac_mean, 1), "%")),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.8, colour = "#333333") +
  scale_x_continuous(breaks = 1:10,
                     labels = c("1\n(most\ndeprived)", 2:9,
                                "10\n(least\ndeprived)")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     labels = function(x) paste0(x, "%")) +
  scale_fill_manual(
    values = c("CF1: decile 10 \u03c0\u2090" = "#2166ac",
               "CF2: decile 5 \u03c0\u2090"  = "#74add1"),
    name = "Counterfactual"
  ) +
  labs(
    title    = "Proportion of COVID-19 admissions attributable to health inequality",
    subtitle = paste0("CF1: all deciles assigned decile 10 \u03c0\u2090 | ",
                      "CF2: all deciles assigned decile 5 \u03c0\u2090"),
    x        = "IMD deprivation decile",
    y        = "Attributable fraction (%)",
    caption  = cf_caption
  ) +
  theme_pub + theme(legend.position = "top")

ggsave("output/plots/counterfactual/16_fig1_attributable_fraction.png",
       p1, width = 11, height = 6.5, dpi = 200)
cat("  Saved: 16_fig1_attributable_fraction.png\n")
print(p1)

# Figure 2: Region comparison
cat("Plot 2: Region comparison...\n")

region_pop <- read_csv("data/processed/population_age_imd_region.csv",
                       show_col_types = FALSE)

region_pop_agg <- region_pop %>%
  mutate(age_band_model = case_when(
    age_band %in% c("75-79","80+") ~ "75+", TRUE ~ age_band)) %>%
  group_by(itl1_name, lad_imd_decile, age_band_model) %>%
  summarise(population = sum(population), .groups = "drop")

region_total <- region_pop_agg %>%
  group_by(itl1_name) %>%
  summarise(total_pop = sum(population), .groups = "drop")

age_map <- tibble(
  age_group      = age_labels,
  age_band_model = c("0-4","0-4","5-9","10-14","15-19","20-24",
                     "25-29","30-34","35-39","40-44","45-49","50-54",
                     "55-59","60-64","65-69","70-74","75+"))

pop_w <- rural_age %>%
  filter(rural == "Urban", Age %in% c("Under 1","1 to 4")) %>%
  group_by(Age) %>%
  summarise(pop = sum(Population), .groups = "drop") %>%
  mutate(w = pop / sum(pop))
w_u1 <- pop_w$w[pop_w$Age == "Under 1"]
w_14 <- pop_w$w[pop_w$Age == "1 to 4"]

agg_to_bands <- function(df, adm_col) {
  df %>%
    left_join(age_map, by = "age_group") %>%
    group_by(imd_decile, age_band_model) %>%
    summarise(
      adm_per1000 = if (n() == 2)
        weighted.mean(.data[[adm_col]],
                      c(w_u1,w_14)[match(age_group,c("Under 1","1 to 4"))])
      else first(.data[[adm_col]]),
      .groups = "drop")
}

compute_region <- function(rates, suffix) {
  region_pop_agg %>%
    left_join(rates, by = c("lad_imd_decile"="imd_decile","age_band_model")) %>%
    mutate(abs_adm = (adm_per1000/1000)*population) %>%
    group_by(itl1_name, lad_imd_decile) %>%
    summarise(abs_adm = sum(abs_adm, na.rm=TRUE), .groups="drop") %>%
    left_join(region_total, by="itl1_name") %>%
    mutate(!!paste0("adm_",suffix) := abs_adm/total_pop*1000) %>%
    select(itl1_name, lad_imd_decile, !!paste0("adm_",suffix))
}

region_compare <- left_join(
  compute_region(agg_to_bands(results_cf1, "adm_base"), "base"),
  compute_region(agg_to_bands(results_cf1, "adm_cf"),   "cf1"),
  by = c("itl1_name","lad_imd_decile")
) %>%
  left_join(
    compute_region(agg_to_bands(results_cf2, "adm_cf"), "cf2"),
    by = c("itl1_name","lad_imd_decile")
  ) %>%
  pivot_longer(cols = c(adm_base, adm_cf1, adm_cf2),
               names_to = "scenario", values_to = "adm_per1000") %>%
  mutate(
    scenario = recode(scenario,
                      "adm_base" = "Baseline",
                      "adm_cf1"  = "CF1: decile 10 \u03c0\u2090",
                      "adm_cf2"  = "CF2: decile 5 \u03c0\u2090"),
    region_short = str_remove(itl1_name, " \\(England\\)")
  )

p2 <- ggplot(region_compare,
             aes(x = factor(lad_imd_decile), y = adm_per1000,
                 colour = scenario,
                 group  = interaction(region_short, scenario),
                 linetype = scenario)) +
  geom_line(linewidth = 0.85, alpha = 0.85) +
  geom_point(size = 1.4) +
  facet_wrap(~ region_short, nrow = 3) +
  scale_colour_manual(
    values = c("Baseline"                    = "#2166ac",
               "CF1: decile 10 \u03c0\u2090" = "#d73027",
               "CF2: decile 5 \u03c0\u2090"  = "#f4a582"),
    name = NULL) +
  scale_linetype_manual(
    values = c("Baseline"                    = "solid",
               "CF1: decile 10 \u03c0\u2090" = "solid",
               "CF2: decile 5 \u03c0\u2090"  = "dashed"),
    name = NULL) +
  scale_x_discrete(labels = c("1","","","","5","","","","","10")) +
  labs(
    title    = "COVID-19 admissions: baseline vs counterfactuals by ITL1 region",
    subtitle = "CF1 = decile 10 \u03c0\u2090 (red solid) | CF2 = decile 5 \u03c0\u2090 (orange dashed)",
    x        = "IMD deprivation decile",
    y        = "Cumulative admissions per 1,000",
    caption  = cf_caption
  ) +
  theme_pub +
  theme(strip.text      = element_text(face="bold", size=8.5),
        legend.position = "top",
        axis.text.x     = element_text(size=7))

ggsave("output/plots/counterfactual/16_fig2_burden_comparison_region.png",
       p2, width = 14, height = 10, dpi = 200)
cat("  Saved: 16_fig2_burden_comparison_region.png\n")

# Figure 3: National avoided burden by age group
cat("Plot 3: National summary by age group...\n")

age_summary <- results_all %>%
  group_by(scenario, age_group, age_idx) %>%
  summarise(
    adm_avoided_abs   = sum(adm_avoided_abs),
    death_avoided_abs = sum(death_avoided_abs),
    .groups = "drop"
  ) %>%
  mutate(age_group = factor(age_group, levels = age_labels),
         scenario  = factor(scenario,
                            levels = c("CF1: decile 10 \u03c0\u2090",
                                       "CF2: decile 5 \u03c0\u2090")))

p3 <- ggplot(age_summary,
             aes(x = adm_avoided_abs,
                 y = fct_reorder(age_group, age_idx),
                 fill = scenario)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.9) +
  geom_text(aes(label = comma(round(adm_avoided_abs))),
            position = position_dodge(width = 0.7),
            hjust = -0.1, size = 2.5, colour = "#333333") +
  scale_fill_manual(
    values = c("CF1: decile 10 \u03c0\u2090" = "#2166ac",
               "CF2: decile 5 \u03c0\u2090"  = "#74add1"),
    name = "Counterfactual"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.22)), labels = comma) +
  labs(
    title    = "Avoided hospital admissions by age group",
    subtitle = "CF1 (decile 10 \u03c0\u2090) vs CF2 (decile 5 \u03c0\u2090)",
    x        = "Avoided admissions (absolute count, blended population)",
    y        = "Age group",
    caption  = cf_caption
  ) +
  theme_pub +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(colour="#eeeeee"),
        legend.position    = "top")

ggsave("output/plots/counterfactual/16_fig3_national_summary.png",
       p3, width = 13, height = 7, dpi = 200)
cat("  Saved: 16_fig3_national_summary.png\n")
print(p3)

# Save results
write_csv(
  results_all %>%
    select(scenario, imd_decile, age_group, age_idx,
           adm_base, adm_cf, adm_diff, adm_attr_frac,
           death_base, death_cf, death_diff, death_attr_frac,
           adm_avoided_abs, death_avoided_abs) %>%
    arrange(scenario, imd_decile, age_idx),
  "output/plots/counterfactual/16_counterfactual_full.csv"
)

write_csv(decile_summary,
          "output/plots/counterfactual/16_counterfactual_decile_summary.csv")

cat("\n============================================================\n")
cat("National totals:\n")
for (sc in unique(results_all$scenario)) {
  r <- results_all %>% filter(scenario == sc)
  cat(sprintf("  %s\n", sc))
  cat(sprintf("    Avoided admissions: %s | deaths: %s | mean attr frac: %.1f%%\n",
              comma(round(sum(r$adm_avoided_abs))),
              comma(round(sum(r$death_avoided_abs))),
              mean(r$adm_attr_frac)))
}
cat("\nAll outputs saved to output/plots/counterfactual/\n")
