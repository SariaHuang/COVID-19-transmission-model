# ==============================================================
# Script: 15_region_age_imd_plots.R
#
# Purpose: Region x age x IMD deprivation plots.
#          Produces BOTH main analysis and sensitivity analysis
#          plots with unified color scales for direct comparison.
#
# Main analysis:        15_main_*.png    (fitted beta per decile)
# Sensitivity analysis: 15_sensitivity_* (fixed beta = decile 1)
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

if (!exists("age_seird_hosp")) {
  SKIP09_RUN <- TRUE
  source("Rscripts/08_refine_gamma_hd_hr.R")
  source("Rscripts/09_age_imd_stratified_odin.R")
  rm(SKIP09_RUN)
  cat("odin2 model loaded.\n")
}

if (!exists("age_labels")) {
  age_labels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                  "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                  "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                  "70 to 74","75+")
}

if (!exists("run_epidemic_fit_age")) {
  source("Rscripts/14_age_imd_gradient_plots.R")
}

if (!exists("beta_posterior")) {
  burnin <- 500
  beta_posterior <- sapply(1:10, function(d) {
    fit <- readRDS(paste0("output/fitting/fitted_samples_imd", d, ".rds"))
    median(sapply(1:3, function(ch)
      median(exp(fit$pars[1, (burnin+1):2000, ch]))))
  })
}

dir.create("output/plots/region", recursive = TRUE, showWarnings = FALSE)

# Load and prepare regional population data
region_pop <- read_csv("data/processed/population_age_imd_region.csv",
                       show_col_types = FALSE)

region_pop_clean <- region_pop %>%
  mutate(age_band_model = case_when(
    age_band %in% c("75-79", "80+") ~ "75+",
    TRUE ~ age_band
  )) %>%
  group_by(itl1_name, lad_imd_decile, age_band_model) %>%
  summarise(population = sum(population), .groups = "drop")

region_total_pop <- region_pop_clean %>%
  group_by(itl1_name) %>%
  summarise(total_pop = sum(population), .groups = "drop")

age_band_map <- tibble(
  model_age = age_labels,
  age_band_model = c("0-4","0-4","5-9","10-14","15-19","20-24",
                     "25-29","30-34","35-39","40-44","45-49","50-54",
                     "55-59","60-64","65-69","70-74","75+")
)

pop_w <- rural_age %>%
  filter(rural == "Urban", Age %in% c("Under 1","1 to 4")) %>%
  group_by(Age) %>%
  summarise(pop = sum(Population), .groups = "drop") %>%
  mutate(w = pop / sum(pop))
w_u1 <- pop_w$w[pop_w$Age == "Under 1"]
w_14 <- pop_w$w[pop_w$Age == "1 to 4"]

# Run model + compute regional burden
compute_region_burden <- function(beta_values) {
  
  all_age <- map_dfr(1:10, function(d) {
    run_epidemic_fit_age(imd_decile = d, beta = beta_values[d])
  })
  
  fd <- all_age %>%
    group_by(imd_decile, age_group, age_idx) %>%
    slice_max(day, n = 1) %>%
    ungroup()
  
  rates_agg <- fd %>%
    select(imd_decile, model_age = age_group,
           cum_adm_per1000, cum_death_per1000) %>%
    left_join(age_band_map, by = "model_age") %>%
    group_by(imd_decile, age_band_model) %>%
    summarise(
      cum_adm_per1000 = if (n() == 2)
        weighted.mean(cum_adm_per1000,
                      c(w_u1, w_14)[match(model_age, c("Under 1","1 to 4"))])
      else first(cum_adm_per1000),
      cum_death_per1000 = if (n() == 2)
        weighted.mean(cum_death_per1000,
                      c(w_u1, w_14)[match(model_age, c("Under 1","1 to 4"))])
      else first(cum_death_per1000),
      .groups = "drop"
    )
  
  burden <- region_pop_clean %>%
    left_join(rates_agg,
              by = c("lad_imd_decile" = "imd_decile",
                     "age_band_model")) %>%
    mutate(abs_adm   = (cum_adm_per1000   / 1000) * population,
           abs_death = (cum_death_per1000 / 1000) * population)
  
  region_decile <- burden %>%
    group_by(itl1_name, lad_imd_decile) %>%
    summarise(abs_adm   = sum(abs_adm,   na.rm = TRUE),
              abs_death = sum(abs_death, na.rm = TRUE),
              .groups = "drop") %>%
    left_join(region_total_pop, by = "itl1_name") %>%
    mutate(adm_per1000   = abs_adm   / total_pop * 1000,
           death_per1000 = abs_death / total_pop * 1000)
  
  region_decile_age <- burden %>%
    left_join(region_total_pop, by = "itl1_name") %>%
    mutate(adm_per1000   = abs_adm   / total_pop * 1000,
           death_per1000 = abs_death / total_pop * 1000)
  
  list(region_decile     = region_decile,
       region_decile_age = region_decile_age)
}

# Run both analyses
cat("Running main analysis (fitted beta)...\n")
results_main <- compute_region_burden(beta_posterior)

cat("Running sensitivity analysis (fixed beta)...\n")
results_sens <- compute_region_burden(rep(beta_posterior[1], 10))

# Unified color scale limits (from combined data)
combined_age <- bind_rows(
  results_main$region_decile_age,
  results_sens$region_decile_age
)

adm_cap   <- quantile(combined_age$adm_per1000,   0.95, na.rm = TRUE)
death_cap <- quantile(combined_age$death_per1000, 0.95, na.rm = TRUE)

cat(sprintf("Unified color scale caps: admissions = %.4f, deaths = %.5f\n",
            adm_cap, death_cap))

# Shared plot settings
region_colours <- c(
  "East"                     = "#1b9e77",
  "East Midlands (England)"  = "#d95f02",
  "London"                   = "#7570b3",
  "North East (England)"     = "#e7298a",
  "North West (England)"     = "#66a61e",
  "South East (England)"     = "#e6ab02",
  "South West (England)"     = "#a6761d",
  "West Midlands (England)"  = "#666666",
  "Yorkshire and The Humber" = "#1f78b4"
)

age_band_order <- c("0-4","5-9","10-14","15-19","20-24","25-29",
                    "30-34","35-39","40-44","45-49","50-54","55-59",
                    "60-64","65-69","70-74","75+")

theme_publication <- theme_minimal(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 12,
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 8.5, colour = "#444444",
                                    margin = margin(b = 8)),
    plot.caption     = element_text(size = 7, colour = "#888888",
                                    hjust = 0, margin = margin(t = 8)),
    strip.text       = element_text(face = "bold", size = 8.5),
    strip.background = element_rect(fill = "#f5f5f5", colour = NA),
    axis.text        = element_text(size = 7.5, colour = "#333333"),
    axis.title       = element_text(size = 9),
    legend.title     = element_text(size = 8, face = "bold"),
    legend.text      = element_text(size = 7.5),
    legend.key.height = unit(1.1, "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "#eeeeee"),
    plot.margin      = margin(12, 12, 8, 12)
  )

# Plot function (called once per analysis)
make_plots <- function(results, file_prefix, title_tag) {
  
  rd  <- results$region_decile
  rda <- results$region_decile_age
  
  caption_txt <- paste0(
    "Model: age \u00d7 IMD SEIRD + hospital pathway (odin2). ",
    "Parameters: Goodfellow et al. (2024) + Knock et al. (2021). ",
    title_tag, "."
  )
  
  # Plot 1: Gradient admissions
  p1 <- ggplot(rd, aes(x = lad_imd_decile, y = adm_per1000,
                       colour = itl1_name, group = itl1_name)) +
    geom_line(linewidth = 1.0, alpha = 0.9) +
    geom_point(size = 2.0) +
    scale_x_continuous(breaks = 1:10,
                       labels = c("1\n(most\ndeprived)", 2:9,
                                  "10\n(least\ndeprived)")) +
    scale_colour_manual(values = region_colours, name = "ITL1 Region") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = "Cumulative COVID-19 hospital admissions by IMD decile and ITL1 region",
      subtitle = title_tag,
      x        = "IMD deprivation decile",
      y        = "Cumulative admissions per 1,000 population",
      caption  = caption_txt
    ) +
    theme_publication +
    theme(legend.position = "right",
          legend.key.height = unit(0.5, "cm"))
  
  ggsave(paste0("output/plots/region/", file_prefix, "_gradient_adm.png"),
         p1, width = 12, height = 6.5, dpi = 200)
  
  # Plot 2: Gradient deaths
  p2 <- ggplot(rd, aes(x = lad_imd_decile, y = death_per1000,
                       colour = itl1_name, group = itl1_name)) +
    geom_line(linewidth = 1.0, alpha = 0.9) +
    geom_point(size = 2.0) +
    scale_x_continuous(breaks = 1:10,
                       labels = c("1\n(most\ndeprived)", 2:9,
                                  "10\n(least\ndeprived)")) +
    scale_colour_manual(values = region_colours, name = "ITL1 Region") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(
      title    = "Cumulative COVID-19 in-hospital deaths by IMD decile and ITL1 region",
      subtitle = title_tag,
      x        = "IMD deprivation decile",
      y        = "Cumulative deaths per 1,000 population",
      caption  = caption_txt
    ) +
    theme_publication +
    theme(legend.position = "right",
          legend.key.height = unit(0.5, "cm"))
  
  ggsave(paste0("output/plots/region/", file_prefix, "_gradient_death.png"),
         p2, width = 12, height = 6.5, dpi = 200)
  
  # Plot 3: Heatmap admissions (unified color scale)
  p3 <- rda %>%
    mutate(
      age_band_model = factor(age_band_model, levels = age_band_order),
      region_short   = str_remove(itl1_name, " \\(England\\)")
    ) %>%
    ggplot(aes(x = factor(lad_imd_decile),
               y = age_band_model,
               fill = adm_per1000)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    facet_wrap(~ region_short, nrow = 3) +
    scale_fill_distiller(
      palette   = "RdYlBu",
      direction = -1,
      name      = "Admissions\nper 1,000",
      limits    = c(0, adm_cap),
      oob       = scales::squish,
      breaks    = pretty(c(0, adm_cap), n = 5),
      labels    = number_format(accuracy = 0.001)
    ) +
    scale_x_discrete(labels = c("1", 2:9, "10")) +
    labs(
      title    = "Cumulative COVID-19 admissions per 1,000 by region, age and IMD decile",
      subtitle = title_tag,
      x        = "IMD decile (1 = most deprived)",
      y        = "Age band",
      caption  = caption_txt
    ) +
    theme_publication +
    theme(axis.text.x      = element_text(size = 6.5),
          axis.text.y      = element_text(size = 7),
          panel.grid.major = element_blank())
  
  ggsave(paste0("output/plots/region/", file_prefix, "_heatmap_adm.png"),
         p3, width = 16, height = 10, dpi = 200)
  
  # Plot 4: Heatmap deaths (unified color scale)
  p4 <- rda %>%
    mutate(
      age_band_model = factor(age_band_model, levels = age_band_order),
      region_short   = str_remove(itl1_name, " \\(England\\)")
    ) %>%
    ggplot(aes(x = factor(lad_imd_decile),
               y = age_band_model,
               fill = death_per1000)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    facet_wrap(~ region_short, nrow = 3) +
    scale_fill_distiller(
      palette   = "YlOrRd",
      direction = 1,
      name      = "Deaths\nper 1,000",
      limits    = c(0, death_cap),
      oob       = scales::squish,
      breaks    = pretty(c(0, death_cap), n = 5),
      labels    = number_format(accuracy = 0.0001)
    ) +
    scale_x_discrete(labels = c("1", 2:9, "10")) +
    labs(
      title    = "Cumulative COVID-19 in-hospital deaths per 1,000 by region, age and IMD decile",
      subtitle = title_tag,
      x        = "IMD decile (1 = most deprived)",
      y        = "Age band",
      caption  = caption_txt
    ) +
    theme_publication +
    theme(axis.text.x      = element_text(size = 6.5),
          axis.text.y      = element_text(size = 7),
          panel.grid.major = element_blank())
  
  ggsave(paste0("output/plots/region/", file_prefix, "_heatmap_death.png"),
         p4, width = 16, height = 10, dpi = 200)
  
  # Plot 5: Stacked bar
  p5 <- rd %>%
    mutate(
      region_short = str_remove(itl1_name, " \\(England\\)"),
      decile_f     = factor(lad_imd_decile,
                            labels = paste0("Decile ", 1:10))
    ) %>%
    group_by(region_short) %>%
    mutate(total = sum(adm_per1000)) %>%
    ungroup() %>%
    ggplot(aes(x = reorder(region_short, -total),
               y = adm_per1000,
               fill = factor(lad_imd_decile))) +
    geom_col(position = "stack", width = 0.65) +
    geom_text(aes(label = ifelse(adm_per1000 > 0.02,
                                 round(adm_per1000, 2), "")),
              position = position_stack(vjust = 0.5),
              size = 2.2, colour = "white", fontface = "bold") +
    scale_fill_brewer(palette = "RdYlBu", direction = -1,
                      name = "IMD decile\n(1 = most deprived)") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(
      title    = "Cumulative COVID-19 hospital admissions per 1,000 by ITL1 region",
      subtitle = paste0("Stacked by IMD deprivation decile \u2014 ", title_tag),
      x        = "ITL1 Region",
      y        = "Cumulative admissions per 1,000 population",
      caption  = caption_txt
    ) +
    theme_publication +
    theme(axis.text.x       = element_text(angle = 28, hjust = 1, size = 8.5),
          legend.key.height = unit(0.45, "cm"))
  
  ggsave(paste0("output/plots/region/", file_prefix, "_bar_stacked.png"),
         p5, width = 12, height = 7, dpi = 200)
  
  write_csv(rd, paste0("output/plots/region/", file_prefix, "_summary.csv"))
  cat("  Saved 5 plots + CSV with prefix:", file_prefix, "\n")
}

# Generate all plots
cat("\n--- Main analysis plots ---\n")
make_plots(
  results     = results_main,
  file_prefix = "15_main",
  title_tag   = "Main analysis: decile-specific fitted \u03b2"
)

cat("\n--- Sensitivity analysis plots ---\n")
make_plots(
  results     = results_sens,
  file_prefix = "15_sensitivity",
  title_tag   = paste0("Sensitivity analysis: uniform \u03b2 = ",
                       round(beta_posterior[1], 4),
                       " (decile 1 posterior median)")
)

cat("\nDone. All 10 plots saved to output/plots/region/\n")
