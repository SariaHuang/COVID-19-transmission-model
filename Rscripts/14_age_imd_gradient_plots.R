# ==============================================================
# Script: 14_age_imd_gradient_plots.R
#
# Purpose: Visualise model outputs by age group x IMD decile,
#          using posterior median beta from script 13.
#
# Produces:
#   (1) Heatmap: cumulative admissions per 1,000 (age x decile)
#   (2) Heatmap: cumulative deaths per 1,000 (age x decile)
#   (3) Gradient lines: admissions per 1,000 by IMD decile
#   (4) Gradient lines: deaths per 1,000 by IMD decile
#   (5) Age-stratified admission curves by decile (faceted)
#
# Inputs:
#   output/fitting/fitted_samples_imd{1-10}.rds  (script 13)
#   All model inputs loaded by scripts 08-09
# ==============================================================

library(odin2)
library(dust2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(purrr)

if (!exists("age_seird_hosp")) {
  SKIP09_RUN <- TRUE
  source("Rscripts/08_refine_gamma_hd_hr.R")
  source("Rscripts/09_age_imd_stratified_odin.R")
  rm(SKIP09_RUN)
}
stopifnot(exists("age_seird_hosp"))
cat("odin2 model loaded.\n")

age_labels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                "70 to 74","75+")

dir.create("output/plots", recursive = TRUE, showWarnings = FALSE)

# State indices in dust2 output (170 states: 10 compartments x 17 ages)
# S[1:17], E[18:34], Ip[35:51], Ic[52:68], Is[69:85],
# HD[86:102], HR[103:119], D[120:136], R[137:153], Adm[154:170]
adm_idx <- 154:170
d_idx   <- 120:136

# Run odin2 model, returning age-specific cumulative admissions and deaths
run_epidemic_fit_age <- function(imd_decile, beta, I0_frac = 1e-4,
                                 urban = TRUE) {
  contact <- as.matrix(read.csv(
    paste0("data/parameters/contact_matrix_imd", imd_decile, ".csv"),
    header = FALSE
  ))
  
  pi_a    <- pi_matrix[[paste0("imd_", imd_decile)]]
  h_a     <- h_mu$h_a
  mu_ca_h <- h_mu$mu_ca_h
  
  proportion <- rural_age %>%
    filter(IMD == imd_decile,
           rural == if (urban) "Urban" else "Rural") %>%
    arrange(Age) %>%
    pull(Proportion)
  
  pop_by_age <- rural_age %>%
    filter(IMD == imd_decile,
           rural == if (urban) "Urban" else "Rural") %>%
    arrange(Age) %>%
    pull(Population)
  
  S0    <- proportion
  S0[8] <- S0[8] - I0_frac
  Ip0   <- c(rep(0, 7), I0_frac, rep(0, 9))
  
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
  
  # Extract age-specific Adm and D from dust2 output
  # out[adm_idx[a], ] = cumulative admissions for age group a
  # out[d_idx[a], ]   = cumulative deaths for age group a
  map_dfr(seq_along(age_labels), function(a) {
    data.frame(
      day               = 0:365,
      imd_decile        = imd_decile,
      age_group         = age_labels[a],
      age_idx           = a,
      pop_urban         = pop_by_age[a],
      cum_adm_per1000   = out[adm_idx[a], ] * 1000,
      cum_death_per1000 = out[d_idx[a],   ] * 1000
    )
  }) %>%
    mutate(
      new_adm_per1000   = pmax(c(0, diff(cum_adm_per1000)),   0),
      new_death_per1000 = pmax(c(0, diff(cum_death_per1000)), 0)
    )
}

# Extract posterior median beta for each decile
burnin <- 500
cat("Extracting posterior median beta...\n")

beta_posterior <- sapply(1:10, function(d) {
  fit <- readRDS(paste0("output/fitting/fitted_samples_imd", d, ".rds"))
  median(sapply(1:3, function(ch)
    median(exp(fit$pars[1, (burnin+1):2000, ch]))))
})

for (d in 1:10) cat(sprintf("  Decile %2d: beta = %.5f\n", d, beta_posterior[d]))

# Run model for all 10 deciles
cat("\nRunning model for all 10 deciles...\n")
all_age_results <- map_dfr(1:10, function(d) {
  cat("  Decile", d, "\n")
  run_epidemic_fit_age(imd_decile = d, beta = beta_posterior[d])
})
cat("Total rows:", nrow(all_age_results), "\n")

# Final-day cumulative totals
final_day <- all_age_results %>%
  group_by(imd_decile, age_group, age_idx) %>%
  slice_max(day, n = 1) %>%
  ungroup()

# Plot 1: Heatmap -- cumulative admissions per 1,000
cat("\nPlot 1: Admission heatmap...\n")

p1 <- final_day %>%
  mutate(age_group = factor(age_group, levels = age_labels)) %>%
  ggplot(aes(x = factor(imd_decile), y = age_group,
             fill = cum_adm_per1000)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(cum_adm_per1000, 1),
                colour = cum_adm_per1000 > 5),
            size = 3) +
  scale_colour_manual(values = c("TRUE"="white","FALSE"="grey20"),
                      guide = "none") +
  scale_fill_distiller(palette = "RdYlBu", direction = -1,
                       name = "Cumulative\nadmissions\nper 1,000") +
  scale_x_discrete(labels = c("1\n(most\ndeprived)", 2:9,
                              "10\n(least\ndeprived)")) +
  labs(
    title    = "Cumulative COVID-19 hospital admissions per 1,000 population",
    subtitle = "By age group and IMD decile — fitted model (posterior median \u03b2)",
    x        = "IMD deprivation decile",
    y        = "Age group",
    caption  = "Model: age \u00d7 IMD SEIRD + hospital (odin2). Parameters: Goodfellow (2024) + Knock (2021)."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, colour = "#555555"),
    axis.text.y       = element_text(size = 8),
    legend.key.height = unit(1.2, "cm")
  )

ggsave("output/plots/14_heatmap_admissions.png", p1,
       width = 11, height = 7, dpi = 200)
cat("  Saved: 14_heatmap_admissions.png\n")
print(p1)

# Plot 2: Heatmap -- cumulative deaths per 1,000
cat("Plot 2: Death heatmap...\n")

p2 <- final_day %>%
  mutate(age_group = factor(age_group, levels = age_labels)) %>%
  ggplot(aes(x = factor(imd_decile), y = age_group,
             fill = cum_death_per1000)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = round(cum_death_per1000, 2),
                colour = cum_death_per1000 > 1.3),
            size = 3) +
  scale_colour_manual(values = c("TRUE"="white","FALSE"="grey20"),
                      guide = "none") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Cumulative\ndeaths\nper 1,000") +
  scale_x_discrete(labels = c("1\n(most\ndeprived)", 2:9,
                              "10\n(least\ndeprived)")) +
  labs(
    title    = "Cumulative COVID-19 in-hospital deaths per 1,000 population",
    subtitle = "By age group and IMD decile — fitted model (posterior median \u03b2)",
    x        = "IMD deprivation decile",
    y        = "Age group",
    caption  = "Model: age \u00d7 IMD SEIRD + hospital (odin2). Parameters: Goodfellow (2024) + Knock (2021)."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, colour = "#555555"),
    axis.text.y       = element_text(size = 8),
    legend.key.height = unit(1.2, "cm")
  )

ggsave("output/plots/14_heatmap_deaths.png", p2,
       width = 11, height = 7, dpi = 200)
cat("  Saved: 14_heatmap_deaths.png\n")

# Plot 3: Gradient lines -- admissions
cat("Plot 3: IMD gradient lines (admissions)...\n")

highlight_ages <- c("65 to 69","70 to 74","75+","55 to 59","60 to 64")

gradient_df <- final_day %>%
  mutate(
    age_group  = factor(age_group, levels = age_labels),
    highlight  = age_group %in% highlight_ages,
    line_alpha = ifelse(highlight, 1, 0.35),
    line_size  = ifelse(highlight, 1.1, 0.5)
  )

p3 <- ggplot(gradient_df,
             aes(x = imd_decile, y = cum_adm_per1000,
                 group = age_group, colour = age_group)) +
  geom_line(aes(alpha = line_alpha, linewidth = line_size)) +
  geom_point(data = filter(gradient_df, highlight), size = 1.8) +
  scale_x_continuous(breaks = 1:10,
                     labels = c("1\n(most\ndeprived)", 2:9,
                                "10\n(least\ndeprived)")) +
  scale_colour_viridis_d(option = "plasma", name = "Age group",
                         direction = -1) +
  scale_alpha_identity() +
  scale_linewidth_identity() +
  labs(
    title    = "Deprivation gradient in cumulative hospital admissions",
    subtitle = "Per 1,000 population by age group — posterior median \u03b2",
    x        = "IMD deprivation decile (1 = most deprived)",
    y        = "Cumulative admissions per 1,000",
    caption  = "Fitted model, posterior median \u03b2 per decile."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, colour = "#555555"),
    legend.text       = element_text(size = 8),
    legend.key.height = unit(0.4, "cm")
  )

ggsave("output/plots/14_gradient_admissions.png", p3,
       width = 11, height = 6.5, dpi = 200)
cat("  Saved: 14_gradient_admissions.png\n")

# Plot 4: Gradient lines -- deaths
cat("Plot 4: IMD gradient lines (deaths)...\n")

p4 <- ggplot(gradient_df,
             aes(x = imd_decile, y = cum_death_per1000,
                 group = age_group, colour = age_group)) +
  geom_line(aes(alpha = line_alpha, linewidth = line_size)) +
  geom_point(data = filter(gradient_df, highlight), size = 1.8) +
  scale_x_continuous(breaks = 1:10,
                     labels = c("1\n(most\ndeprived)", 2:9,
                                "10\n(least\ndeprived)")) +
  scale_colour_viridis_d(option = "inferno", name = "Age group",
                         direction = -1) +
  scale_alpha_identity() +
  scale_linewidth_identity() +
  labs(
    title    = "Deprivation gradient in cumulative in-hospital deaths",
    subtitle = "Per 1,000 population by age group — posterior median \u03b2",
    x        = "IMD deprivation decile (1 = most deprived)",
    y        = "Cumulative deaths per 1,000",
    caption  = "Fitted model, posterior median \u03b2 per decile."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, colour = "#555555"),
    legend.text       = element_text(size = 8),
    legend.key.height = unit(0.4, "cm")
  )

ggsave("output/plots/14_gradient_deaths.png", p4,
       width = 11, height = 6.5, dpi = 200)
cat("  Saved: 14_gradient_deaths.png\n")

# Plot 5: Age-stratified admission curves by decile
cat("Plot 5: Age-stratified admission curves...\n")

decile_labels <- setNames(
  paste0("IMD decile ", 1:10,
         ifelse(1:10 == 1, " (most deprived)",
                ifelse(1:10 == 10, " (least deprived)", ""))),
  1:10
)

p5 <- all_age_results %>%
  mutate(
    age_group  = factor(age_group, levels = age_labels),
    decile_lab = factor(decile_labels[as.character(imd_decile)],
                        levels = decile_labels)
  ) %>%
  ggplot(aes(x = day, y = new_adm_per1000,
             group = age_group, colour = age_group)) +
  geom_line(linewidth = 0.6, alpha = 0.85) +
  facet_wrap(~ decile_lab, nrow = 2) +
  coord_cartesian(ylim = c(0, NA)) +
  scale_colour_viridis_d(option = "plasma", name = "Age group",
                         direction = -1) +
  labs(
    title    = "Age-stratified daily hospital admissions by IMD deprivation decile",
    subtitle = "Per 1,000 population — fitted model, posterior median \u03b2",
    x        = "Day of epidemic",
    y        = "New admissions per 1,000",
    caption  = "Model: age \u00d7 IMD SEIRD + hospital (odin2). Parameters: Goodfellow (2024) + Knock (2021)."
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "#555555"),
    strip.text       = element_text(face = "bold", size = 8),
    legend.text      = element_text(size = 7),
    legend.key.height = unit(0.35, "cm"),
    panel.grid.minor = element_blank()
  )

ggsave("output/plots/14_age_curves_by_decile.png", p5,
       width = 16, height = 9, dpi = 200)
cat("  Saved: 14_age_curves_by_decile.png\n")

# Summary table
cat("\n--- Summary: cumulative burden by IMD decile ---\n")
summary_table <- final_day %>%
  group_by(imd_decile) %>%
  summarise(
    total_adm_per1000   = round(sum(cum_adm_per1000),   2),
    total_death_per1000 = round(sum(cum_death_per1000), 3),
    .groups = "drop"
  )
print(summary_table, n = 10)

write_csv(summary_table, "output/plots/14_summary_burden_by_decile.csv")

cat("\nAll outputs saved to output/plots/\n")
