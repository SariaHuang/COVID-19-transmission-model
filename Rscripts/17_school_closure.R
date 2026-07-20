# ==============================================================
# Script: 17_school_closure.R
#
# Purpose: School closure policy scenario analysis.
#          Extends Goodfellow et al. (2024) by applying school
#          closure to a model fitted to real NHS admissions data.
#          Results disaggregated by IMD decile, age group,
#          and ITL1 region.
#
# Method:
#   School contact fraction identified from POLYMOD UK data
#   (socialmixr, Mossong et al. 2008) using Goodfellow's
#   density-correction method. Fraction removed from decile-
#   specific model contact matrices.
#
#   M_closed[i,j] = M_total[i,j] * (1 - x * school_frac[i,j])
#   R0 = spectral radius of Next Generation Matrix using
#        decile-specific model contact matrices.
#
# Beta: fixed at decile 1 posterior median (0.031).
# Population: blended urban/rural (consistent with scripts 09/13).
#
# Outputs: output/plots/school_closure/
#   17_fig1_R0_reduction.png
#   17_fig2_admission_reduction.png
#   17_fig3_pct_reduction_by_decile.png
#   17_fig4_reduction_by_age.png
#   17_fig5_region_decile_facet.png
#   17_fig6_region_avoided.png
#   17_school_closure_decile_summary.csv
#   17_school_closure_full.csv
# ==============================================================

library(odin2)
library(dust2)
library(socialmixr)
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
cat(sprintf("Fixed beta: %.5f\n", beta_fixed))

dir.create("output/plots/school_closure", recursive = TRUE,
           showWarnings = FALSE)

adm_idx <- 154:170
d_idx   <- 120:136

# POLYMOD school fraction matrix
cat("\nGenerating POLYMOD school fraction matrix...\n")

data(polymod)
age_breaks <- c(0, 1, seq(5, 75, 5))

polymod_total <- socialmixr::contact_matrix(
  polymod, countries = "United Kingdom",
  age.limits = age_breaks, return.demography = TRUE, symmetric = TRUE
)
M_polymod <- polymod_total$matrix
N         <- polymod_total$demography$population

G <- matrix(nrow=17, ncol=17)
for (i in 1:17) for (j in 1:17) G[i,j] <- M_polymod[i,j] * sum(N) / N[j]

polymod_school <- socialmixr::contact_matrix(
  polymod, countries = "United Kingdom",
  age.limits = age_breaks, return.demography = TRUE, symmetric = TRUE,
  filter = list(cnt_school = 1)
)
M_school_polymod <- polymod_school$matrix

G_school <- matrix(nrow=17, ncol=17)
for (i in 1:17) for (j in 1:17) {
  G_school[i,j] <- M_school_polymod[i,j] * sum(N) / N[j]
}

school_frac <- pmin(G_school / (G + 1e-10), 1)
cat(sprintf("  School contacts: %.1f%% of total\n", sum(G_school)/sum(G)*100))

make_school_closure_matrix <- function(imd_decile, x = 1) {
  M_total <- as.matrix(read.csv(
    paste0("data/parameters/contact_matrix_imd", imd_decile, ".csv"),
    header = FALSE))
  pmax(M_total * (1 - x * school_frac), 0)
}

compute_R0 <- function(M, pi_a, beta,
                       sympt=1/2.1, rec_c=2.9, rec_s=5, xi=0.5) {
  ip  <- pi_a*(1/sympt+rec_c) + xi*(1-pi_a)*rec_s
  NGM <- beta * M * outer(rep(1,17), ip)
  max(Re(eigen(NGM, only.values=TRUE)$values))
}

cat("\nComputing R0 for all deciles...\n")
r0_results <- map_dfr(1:10, function(d) {
  pi_a    <- pi_matrix[[paste0("imd_", d)]]
  M_base  <- as.matrix(read.csv(
    paste0("data/parameters/contact_matrix_imd", d, ".csv"), header=FALSE))
  M_close <- make_school_closure_matrix(d)
  r0_b <- compute_R0(M_base,  pi_a, beta_fixed)
  r0_c <- compute_R0(M_close, pi_a, beta_fixed)
  data.frame(imd_decile   = d,
             R0_baseline  = round(r0_b, 3),
             R0_closed    = round(r0_c, 3),
             R0_reduction = round(r0_b - r0_c, 3),
             R0_pct_red   = round((r0_b - r0_c)/r0_b*100, 1))
})
cat("\n--- R0 results ---\n")
print(r0_results)

# Run odin2 model using blended population
run_epidemic_school <- function(imd_decile, beta, school_closed = FALSE) {
  contact <- if (school_closed) make_school_closure_matrix(imd_decile) else
    as.matrix(read.csv(
      paste0("data/parameters/contact_matrix_imd", imd_decile, ".csv"),
      header = FALSE))
  
  pi_a    <- pi_matrix[[paste0("imd_", imd_decile)]]
  h_a     <- h_mu$h_a
  mu_ca_h <- h_mu$mu_ca_h
  
  blended    <- get_blended_inputs(imd_decile)
  proportion <- blended$proportion
  pop_by_age <- proportion * blended$pop_size
  
  S0 <- proportion; S0[8] <- S0[8] - 1e-4
  Ip0 <- c(rep(0,7), 1e-4, rep(0,9))
  
  sys <- dust2::dust_system_create(age_seird_hosp, list(
    S0=S0, Ip0=Ip0, proportion=proportion,
    pi_a=pi_a, h_a=h_a, mu_ca_h=mu_ca_h,
    contact=contact, gam_hd=gamma_hd_vec, gam_hr=gamma_hr_vec,
    susc=beta))
  dust2::dust_system_set_state_initial(sys)
  out <- dust2::dust_system_simulate(sys, seq(0, 365, by=1))
  
  map_dfr(seq_along(age_labels), function(a) {
    data.frame(day=0:365, imd_decile=imd_decile,
               age_group=age_labels[a], age_idx=a,
               pop_blended   = pop_by_age[a],
               cum_adm_per1000   = out[adm_idx[a], ] * 1000,
               cum_death_per1000 = out[d_idx[a],   ] * 1000,
               school_closed=school_closed)
  })
}

cat("\nRunning ODE: baseline...\n")
baseline_ode <- map_dfr(1:10, function(d) {
  cat("  Decile", d, "\n")
  run_epidemic_school(d, beta_fixed, school_closed=FALSE)
})

cat("Running ODE: school closure...\n")
closure_ode <- map_dfr(1:10, function(d) {
  cat("  Decile", d, "\n")
  run_epidemic_school(d, beta_fixed, school_closed=TRUE)
})

final_base_sc <- baseline_ode %>%
  group_by(imd_decile, age_group, age_idx, pop_blended) %>%
  slice_max(day, n=1) %>% ungroup() %>%
  rename(adm_base=cum_adm_per1000, death_base=cum_death_per1000) %>%
  select(-school_closed)

final_closed_sc <- closure_ode %>%
  group_by(imd_decile, age_group, age_idx) %>%
  slice_max(day, n=1) %>% ungroup() %>%
  rename(adm_closed=cum_adm_per1000, death_closed=cum_death_per1000) %>%
  select(-school_closed, -pop_blended)

sc_results <- left_join(final_base_sc, final_closed_sc,
                        by=c("imd_decile","age_group","age_idx")) %>%
  mutate(
    adm_reduction     = adm_base   - adm_closed,
    death_reduction   = death_base - death_closed,
    adm_pct_red       = pmax(adm_reduction   / adm_base   * 100, 0),
    death_pct_red     = pmax(death_reduction / death_base * 100, 0),
    adm_avoided_abs   = adm_reduction   / 1000 * pop_blended,
    death_avoided_abs = death_reduction / 1000 * pop_blended
  )

decile_sc_summary <- sc_results %>%
  group_by(imd_decile) %>%
  summarise(
    adm_base_total     = sum(adm_base),
    adm_closed_total   = sum(adm_closed),
    adm_pct_red_mean   = round(mean(adm_pct_red), 1),
    death_pct_red_mean = round(mean(death_pct_red), 1),
    adm_avoided_abs    = round(sum(adm_avoided_abs)),
    death_avoided_abs  = round(sum(death_avoided_abs)),
    .groups = "drop"
  ) %>%
  left_join(r0_results, by="imd_decile")

# Regional population
cat("\nLoading regional population data...\n")
region_pop <- read_csv("data/processed/population_age_imd_region.csv",
                       show_col_types=FALSE)

region_pop_agg <- region_pop %>%
  mutate(age_band_model = case_when(
    age_band %in% c("75-79","80+") ~ "75+", TRUE ~ age_band)) %>%
  group_by(itl1_name, lad_imd_decile, age_band_model) %>%
  summarise(population=sum(population), .groups="drop")

region_total <- region_pop_agg %>%
  group_by(itl1_name) %>%
  summarise(total_pop=sum(population), .groups="drop")

age_map <- tibble(
  age_group      = age_labels,
  age_band_model = c("0-4","0-4","5-9","10-14","15-19","20-24",
                     "25-29","30-34","35-39","40-44","45-49","50-54",
                     "55-59","60-64","65-69","70-74","75+"))

pop_w <- rural_age %>%
  filter(rural=="Urban", Age %in% c("Under 1","1 to 4")) %>%
  group_by(Age) %>% summarise(pop=sum(Population), .groups="drop") %>%
  mutate(w=pop/sum(pop))
w_u1 <- pop_w$w[pop_w$Age == "Under 1"]
w_14 <- pop_w$w[pop_w$Age == "1 to 4"]

agg_sc <- function(df, col) {
  df %>%
    left_join(age_map, by="age_group") %>%
    group_by(imd_decile, age_band_model) %>%
    summarise(
      value = if (n()==2)
        weighted.mean(.data[[col]],
                      c(w_u1,w_14)[match(age_group,c("Under 1","1 to 4"))])
      else first(.data[[col]]),
      .groups="drop") %>%
    rename(!!col := value)
}

rates_base_sc   <- agg_sc(sc_results, "adm_base")
rates_closed_sc <- agg_sc(sc_results, "adm_closed")

compute_region_sc <- function(rates, col) {
  region_pop_agg %>%
    left_join(rates, by=c("lad_imd_decile"="imd_decile","age_band_model")) %>%
    mutate(abs_adm=(.data[[col]]/1000)*population) %>%
    group_by(itl1_name, lad_imd_decile) %>%
    summarise(abs_adm=sum(abs_adm, na.rm=TRUE), .groups="drop") %>%
    left_join(region_total, by="itl1_name") %>%
    mutate(adm_per1000=abs_adm/total_pop*1000) %>%
    select(itl1_name, lad_imd_decile, adm_per1000)
}

region_sc_compare <- left_join(
  compute_region_sc(rates_base_sc,   "adm_base")   %>% rename(adm_base   = adm_per1000),
  compute_region_sc(rates_closed_sc, "adm_closed") %>% rename(adm_closed = adm_per1000),
  by=c("itl1_name","lad_imd_decile")) %>%
  mutate(region_short = str_remove(itl1_name," \\(England\\)"))

# Shared theme and colours
theme_pub <- theme_minimal(base_size=11) +
  theme(
    plot.title    = element_text(face="bold", size=12, margin=margin(b=4)),
    plot.subtitle = element_text(size=9, colour="#444444", margin=margin(b=8)),
    plot.caption  = element_text(size=7.5, colour="#888888",
                                 hjust=0, margin=margin(t=8)),
    axis.title    = element_text(size=9.5),
    axis.text     = element_text(size=8.5, colour="#333333"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour="#eeeeee"),
    legend.title  = element_text(size=8.5, face="bold"),
    legend.text   = element_text(size=8),
    plot.margin   = margin(12,12,8,12))

region_colours <- c(
  "East"                     = "#1b9e77",
  "East Midlands (England)"  = "#d95f02",
  "London"                   = "#7570b3",
  "North East (England)"     = "#e7298a",
  "North West (England)"     = "#66a61e",
  "South East (England)"     = "#e6ab02",
  "South West (England)"     = "#a6761d",
  "West Midlands (England)"  = "#666666",
  "Yorkshire and The Humber" = "#1f78b4")

sc_caption <- paste0(
  "School closure: POLYMOD UK school contact fraction removed from ",
  "decile-specific contact matrices (Mossong et al. 2008, via socialmixr), ",
  "following Goodfellow et al. (2024). \u03b2 fixed at decile 1 posterior ",
  "median (0.031). Population: blended urban/rural per decile. ",
  "Model: age \u00d7 IMD SEIRD + hospital (odin2).")

# Figure 1: R0
cat("\nPlot 1: R0 by decile...\n")
fig1_data <- decile_sc_summary %>%
  select(imd_decile, R0_baseline, R0_closed) %>%
  pivot_longer(cols=c(R0_baseline,R0_closed),
               names_to="scenario", values_to="R0") %>%
  mutate(scenario=recode(scenario,"R0_baseline"="Baseline",
                         "R0_closed"="Schools closed"))

p1 <- ggplot(fig1_data, aes(x=imd_decile, y=R0, colour=scenario, group=scenario)) +
  geom_line(linewidth=1.1) + geom_point(size=2.5) +
  geom_hline(yintercept=1, linetype="dashed", colour="#cc0000", alpha=0.6) +
  annotate("text", x=0.7, y=1.05, label="R\u2080 = 1",
           colour="#cc0000", size=3, hjust=0) +
  scale_x_continuous(breaks=1:10,
                     labels=c("1\n(most\ndeprived)",2:9,"10\n(least\ndeprived)")) +
  scale_colour_manual(values=c("Baseline"="steelblue","Schools closed"="#d73027"),
                      name=NULL) +
  labs(title="Effect of school closure on R\u2080 by IMD deprivation decile",
       subtitle="POLYMOD UK school contact fraction removed \u2014 \u03b2 fixed at 0.031",
       x="IMD deprivation decile", y="Basic reproduction number (R\u2080)",
       caption=sc_caption) +
  theme_pub + theme(legend.position="top")

ggsave("output/plots/school_closure/17_fig1_R0_reduction.png",
       p1, width=11, height=6.5, dpi=200)
cat("  Saved: 17_fig1_R0_reduction.png\n")
print(p1)

# Figure 2: Admissions
cat("Plot 2: Admissions by decile...\n")
fig2_data <- decile_sc_summary %>%
  select(imd_decile, adm_base_total, adm_closed_total) %>%
  pivot_longer(cols=c(adm_base_total,adm_closed_total),
               names_to="scenario", values_to="adm_per1000") %>%
  mutate(scenario=recode(scenario,"adm_base_total"="Baseline",
                         "adm_closed_total"="Schools closed"))

p2 <- ggplot(fig2_data, aes(x=imd_decile, y=adm_per1000,
                            colour=scenario, group=scenario)) +
  geom_line(linewidth=1.1) + geom_point(size=2.5) +
  scale_x_continuous(breaks=1:10,
                     labels=c("1\n(most\ndeprived)",2:9,"10\n(least\ndeprived)")) +
  scale_colour_manual(values=c("Baseline"="steelblue","Schools closed"="#d73027"),
                      name=NULL) +
  labs(title="Cumulative hospital admissions per 1,000: baseline vs school closure",
       subtitle="By IMD deprivation decile \u2014 \u03b2 fixed at 0.031",
       x="IMD deprivation decile", y="Cumulative admissions per 1,000 population",
       caption=sc_caption) +
  theme_pub + theme(legend.position="top")

ggsave("output/plots/school_closure/17_fig2_admission_reduction.png",
       p2, width=11, height=6.5, dpi=200)
cat("  Saved: 17_fig2_admission_reduction.png\n")

# Figure 3: % reduction by decile
cat("Plot 3: % reduction by decile...\n")
p3 <- decile_sc_summary %>%
  ggplot(aes(x=imd_decile, y=adm_pct_red_mean)) +
  geom_col(fill="#d73027", alpha=0.85, width=0.65) +
  geom_text(aes(label=paste0(adm_pct_red_mean,"%")),
            vjust=-0.4, size=3.2, colour="#333333") +
  scale_x_continuous(breaks=1:10,
                     labels=c("1\n(most\ndeprived)",2:9,"10\n(least\ndeprived)")) +
  scale_y_continuous(expand=expansion(mult=c(0,0.12)),
                     labels=function(x) paste0(x,"%")) +
  labs(title="Reduction in hospital admissions from school closure by IMD decile",
       subtitle="Mean % reduction across age groups \u2014 \u03b2 fixed at 0.031",
       x="IMD deprivation decile", y="Reduction in cumulative admissions (%)",
       caption=sc_caption) +
  theme_pub

ggsave("output/plots/school_closure/17_fig3_pct_reduction_by_decile.png",
       p3, width=11, height=6.5, dpi=200)
cat("  Saved: 17_fig3_pct_reduction_by_decile.png\n")

# Figure 4: % reduction by age group
cat("Plot 4: % reduction by age group...\n")
age_sc_summary <- sc_results %>%
  group_by(age_group, age_idx) %>%
  summarise(adm_pct_red_mean=mean(adm_pct_red), .groups="drop") %>%
  mutate(age_group=factor(age_group, levels=age_labels))

p4 <- ggplot(age_sc_summary, aes(x=adm_pct_red_mean,
                                 y=reorder(age_group, age_idx))) +
  geom_col(fill="#4393c3", alpha=0.85) +
  geom_text(aes(label=paste0(round(adm_pct_red_mean,1),"%")),
            hjust=-0.1, size=2.8, colour="#333333") +
  scale_x_continuous(expand=expansion(mult=c(0,0.15)),
                     labels=function(x) paste0(x,"%")) +
  labs(title="School closure: reduction in admissions by age group",
       subtitle="National average across all IMD deciles \u2014 \u03b2 fixed at 0.031",
       x="Reduction in cumulative admissions (%)", y="Age group",
       caption=sc_caption) +
  theme_pub +
  theme(panel.grid.major.y=element_blank(),
        panel.grid.major.x=element_line(colour="#eeeeee"))

ggsave("output/plots/school_closure/17_fig4_reduction_by_age.png",
       p4, width=11, height=7, dpi=200)
cat("  Saved: 17_fig4_reduction_by_age.png\n")

# Figure 5: Region x decile facet
cat("Plot 5: Region x decile (facet)...\n")
p5_facet <- region_sc_compare %>%
  left_join(region_pop_agg %>%
              group_by(itl1_name, lad_imd_decile) %>%
              summarise(cell_pop=sum(population), .groups="drop"),
            by=c("itl1_name","lad_imd_decile")) %>%
  filter(cell_pop >= 10000) %>%
  pivot_longer(cols=c(adm_base,adm_closed),
               names_to="scenario", values_to="adm_per1000") %>%
  mutate(scenario=recode(scenario,"adm_base"="Baseline",
                         "adm_closed"="Schools closed"),
         region_short=str_remove(itl1_name," \\(England\\)")) %>%
  ggplot(aes(x=factor(lad_imd_decile), y=adm_per1000,
             colour=scenario, group=scenario)) +
  geom_line(linewidth=0.9, alpha=0.85) + geom_point(size=1.5) +
  facet_wrap(~ region_short, nrow=3) +
  scale_colour_manual(values=c("Baseline"="steelblue","Schools closed"="#d73027"),
                      name=NULL) +
  scale_x_discrete(labels=c("1","","","","5","","","","","10")) +
  labs(title="School closure effect on admissions by ITL1 region and IMD decile",
       subtitle="Blue = baseline, red = schools closed \u2014 \u03b2 fixed at 0.031",
       x="IMD deprivation decile (1 = most deprived)",
       y="Cumulative admissions per 1,000", caption=sc_caption) +
  theme_pub +
  theme(strip.text=element_text(face="bold", size=8.5), legend.position="top")

ggsave("output/plots/school_closure/17_fig5_region_decile_facet.png",
       p5_facet, width=16, height=10, dpi=200)
cat("  Saved: 17_fig5_region_decile_facet.png\n")
print(p5_facet)

# Figure 6: Avoided admissions by region
cat("Plot 6: Avoided admissions by ITL1 region...\n")
region_avoided <- region_sc_compare %>%
  left_join(region_total, by="itl1_name") %>%
  mutate(abs_base    = adm_base   / 1000 * total_pop,
         abs_closed  = adm_closed / 1000 * total_pop,
         abs_avoided = abs_base - abs_closed) %>%
  group_by(itl1_name) %>%
  summarise(total_avoided = sum(abs_avoided, na.rm=TRUE),
            total_base    = sum(abs_base,    na.rm=TRUE),
            pct_red       = total_avoided/total_base*100,
            region_short  = first(str_remove(itl1_name," \\(England\\)")),
            .groups="drop")

p6 <- ggplot(region_avoided,
             aes(x=reorder(region_short,-total_avoided),
                 y=total_avoided, fill=itl1_name)) +
  geom_col(width=0.65, alpha=0.9) +
  geom_text(aes(label=paste0(round(pct_red,1),"% reduction")),
            vjust=-0.4, size=2.8, colour="#333333") +
  scale_fill_manual(values=region_colours, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.12)), labels=comma) +
  labs(title="Avoided hospital admissions from school closure by ITL1 region",
       subtitle="Absolute count and % reduction \u2014 \u03b2 fixed at 0.031",
       x="ITL1 Region", y="Avoided admissions (absolute count)",
       caption=sc_caption) +
  theme_pub + theme(axis.text.x=element_text(angle=28, hjust=1, size=9))

ggsave("output/plots/school_closure/17_fig6_region_avoided.png",
       p6, width=12, height=7, dpi=200)
cat("  Saved: 17_fig6_region_avoided.png\n")
print(p6)

# Save tables
write_csv(decile_sc_summary,
          "output/plots/school_closure/17_school_closure_decile_summary.csv")
write_csv(sc_results %>%
            select(imd_decile, age_group, age_idx,
                   adm_base, adm_closed, adm_reduction, adm_pct_red,
                   death_base, death_closed, death_reduction, death_pct_red,
                   adm_avoided_abs, death_avoided_abs) %>%
            arrange(imd_decile, age_idx),
          "output/plots/school_closure/17_school_closure_full.csv")

cat("\n============================================================\n")
cat("School closure national totals:\n")
cat(sprintf("  Total avoided admissions: %s\n",
            comma(round(sum(sc_results$adm_avoided_abs)))))
cat(sprintf("  Total avoided deaths:     %s\n",
            comma(round(sum(sc_results$death_avoided_abs)))))
cat(sprintf("  Mean R0 reduction:        %.3f (%.1f%%)\n",
            mean(r0_results$R0_reduction), mean(r0_results$R0_pct_red)))
cat(sprintf("  Mean admission reduction: %.1f%%\n",
            mean(decile_sc_summary$adm_pct_red_mean)))
cat("\nAll outputs saved to output/plots/school_closure/\n")
