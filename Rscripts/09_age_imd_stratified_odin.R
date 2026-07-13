# ==============================================================
# Script: 09_age_imd_stratified_odin.R
#
# Purpose: Age-stratified SEIRD + hospital model in odin2,
#          extending Goodfellow et al. (2024).
#          Runs one epidemic per IMD deprivation decile.
#
# Inputs:
#   data/parameters/pi_matrix.csv
#   data/parameters/h_mu_by_age.csv
#   data/parameters/rural_age.csv
#   data/parameters/gamma_hd_hr_by_age.csv
#   data/parameters/contact_matrix_imd{1-10}.csv
#
# Output:
#   output/imd_hospital_gradient_odin.png

# ==============================================================

library(odin2)
library(dust2)
library(dplyr)
library(readr)
library(ggplot2)

age_levels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                "70 to 74","75+")

# Load parameters
pi_matrix <- read_csv("data/parameters/pi_matrix.csv", show_col_types = FALSE)
h_mu      <- read_csv("data/parameters/h_mu_by_age.csv", show_col_types = FALSE)
rural_age <- read_csv("data/parameters/rural_age.csv", show_col_types = FALSE) %>%
  mutate(Age = factor(Age, levels = age_levels))

gamma_hd_hr  <- read_csv("data/parameters/gamma_hd_hr_by_age.csv",
                         show_col_types = FALSE)
gamma_hd_vec <- gamma_hd_hr$gamma_hd
gamma_hr_vec <- gamma_hd_hr$gamma_hr
cat("gamma_hd range:", round(range(gamma_hd_vec), 4), "\n")
cat("gamma_hr range:", round(range(gamma_hr_vec), 4), "\n")

# odin2 model
# States (each length 17, ordered as per initial() declarations):
#   S, E, Ip, Ic, Is, HD, HR, D, R, Adm
# rec_c is passed as duration (days); used as 1/rec_c in ODEs
# rec_s is passed as rate (1/5 per day)
age_seird_hosp <- odin2::odin({
  
  # Force of infection
  inf_weighted[]  <- (Ip[i] + Ic[i] + xi * Is[i]) / proportion[i]
  dim(inf_weighted) <- 17
  
  weighted_contact[,] <- contact[i, j] * inf_weighted[j]
  dim(weighted_contact) <- c(17, 17)
  
  lambda[] <- susc * sum(weighted_contact[i,])
  dim(lambda) <- 17
  
  # ODEs
  deriv(S[])   <- -lambda[i] * S[i]
  deriv(E[])   <-  lambda[i] * S[i] - infec * E[i]
  deriv(Ip[])  <-  pi_a[i] * infec * E[i] - sympt * Ip[i]
  deriv(Ic[])  <-  sympt * Ip[i] - (1 / rec_c) * Ic[i]
  deriv(Is[])  <-  (1 - pi_a[i]) * infec * E[i] - rec_s * Is[i]
  deriv(HD[])  <-  (h_a[i] * mu_ca_h[i] / rec_c) * Ic[i] - gam_hd[i] * HD[i]
  deriv(HR[])  <-  (h_a[i] * (1 - mu_ca_h[i]) / rec_c) * Ic[i] - gam_hr[i] * HR[i]
  deriv(D[])   <-  gam_hd[i] * HD[i]
  deriv(R[])   <-  gam_hr[i] * HR[i] + ((1 - h_a[i]) / rec_c) * Ic[i] + rec_s * Is[i]
  deriv(Adm[]) <- (h_a[i] / rec_c) * Ic[i]
  
  dim(S)   <- 17
  dim(E)   <- 17
  dim(Ip)  <- 17
  dim(Ic)  <- 17
  dim(Is)  <- 17
  dim(HD)  <- 17
  dim(HR)  <- 17
  dim(D)   <- 17
  dim(R)   <- 17
  dim(Adm) <- 17
  
  # Initial conditions
  initial(S[])   <- S0[i]
  initial(E[])   <- 0
  initial(Ip[])  <- Ip0[i]
  initial(Ic[])  <- 0
  initial(Is[])  <- 0
  initial(HD[])  <- 0
  initial(HR[])  <- 0
  initial(D[])   <- 0
  initial(R[])   <- 0
  initial(Adm[]) <- 0
  
  # Array parameters
  S0[]         <- parameter()
  Ip0[]        <- parameter()
  proportion[] <- parameter()
  pi_a[]       <- parameter()
  h_a[]        <- parameter()
  mu_ca_h[]    <- parameter()
  gam_hd[]     <- parameter()
  gam_hr[]     <- parameter()
  contact[,]   <- parameter()
  
  dim(S0)         <- 17
  dim(Ip0)        <- 17
  dim(proportion) <- 17
  dim(pi_a)       <- 17
  dim(h_a)        <- 17
  dim(mu_ca_h)    <- 17
  dim(gam_hd)     <- 17
  dim(gam_hr)     <- 17
  dim(contact)    <- c(17, 17)
  
  # Scalar parameters
  susc  <- parameter(0.06)
  infec <- parameter(0.3333333)
  sympt <- parameter(0.4761905)
  rec_c <- parameter(2.9)
  rec_s <- parameter(0.2)
  xi    <- parameter(0.5)
})

# Run for one IMD decile
run_epidemic_odin <- function(imd_decile, urban = TRUE) {
  
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
  
  S0    <- proportion
  S0[8] <- S0[8] - 1e-3
  Ip0   <- c(rep(0, 7), 1e-3, rep(0, 9))
  
  sys <- dust2::dust_system_create(age_seird_hosp, list(
    S0         = S0,
    Ip0        = Ip0,
    proportion = proportion,
    pi_a       = pi_a,
    h_a        = h_a,
    mu_ca_h    = mu_ca_h,
    contact    = contact,
    gam_hd     = gamma_hd_vec,
    gam_hr     = gamma_hr_vec
  ))
  dust2::dust_system_set_state_initial(sys)
  
  times <- seq(0, 365, by = 1)
  out   <- dust2::dust_system_simulate(sys, times)
  
  # States ordered as initial() declarations:
  # S[1:17]=1:17, E=18:34, Ip=35:51, Ic=52:68, Is=69:85,
  # HD=86:102, HR=103:119, D=120:136, R=137:153, Adm=154:170
  n_age   <- 17
  s_idx   <- 1:n_age
  adm_idx <- (9 * n_age + 1):(10 * n_age)
  
  data.frame(
    day            = times,
    imd_decile     = imd_decile,
    cum_admissions = colSums(out[adm_idx, ]),
    S_remaining    = colSums(out[s_idx, ])
  ) %>%
    mutate(
      new_adm_per1000 = c(0, diff(cum_admissions)) * 1000,
      attack_rate     = 1 - S_remaining
    )
}

if (!exists("SKIP09_RUN") || !isTRUE(SKIP09_RUN)) {
  
  cat("Running odin2 model for 10 IMD deciles...\n")
  results_df <- bind_rows(lapply(1:10, function(d) {
    cat("  IMD decile", d, "\n")
    run_epidemic_odin(imd_decile = d)
  }))
  
  summary_stats <- results_df %>%
    group_by(imd_decile) %>%
    summarise(
      peak_adm_per1000  = max(new_adm_per1000, na.rm = TRUE),
      total_adm_per1000 = max(cum_admissions,  na.rm = TRUE) * 1000,
      attack_rate       = max(attack_rate,      na.rm = TRUE),
      .groups = "drop"
    )
  
  cat("\n--- Results by IMD decile (1 = most deprived) ---\n")
  print(summary_stats)
  
  decile_colours <- c(
    "#67001f","#b2182b","#d6604d","#f4a582","#fddbc7",
    "#d1e5f0","#92c5de","#4393c3","#2166ac","#053061"
  )
  
  p <- ggplot(results_df,
              aes(x = day, y = new_adm_per1000,
                  colour = factor(imd_decile, levels = 10:1),
                  group  = imd_decile)) +
    geom_line(linewidth = 0.9, alpha = 0.9) +
    scale_colour_manual(values = rev(decile_colours),
                        name = "IMD decile\n(1 = most\ndeprived)") +
    labs(
      title    = "Modelled daily hospital admissions by IMD deprivation decile",
      subtitle = "Age-stratified SEIRD + hospital model (odin2), England (urban)",
      x        = "Day of epidemic",
      y        = "New hospital admissions per 1,000 population",
      caption  = "Parameters: Goodfellow et al. (2024) + Knock et al. (2021)."
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_text(face = "bold", size = 12),
      plot.subtitle    = element_text(size = 9, colour = "#555555"),
      plot.caption     = element_text(size = 7.5, colour = "#888888", hjust = 0),
      panel.grid.minor = element_blank()
    )
  
  print(p)
  dir.create("output", showWarnings = FALSE)
  ggsave("output/imd_hospital_gradient_odin.png", p,
         width = 11, height = 6.5, dpi = 200)
  cat("Plot saved: output/imd_hospital_gradient_odin.png\n")
}

