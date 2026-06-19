# ==============================================================
# Script: 08_age_imd_stratified_odin.R
#
# Purpose: Extend Goodfellow et al. (2024)'s age-stratified SEIRD
#          model to add hospital compartments (HD, HR), and run
#          separately for each of the 10 IMD deprivation deciles.
#
# Structure follows Goodfellow's SEIRD_model.R, with two changes:
#   (1) Hospital compartments HD and HR inserted between Ic and
#       final outcomes, using h_a and mu_ca_h from Knock (2021).
#   (2) ODE written without with() to avoid scoping issues.
#
# Inputs:
#   data/parameters/pi_matrix.csv
#   data/parameters/h_mu_by_age.csv
#   data/parameters/contact_matrix_imd1.csv ... imd10.csv
#   data/parameters/rural_age.csv
# ==============================================================

library(odin)
library(dplyr)
library(readr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load parameters (same as script 08)
# ------------------------------------------------------------
age_levels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                "70 to 74","75+")

pi_matrix <- read_csv("data/parameters/pi_matrix.csv", show_col_types = FALSE)
h_mu      <- read_csv("data/parameters/h_mu_by_age.csv", show_col_types = FALSE)
rural_age <- read_csv("data/parameters/rural_age.csv", show_col_types = FALSE) %>%
  mutate(Age = factor(Age, levels = age_levels))

gamma_hd <- 1 / 11.3
gamma_hr <- 1 / 14.1

# ------------------------------------------------------------
# 2. odin model definition
#
#    State order doesn't matter in odin (unlike deSolve) -- each
#    compartment is its own named array, indexed by age 1:17.
#    Force of infection uses sum(contact[i,] * inf_weighted[])
#    as odin's equivalent of matrix %*% vector.
# ------------------------------------------------------------
age_seird_hosp <- odin::odin({
  
  n_age <- 17
  
  # ---- Force of infection ----
  inf_weighted[] <- (Ip[i] + Ic[i] + xi * Is[i]) / proportion[i]
  dim(inf_weighted) <- n_age
  
  # sum() cannot take a compound expression directly -- compute the
  # elementwise product matrix first, then sum each row of it
  weighted_contact[, ] <- contact[i, j] * inf_weighted[j]
  dim(weighted_contact) <- c(n_age, n_age)
  
  lambda[] <- susc * sum(weighted_contact[i, ])
  dim(lambda) <- n_age
  
  # ---- Derivatives (one line per compartment, vectorised over age) ----
  deriv(S[])  <- -lambda[i] * S[i]
  deriv(E[])  <-  lambda[i] * S[i] - infec * E[i]
  deriv(Ip[]) <-  pi_a[i] * infec * E[i] - sympt * Ip[i]
  deriv(Ic[]) <-  sympt * Ip[i] - (1 / rec_c) * Ic[i]
  deriv(Is[]) <-  (1 - pi_a[i]) * infec * E[i] - rec_s * Is[i]
  deriv(HD[]) <-  (h_a[i] * mu_ca_h[i] / rec_c) * Ic[i] - gam_hd * HD[i]
  deriv(HR[]) <-  (h_a[i] * (1 - mu_ca_h[i]) / rec_c) * Ic[i] - gam_hr * HR[i]
  deriv(D[])  <-  gam_hd * HD[i]
  deriv(R[])  <-  gam_hr * HR[i] + ((1 - h_a[i]) / rec_c) * Ic[i] + rec_s * Is[i]
  deriv(Adm[]) <- (h_a[i] / rec_c) * Ic[i]
  
  dim(S)   <- n_age
  dim(E)   <- n_age
  dim(Ip)  <- n_age
  dim(Ic)  <- n_age
  dim(Is)  <- n_age
  dim(HD)  <- n_age
  dim(HR)  <- n_age
  dim(D)   <- n_age
  dim(R)   <- n_age
  dim(Adm) <- n_age
  
  # ---- Initial conditions ----
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
  
  # ---- User-supplied data (passed in from R at $new()) ----
  S0[]  <- user()
  dim(S0) <- n_age
  Ip0[] <- user()
  dim(Ip0) <- n_age
  proportion[] <- user()
  dim(proportion) <- n_age
  pi_a[] <- user()
  dim(pi_a) <- n_age
  h_a[] <- user()
  dim(h_a) <- n_age
  mu_ca_h[] <- user()
  dim(mu_ca_h) <- n_age
  contact[, ] <- user()
  dim(contact) <- c(n_age, n_age)
  
  # ---- Scalar parameters (same values as Goodfellow / script 08) ----
  susc   <- user(0.06)    # transmission probability per contact
  infec  <- user(0.3333333) # 1/3, rate out of E
  sympt  <- user(0.4761905) # 1/2.1, rate Ip -> Ic
  rec_c  <- user(2.9)      # NOTE: used as denominator (1/rec_c), so pass duration not rate
  rec_s  <- user(0.2)      # 1/5, rate out of Is
  xi     <- user(0.5)
  gam_hd <- user()
  gam_hr <- user()
})

# ------------------------------------------------------------
# 3. Run for one IMD decile
# ------------------------------------------------------------
run_epidemic_odin <- function(imd_decile, urban = TRUE) {
  
  setting <- if (urban) "Urban" else "Rural"
  
  contact <- as.matrix(read.csv(
    paste0("data/parameters/contact_matrix_imd", imd_decile, ".csv"),
    header = FALSE
  ))
  
  pi_a <- pi_matrix[[paste0("imd_", imd_decile)]]
  h_a     <- h_mu$h_a
  mu_ca_h <- h_mu$mu_ca_h
  
  proportion <- rural_age %>%
    filter(IMD == imd_decile, rural == setting) %>%
    arrange(Age) %>%
    pull(Proportion)
  
  S0  <- proportion
  S0[8] <- S0[8] - 1e-3
  Ip0 <- c(rep(0, 7), 1e-3, rep(0, 9))
  
  mod <- age_seird_hosp$new(
    S0         = S0,
    Ip0        = Ip0,
    proportion = proportion,
    pi_a       = pi_a,
    h_a        = h_a,
    mu_ca_h    = mu_ca_h,
    contact    = contact,
    gam_hd     = gamma_hd,
    gam_hr     = gamma_hr
  )
  
  times <- seq(0, 365, by = 1)
  out   <- as.data.frame(mod$run(times))
  
  # odin names array outputs like Adm[1], Adm[2], ..., S[1], S[2]...
  adm_cols <- grep("^Adm\\[", names(out))
  S_cols   <- grep("^S\\[",   names(out))
  
  result <- data.frame(
    day            = out$t,
    imd_decile     = imd_decile,
    cum_admissions = rowSums(out[, adm_cols]),
    S_remaining    = rowSums(out[, S_cols])
  ) %>%
    mutate(
      new_adm_per1000 = c(0, diff(cum_admissions)) * 1000,
      attack_rate     = 1 - S_remaining
    )
  
  return(result)
}

# ------------------------------------------------------------
# 4. Run for all 10 IMD deciles
# ------------------------------------------------------------
cat("Running odin model for 10 IMD deciles...\n")
all_results <- lapply(1:10, function(d) {
  cat("  IMD decile", d, "\n")
  run_epidemic_odin(imd_decile = d, urban = TRUE)
})
results_df <- bind_rows(all_results)

# ------------------------------------------------------------
# 5. Sanity checks (should match script 08's deSolve results closely)
# ------------------------------------------------------------
summary_stats <- results_df %>%
  group_by(imd_decile) %>%
  summarise(
    peak_adm_per1000  = max(new_adm_per1000, na.rm = TRUE),
    total_adm_per1000 = max(cum_admissions, na.rm = TRUE) * 1000,
    attack_rate       = max(attack_rate, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n--- Results by IMD decile (1 = most deprived) ---\n")
print(summary_stats)
cat("\nCompare these numbers against script 08's deSolve output --\n")
cat("they should be very close (odin and deSolve solve the same ODEs).\n")

# ------------------------------------------------------------
# 6. Plot (same as script 08)
# ------------------------------------------------------------
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
    subtitle = "Age-stratified SEIRD + hospital model (odin), England (urban)",
    x        = "Day of epidemic",
    y        = "New hospital admissions per 1,000 population",
    caption  = "Parameters: Goodfellow et al. (2024) + Knock et al. (2021)."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "#555555"),
    plot.caption  = element_text(size = 7.5, colour = "#888888", hjust = 0),
    panel.grid.minor = element_blank()
  )

print(p)
dir.create("output", showWarnings = FALSE)
ggsave("output/imd_hospital_gradient_odin.png", p, width = 11, height = 6.5, dpi = 200)
cat("Plot saved to output/imd_hospital_gradient_odin.png\n")
