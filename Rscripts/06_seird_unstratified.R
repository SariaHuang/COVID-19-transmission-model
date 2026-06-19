# ==============================================================
# Script: 06_test_base_model.R
#
# Purpose: Single-population (no age/IMD/region stratification)
#          SEIR-HD skeleton model in odin, to verify the model
#          structure compiles and produces a sensible epidemic
#          curve before adding age/IMD/region stratification.
#
# Parameters: see docs/model_parameter_table.md
#   - p, M, xi, pi_c, f, gamma, r_c, r_s: Goodfellow et al. (2024)
#   - h, mu_ca_h: derived from Knock et al. (2021) Table S9,
#     single representative value for this skeleton (not yet
#     age/IMD-stratified)
#   - gamma_hd, gamma_hr: provisional placeholders, pending proper
#     probability-weighted aggregation from Knock et al. (2021)
#     Table S2/S6/S8
#
# This logic was validated in base R before translating to odin:
# population is exactly conserved, and final deaths/admissions
# match pi_c * h * mu_ca_h and pi_c * h respectively.
# ==============================================================

library(odin)

skeleton_gen <- odin::odin({
  
  # ---- Parameters ----
  p        <- user(0.06)   # transmission probability per contact
  M        <- user(11)     # average daily contacts (single population)
  xi       <- user(0.5)    # relative infectiousness of subclinical cases
  pi_c     <- user(0.55)   # clinical fraction (single representative value)
  f        <- user(3)      # latent period duration (days)
  gamma    <- user(2.1)    # preclinical infectious period duration (days)
  r_c      <- user(2.9)    # clinical infectious period duration (days)
  r_s      <- user(5)      # subclinical infectious period duration (days)
  h        <- user(0.036)  # P(hospitalised | clinical), placeholder
  mu_ca_h  <- user(0.1)    # P(death | hospitalised), placeholder
  mu_ca_g  <- user(0)      # P(death | clinical, not hospitalised)
  gamma_hd <- user(0.0885) # rate HD -> D, provisional placeholder
  gamma_hr <- user(0.0709) # rate HR -> R, provisional placeholder
  N        <- user(100000) # total population
  dt       <- user(1)      # time step (days)
  
  # ---- Force of infection ----
  lambda <- p * M * (Ip + Ic + xi * Is) / N
  
  # ---- Flows (number of individuals moving per time step) ----
  n_SE    <- S * lambda * dt
  n_E_Ip  <- E * (pi_c / f) * dt
  n_E_Is  <- E * ((1 - pi_c) / f) * dt
  n_Ip_Ic <- Ip * (1 / gamma) * dt
  n_Is_R  <- Is * (1 / r_s) * dt
  n_Ic_HD <- Ic * (h * mu_ca_h / r_c) * dt
  n_Ic_HR <- Ic * (h * (1 - mu_ca_h) / r_c) * dt
  n_Ic_Dd <- Ic * ((1 - h) * mu_ca_g / r_c) * dt
  n_Ic_Rd <- Ic * ((1 - h) * (1 - mu_ca_g) / r_c) * dt
  n_HD_D  <- HD * gamma_hd * dt
  n_HR_R  <- HR * gamma_hr * dt
  
  # ---- Updates ----
  update(S)  <- S - n_SE
  update(E)  <- E + n_SE - n_E_Ip - n_E_Is
  update(Ip) <- Ip + n_E_Ip - n_Ip_Ic
  update(Is) <- Is + n_E_Is - n_Is_R
  update(Ic) <- Ic + n_Ip_Ic - n_Ic_HD - n_Ic_HR - n_Ic_Dd - n_Ic_Rd
  update(HD) <- HD + n_Ic_HD - n_HD_D
  update(HR) <- HR + n_Ic_HR - n_HR_R
  update(D)  <- D + n_HD_D + n_Ic_Dd
  update(R)  <- R + n_HR_R + n_Is_R + n_Ic_Rd
  
  # New hospital admissions this step (flow, not stock) -- this is
  # what should be compared against real weekly hosp_admissions data,
  # NOT the standing HD+HR stock (which is current occupancy)
  update(new_admissions) <- n_Ic_HD + n_Ic_HR
  update(cum_admissions) <- cum_admissions + n_Ic_HD + n_Ic_HR
  
  # ---- Initial conditions ----
  initial(S)  <- N - 10
  initial(E)  <- 10
  initial(Ip) <- 0
  initial(Is) <- 0
  initial(Ic) <- 0
  initial(HD) <- 0
  initial(HR) <- 0
  initial(D)  <- 0
  initial(R)  <- 0
  initial(new_admissions) <- 0
  initial(cum_admissions) <- 0
})

# ------------------------------------------------------------
# Run the model
# ------------------------------------------------------------
mod <- skeleton_gen$new()
n_days <- 365
out <- as.data.frame(mod$run(0:n_days))

# ------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------
total_pop <- out$S + out$E + out$Ip + out$Is + out$Ic + out$HD + out$HR + out$D + out$R
cat("Population conserved (min/max total):", min(total_pop), max(total_pop), "\n")
cat("Peak Ic day:", out$step[which.max(out$Ic)], "| value:", round(max(out$Ic)), "\n")
cat("Final deaths:", round(tail(out$D, 1)), "\n")
cat("Final cumulative admissions:", round(tail(out$cum_admissions, 1)), "\n")
cat("Attack rate:", round(1 - tail(out$S, 1) / 100000, 3), "\n")

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plot(out$step, out$Ic, type = "l", col = "purple", lwd = 2,
     xlab = "Day", ylab = "Number of people",
     main = "SEIRD unstratified baseline: clinical infectious (Ic) and new admissions")
lines(out$step, out$new_admissions * 7, col = "red", lwd = 2)
legend("topright", legend = c("Ic (clinical infectious)", "New admissions x7"),
       col = c("purple", "red"), lwd = 2)
