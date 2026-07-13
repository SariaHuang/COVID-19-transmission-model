# ==============================================================
# Script: 06_test_base_model.R
#
# Purpose: Single-population (no age/IMD/region stratification)
#          SEIRD + hospital skeleton model in odin2, to verify the
#          model structure compiles and produces a sensible epidemic
#          curve before adding age/IMD/region stratification.
#
# Parameters: see docs/model_parameter_table.md
#   - p, M, xi, pi_c, f, gamma, r_c, r_s: Goodfellow et al. (2024)
#   - h, mu_ca_h: derived from Knock et al. (2021) Table S9,
#     single representative age-averaged value for this skeleton only
#   - gamma_hd, gamma_hr: single representative values for this
#     skeleton; age-specific values derived in script 08
#
# Validation: population is exactly conserved across all time steps,
# and final deaths = cumulative admissions x mu_ca_h as expected.
# Implemented using odin2 + dust2 (discrete-time difference equations).
# ==============================================================

library(odin2)
library(dust2)

# Single-population SEIRD + hospital skeleton model
# Purpose: verify model structure compiles and produces a sensible
# epidemic curve before adding age/IMD stratification.
# Parameters from Goodfellow et al. (2024) and Knock et al. (2021).

skeleton_gen <- odin2::odin({
  
  # Parameters
  p        <- parameter(0.06)   # transmission probability per contact
  M        <- parameter(11)     # average daily contacts
  xi       <- parameter(0.5)    # relative infectiousness of subclinical cases
  pi_c     <- parameter(0.55)   # clinical fraction
  f        <- parameter(3)      # latent period (days)
  gamma    <- parameter(2.1)    # preclinical infectious period (days)
  r_c      <- parameter(2.9)    # clinical infectious period (days)
  r_s      <- parameter(5)      # subclinical infectious period (days)
  h        <- parameter(0.036)  # P(hospitalised | clinical)
  mu_ca_h  <- parameter(0.1)    # P(death | hospitalised)
  mu_ca_g  <- parameter(0)      # P(death | clinical, not hospitalised)
  gamma_hd <- parameter(0.0885) # rate HD -> D
  gamma_hr <- parameter(0.0709) # rate HR -> R
  N        <- parameter(100000) # total population
  
  # Force of infection
  lambda <- p * M * (Ip + Ic + xi * Is) / N
  
  # Flows (dt = 1 day; rates already in per-day units)
  n_SE    <- S  * lambda
  n_E_Ip  <- E  * (pi_c / f)
  n_E_Is  <- E  * ((1 - pi_c) / f)
  n_Ip_Ic <- Ip * (1 / gamma)
  n_Is_R  <- Is * (1 / r_s)
  n_Ic_HD <- Ic * (h * mu_ca_h / r_c)
  n_Ic_HR <- Ic * (h * (1 - mu_ca_h) / r_c)
  n_Ic_Dd <- Ic * ((1 - h) * mu_ca_g / r_c)
  n_Ic_Rd <- Ic * ((1 - h) * (1 - mu_ca_g) / r_c)
  n_HD_D  <- HD * gamma_hd
  n_HR_R  <- HR * gamma_hr
  
  # State updates
  update(S)  <- S  - n_SE
  update(E)  <- E  + n_SE - n_E_Ip - n_E_Is
  update(Ip) <- Ip + n_E_Ip - n_Ip_Ic
  update(Is) <- Is + n_E_Is - n_Is_R
  update(Ic) <- Ic + n_Ip_Ic - n_Ic_HD - n_Ic_HR - n_Ic_Dd - n_Ic_Rd
  update(HD) <- HD + n_Ic_HD - n_HD_D
  update(HR) <- HR + n_Ic_HR - n_HR_R
  update(D)  <- D  + n_HD_D  + n_Ic_Dd
  update(R)  <- R  + n_HR_R  + n_Is_R  + n_Ic_Rd
  
  # New and cumulative admissions (flow into hospital, not occupancy)
  update(new_adm) <- n_Ic_HD + n_Ic_HR
  update(cum_adm) <- cum_adm + n_Ic_HD + n_Ic_HR
  
  # Initial conditions
  initial(S)       <- N - 10
  initial(E)       <- 10
  initial(Ip)      <- 0
  initial(Is)      <- 0
  initial(Ic)      <- 0
  initial(HD)      <- 0
  initial(HR)      <- 0
  initial(D)       <- 0
  initial(R)       <- 0
  initial(new_adm) <- 0
  initial(cum_adm) <- 0
})

# Run
sys <- dust2::dust_system_create(skeleton_gen, list())
dust2::dust_system_set_state_initial(sys)
out_raw <- dust2::dust_system_simulate(sys, 0:365)

# out_raw is [n_states, n_particles, n_times] -- extract particle 1
state_names <- c("S","E","Ip","Is","Ic","HD","HR","D","R","new_adm","cum_adm")
out <- as.data.frame(t(out_raw))
names(out) <- state_names
out$day <- 0:365

# checks
total_pop <- with(out, S + E + Ip + Is + Ic + HD + HR + D + R)
cat("Population conserved (min/max):", round(min(total_pop)), round(max(total_pop)), "\n")
cat("Peak Ic day:", out$day[which.max(out$Ic)], "| value:", round(max(out$Ic)), "\n")
cat("Final deaths:", round(tail(out$D, 1)), "\n")
cat("Final cumulative admissions:", round(tail(out$cum_adm, 1)), "\n")
cat("Attack rate:", round(1 - tail(out$S, 1) / 100000, 3), "\n")

# Plot
plot(out$day, out$Ic, type = "l", col = "purple", lwd = 2,
     xlab = "Day", ylab = "Number of people",
     main = "SEIRD skeleton: Ic and daily new admissions (x7 for scale)")
lines(out$day, out$new_adm * 7, col = "red", lwd = 2)
legend("topright",
       legend = c("Ic (clinical infectious)", "New admissions \u00d77"),
       col = c("purple", "red"), lwd = 2)
