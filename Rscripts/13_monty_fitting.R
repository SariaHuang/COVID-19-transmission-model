# ==============================================================
# Script: 13_monty_fitting.R
#
# Purpose: Fit the age x IMD decile SEIRD + hospital model to
#          observed NHS weekly hospital admissions using monty MCMC.
#
# Approach: WRAPPER approach (not native odin2/dust2/monty).
#   - Reuses the validated odin2 model from script 09 (age_seird_hosp)
#   - monty_model_function wraps a hand-written log-likelihood
#
# Parameters fitted (per IMD decile):
#   log_beta  -- log transmission probability
#   log_size  -- log NegBin overdispersion parameter
#
# I0_frac fixed at 1e-4:
#   Not identifiable separately from beta in a single-wave
#   deterministic ODE model. Fixed at 1e-4 (~1 case per 10,000).
#
# pop_size: blended urban/rural population (consistent with
#   proportion vector and contact matrices in script 09).
#
# Data window:
#   obs_data from 2020-07-27, filtered to autumn/winter wave 1.
#   Wave end manually set to 2021-04-05; auto-detection returned
#   2021-05-10 which extends beyond the genuine inter-wave trough.
# ==============================================================

library(odin2)
library(dust2)
library(monty)
library(dplyr)
library(readr)
library(zoo)
library(ggplot2)

# Load model and parameters from scripts 08 and 09
if (!exists("age_seird_hosp")) {
  SKIP09_RUN <- TRUE
  source("Rscripts/08_refine_gamma_hd_hr.R")
  source("Rscripts/09_age_imd_stratified_odin.R")
  rm(SKIP09_RUN)
}
stopifnot(exists("age_seird_hosp"))
stopifnot(exists("get_blended_inputs"))
cat("odin2 model loaded: age_seird_hosp\n")

# State indices in dust2 output (170 states: 10 compartments x 17 ages)
adm_idx <- 154:170
s_idx   <- 1:17

# Observed weekly admissions
obs_data <- read_csv("data/processed/observed_weekly_admissions.csv",
                     show_col_types = FALSE)

cat("obs_data:", nrow(obs_data), "rows,",
    n_distinct(obs_data$epiweek), "weeks, from",
    format(min(obs_data$epiweek)), "to", format(max(obs_data$epiweek)), "\n")

# Wave 1 end: manually set to 2021-04-05
model_day0 <- min(obs_data$epiweek)

obs_d1_full <- obs_data %>%
  filter(lad_imd_decile == 1) %>%
  arrange(epiweek) %>%
  mutate(smoothed = zoo::rollmean(obs_admissions, k = 3,
                                  fill = NA, align = "center"))

peak_idx       <- which.max(obs_d1_full$smoothed)
first_wave_end <- as.Date("2021-04-05")

cat("\nWave 1 window:\n")
cat("  Model day 0 :", format(model_day0), "\n")
cat("  Peak date   :", format(obs_d1_full$epiweek[peak_idx]), "\n")
cat("  Trough date :", format(first_wave_end), "(manually set)\n\n")

obs_data_wave1 <- obs_data %>%
  filter(epiweek >= model_day0, epiweek <= first_wave_end) %>%
  arrange(epiweek)

n_wave1_weeks <- n_distinct(obs_data_wave1$epiweek)
cat("  Wave 1 window:", n_wave1_weeks, "weeks\n\n")

if (n_wave1_weeks < 8) {
  stop("Wave 1 window is only ", n_wave1_weeks,
       " weeks -- check smoothing or set first_wave_end manually.")
}

# Run odin2 model for one decile using blended population
run_epidemic_fit <- function(imd_decile, beta, I0_frac) {
  
  contact <- as.matrix(read.csv(
    paste0("data/parameters/contact_matrix_imd", imd_decile, ".csv"),
    header = FALSE
  ))
  
  pi_a    <- pi_matrix[[paste0("imd_", imd_decile)]]
  h_a     <- h_mu$h_a
  mu_ca_h <- h_mu$mu_ca_h
  
  blended    <- get_blended_inputs(imd_decile)
  proportion <- blended$proportion
  
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
  
  data.frame(
    day            = 0:365,
    cum_admissions = colSums(out[adm_idx, ]),
    S_remaining    = colSums(out[s_idx,   ])
  ) %>%
    mutate(attack_rate = 1 - S_remaining)
}

# Negative binomial log-likelihood
# pop_size: blended urban/rural population
run_log_likelihood <- function(beta, I0_frac, size,
                               imd_decile, obs_weekly,
                               n_weeks, pop_size) {
  tryCatch({
    out       <- run_epidemic_fit(imd_decile, beta, I0_frac)
    daily_adm <- c(0, diff(out$cum_admissions)) * pop_size
    
    pred_weekly <- sapply(seq_len(n_weeks), function(w) {
      day_start <- 7 * (w - 1) + 1
      day_end   <- min(7 * w, nrow(out))
      if (day_start > nrow(out)) return(0)
      sum(daily_adm[day_start:day_end])
    })
    
    pred_weekly <- pmax(pred_weekly, 0.01)
    sum(dnbinom(obs_weekly, mu = pred_weekly, size = size, log = TRUE))
    
  }, error = function(e) -Inf)
}

# Pre-MCMC structural checks
cat("--- Pre-MCMC structural checks ---\n")

obs_d1_w1 <- obs_data_wave1 %>%
  filter(lad_imd_decile == 1) %>%
  arrange(epiweek)

stopifnot("No wave 1 data for decile 1" = nrow(obs_d1_w1) > 0)

# Blended population for decile 1
pop_size_d1 <- get_blended_inputs(1)$pop_size

check <- run_epidemic_fit(imd_decile = 1, beta = 0.06, I0_frac = 1e-4)
stopifnot(all(is.finite(check$cum_admissions)))
stopifnot(all(diff(check$cum_admissions) >= -1e-9))

cat("  Attack rate at starting values:",
    round(tail(check$attack_rate, 1), 4), "\n")
cat("  Cum admissions (proportion):  ",
    round(tail(check$cum_admissions, 1), 5), "\n")

ll_test <- run_log_likelihood(
  beta = 0.06, I0_frac = 1e-4, size = 10,
  imd_decile = 1,
  obs_weekly = obs_d1_w1$obs_admissions,
  n_weeks    = nrow(obs_d1_w1),
  pop_size   = pop_size_d1
)
cat("  Log-likelihood at starting values:", round(ll_test, 2), "\n")
if (!is.finite(ll_test)) stop("Likelihood not finite at starting values.")
cat("  All pre-MCMC checks passed.\n\n")

# Pred vs obs diagnostic plot at starting values
daily_adm_d1 <- c(0, diff(check$cum_admissions)) * pop_size_d1
n_w1         <- nrow(obs_d1_w1)

pred_start <- sapply(seq_len(n_w1), function(w) {
  day_start <- 7 * (w - 1) + 1
  day_end   <- min(7 * w, nrow(check))
  if (day_start > nrow(check)) return(0)
  sum(daily_adm_d1[day_start:day_end])
})

plot(seq_len(n_w1), obs_d1_w1$obs_admissions,
     type = "l", col = "black",
     xlab = "Week from 2020-07-27",
     ylab = "Weekly admissions",
     main = "Decile 1: observed (black) vs predicted at starting values (red)")
lines(seq_len(n_w1), pred_start, col = "red")
legend("topright", legend = c("Observed","Predicted"),
       col = c("black","red"), lty = 1)

# monty MCMC fitting function
# pop_size: blended urban/rural population per decile
fit_decile <- function(imd_decile, n_samples = 2000, n_chains = 3) {
  
  cat("=== Fitting IMD decile", imd_decile, "===\n")
  
  obs_d <- obs_data_wave1 %>%
    filter(lad_imd_decile == imd_decile) %>%
    arrange(epiweek)
  
  if (nrow(obs_d) == 0) {
    warning("No wave 1 data for decile ", imd_decile, " -- skipping")
    return(NULL)
  }
  
  obs_weekly <- obs_d$obs_admissions
  n_weeks    <- nrow(obs_d)
  
  pop_size <- get_blended_inputs(imd_decile)$pop_size
  
  ll_fn <- function(log_beta, log_size) {
    run_log_likelihood(
      beta       = exp(log_beta),
      I0_frac    = 1e-4,
      size       = exp(log_size),
      imd_decile = imd_decile,
      obs_weekly = obs_weekly,
      n_weeks    = n_weeks,
      pop_size   = pop_size
    )
  }
  
  prior <- monty::monty_dsl({
    log_beta ~ Normal(-2.81, 0.5)
    log_size ~ Normal(1, 1)
  })
  
  likelihood <- monty::monty_model_function(ll_fn)
  posterior  <- likelihood + prior
  sampler    <- monty::monty_sampler_random_walk(diag(2) * c(0.02, 0.04))
  
  samples <- monty::monty_sample(
    posterior, sampler, n_samples,
    initial  = c(log(0.06), log(10)),
    n_chains = n_chains
  )
  
  cat("  Done:", n_samples, "samples x", n_chains, "chains\n")
  return(samples)
}

# Run all 10 deciles
dir.create("output/fitting", recursive = TRUE, showWarnings = FALSE)

cat("Fitting all 10 IMD deciles...\n")
all_fits <- lapply(1:10, function(d) {
  fit <- fit_decile(imd_decile = d, n_samples = 2000, n_chains = 3)
  if (!is.null(fit)) {
    saveRDS(fit, paste0("output/fitting/fitted_samples_imd", d, ".rds"))
  }
  fit
})
cat("All deciles done. Results saved to output/fitting/\n")

# Posterior summary
burnin <- 500
cat("\n--- Posterior summary, all deciles (burn-in 500 discarded) ---\n")
cat(sprintf("%-10s  %-30s  %-30s\n",
            "Decile", "beta median (C1/C2/C3)", "size median (C1/C2/C3)"))

for (d in 1:10) {
  fit <- all_fits[[d]]
  if (is.null(fit)) next
  beta_med <- sapply(1:3, function(ch)
    round(median(exp(fit$pars[1, (burnin+1):2000, ch])), 5))
  size_med <- sapply(1:3, function(ch)
    round(median(exp(fit$pars[2, (burnin+1):2000, ch])), 3))
  cat(sprintf("Decile %2d:  beta = %s   size = %s\n",
              d,
              paste(beta_med, collapse = " / "),
              paste(size_med, collapse = " / ")))
}

# Validation plot: pred vs obs, all 10 deciles
cat("\nDrawing validation plot...\n")
dir.create("output/validation", recursive = TRUE, showWarnings = FALSE)

par(mfrow = c(2, 5), mar = c(4, 4, 2, 1))

for (d in 1:10) {
  fit <- all_fits[[d]]
  if (is.null(fit)) next
  
  beta_post <- median(sapply(1:3, function(ch)
    median(exp(fit$pars[1, (burnin+1):2000, ch]))))
  
  obs_d <- obs_data_wave1 %>%
    filter(lad_imd_decile == d) %>%
    arrange(epiweek)
  
  pop_size <- get_blended_inputs(d)$pop_size
  
  out       <- run_epidemic_fit(imd_decile = d, beta = beta_post, I0_frac = 1e-4)
  daily_adm <- c(0, diff(out$cum_admissions)) * pop_size
  n_w       <- nrow(obs_d)
  
  pred <- sapply(seq_len(n_w), function(w) {
    d1 <- 7 * (w - 1) + 1
    d2 <- min(7 * w, nrow(out))
    if (d1 > nrow(out)) return(0)
    sum(daily_adm[d1:d2])
  })
  
  plot(seq_len(n_w), obs_d$obs_admissions,
       type = "l", col = "black",
       main = paste("Decile", d),
       xlab = "Week", ylab = "Admissions",
       ylim = range(c(obs_d$obs_admissions, pred), na.rm = TRUE))
  lines(seq_len(n_w), pred, col = "red")
  
  if (d == 1) {
    legend("topleft", legend = c("Observed","Predicted"),
           col = c("black","red"), lty = 1, cex = 0.8)
  }
}

dev.copy(png, "output/validation/pred_vs_obs_all_deciles.png",
         width = 1600, height = 700, res = 120)
dev.off()
cat("Validation plot saved: output/validation/pred_vs_obs_all_deciles.png\n")
