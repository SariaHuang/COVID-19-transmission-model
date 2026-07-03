# ==============================================================
# Script: 13_monty_fitting.R
#
# Purpose: Fit the age x IMD decile SEIRD + hospital model to
#          observed NHS weekly hospital admissions using monty MCMC.
#
# Approach: WRAPPER approach (not native odin2/dust2/monty).
#   - Reuses the validated odin model from script 09 (age_seird_hosp)
#   - monty_model_function wraps a hand-written log-likelihood
#
# Parameters fitted (per IMD decile):
#   log_beta  -- log transmission probability
#   log_size  -- log NegBin overdispersion parameter
#
# I0_frac fixed at 1e-4:
#   Not identifiable separately from beta in a single-wave
#   deterministic ODE model (beta and I0_frac are confounded along
#   a likelihood ridge). Fixed at 1e-4 (~1 case per 10,000).
#   Reported as a methodological note in dissertation.
#
# Data window (Path A):
#   obs_data from 2020-07-27, filtered to autumn/winter wave 1
#   (first trough after first peak, auto-detected from decile 1).
#
# Prerequisite: add guard to script 09 section 4:
#   if (!exists("SKIP09_RUN") || !isTRUE(SKIP09_RUN)) { ... }
# ==============================================================

library(odin)
library(monty)
library(dplyr)
library(readr)
library(zoo)
library(ggplot2)

# ------------------------------------------------------------
# 1. Source model definition from scripts 08 and 09
# ------------------------------------------------------------
SKIP09_RUN <- TRUE
source("Rscripts/08_refine_gamma_hd_hr.R")
source("Rscripts/09_age_imd_stratified_odin.R")
rm(SKIP09_RUN)

stopifnot(exists("age_seird_hosp"))
cat("odin model loaded: age_seird_hosp\n")


# ------------------------------------------------------------
# 2. Load observed data
# ------------------------------------------------------------
obs_data <- read_csv("data/processed/observed_weekly_admissions.csv",
                     show_col_types = FALSE)

cat("obs_data:", nrow(obs_data), "rows,",
    n_distinct(obs_data$epiweek), "weeks, from",
    format(min(obs_data$epiweek)), "to", format(max(obs_data$epiweek)), "\n")

# ------------------------------------------------------------
# 3. Define model day 0 and auto-detect wave 1 end
#
#    model day 0 = 2020-07-27 (earliest obs_data date)
#    model day k = date (model_day0 + k days)
#    epiweek 1 = model days 1-7, epiweek 2 = days 8-14, etc.
# ------------------------------------------------------------
model_day0 <- min(obs_data$epiweek)

obs_d1_full <- obs_data %>%
  filter(lad_imd_decile == 1) %>%
  arrange(epiweek) %>%
  mutate(smoothed = zoo::rollmean(obs_admissions, k = 3,
                                  fill = NA, align = "center"))

peak_idx       <- which.max(obs_d1_full$smoothed)
after_peak     <- obs_d1_full[peak_idx:nrow(obs_d1_full), ]
trough_idx_rel <- which.min(after_peak$smoothed)
first_wave_end <- after_peak$epiweek[trough_idx_rel]

cat("\nWave 1 window (auto-detected from decile 1):\n")
cat("  Model day 0 :", format(model_day0), "\n")
cat("  Peak date   :", format(obs_d1_full$epiweek[peak_idx]), "\n")
cat("  Trough date :", format(first_wave_end), "\n")
cat("  VERIFY: does the trough date look like the genuine inter-wave low?\n")
cat("  If not, set first_wave_end manually before continuing.\n\n")

obs_data_wave1 <- obs_data %>%
  filter(epiweek >= model_day0, epiweek <= first_wave_end) %>%
  arrange(epiweek)

n_wave1_weeks <- n_distinct(obs_data_wave1$epiweek)
cat("  Wave 1 window:", n_wave1_weeks, "weeks\n\n")

if (n_wave1_weeks < 8) {
  stop("Wave 1 window is only ", n_wave1_weeks,
       " weeks -- check smoothing or set first_wave_end manually.")
}

# ------------------------------------------------------------
# 4. run_epidemic_fit
#    Wrapper around the validated odin model from script 09.
#    Accepts beta and I0_frac overrides; all other inputs
#    (contact matrices, pi_a, h_a etc.) loaded by script 09.
# ------------------------------------------------------------
run_epidemic_fit <- function(imd_decile, beta, I0_frac, urban = TRUE) {
  
  setting <- if (urban) "Urban" else "Rural"
  
  contact <- as.matrix(read.csv(
    paste0("data/parameters/contact_matrix_imd", imd_decile, ".csv"),
    header = FALSE
  ))
  
  pi_a    <- pi_matrix[[paste0("imd_", imd_decile)]]
  h_a     <- h_mu$h_a
  mu_ca_h <- h_mu$mu_ca_h
  
  proportion <- rural_age %>%
    filter(IMD == imd_decile, rural == setting) %>%
    arrange(Age) %>%
    pull(Proportion)
  
  S0    <- proportion
  S0[8] <- S0[8] - I0_frac
  Ip0   <- c(rep(0, 7), I0_frac, rep(0, 9))
  
  mod <- age_seird_hosp$new(
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
  )
  
  times <- seq(0, 365, by = 1)
  out   <- as.data.frame(mod$run(times))
  
  adm_cols <- grep("^Adm\\[", names(out))
  S_cols   <- grep("^S\\[",   names(out))
  
  data.frame(
    day            = out$t,
    cum_admissions = rowSums(out[, adm_cols]),
    S_remaining    = rowSums(out[, S_cols])
  ) %>%
    mutate(attack_rate = 1 - S_remaining)
}

# ------------------------------------------------------------
# 5. Log-likelihood function
# ------------------------------------------------------------
run_log_likelihood <- function(beta, I0_frac, size,
                               imd_decile, obs_weekly,
                               n_weeks, pop_size) {
  tryCatch({
    
    out <- run_epidemic_fit(imd_decile = imd_decile,
                            beta = beta, I0_frac = I0_frac)
    
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

# ------------------------------------------------------------
# 6. Pre-MCMC structural checks
#    If any check fails, do NOT proceed to MCMC.
# ------------------------------------------------------------
cat("--- Pre-MCMC structural checks ---\n")

obs_d1_w1 <- obs_data_wave1 %>%
  filter(lad_imd_decile == 1) %>%
  arrange(epiweek)

stopifnot("No wave 1 data for decile 1" = nrow(obs_d1_w1) > 0)

pop_size_d1 <- rural_age %>%
  filter(IMD == 1) %>%
  summarise(pop = sum(Population)) %>%
  pull(pop)

check <- run_epidemic_fit(imd_decile = 1, beta = 0.06, I0_frac = 1e-4)
stopifnot(all(is.finite(check$cum_admissions)))
stopifnot(all(diff(check$cum_admissions) >= -1e-9))

cat("  Attack rate at starting values:",
    round(tail(check$attack_rate, 1), 4), "\n")
cat("  Cum admissions (proportion):  ",
    round(tail(check$cum_admissions, 1), 5), "\n")

ll_test <- run_log_likelihood(
  beta       = 0.06,
  I0_frac    = 1e-4,
  size       = 10,
  imd_decile = 1,
  obs_weekly = obs_d1_w1$obs_admissions,
  n_weeks    = nrow(obs_d1_w1),
  pop_size   = pop_size_d1
)
cat("  Log-likelihood at starting values:", round(ll_test, 2), "\n")

if (!is.finite(ll_test)) {
  stop("Likelihood not finite at starting values.\n",
       "Check: wave 1 window dates, pop_size scale, model peak timing.")
}
cat("  All pre-MCMC checks passed.\n\n")

# Pred vs obs diagnostic plot
daily_adm_d1 <- c(0, diff(check$cum_admissions)) * pop_size_d1
n_w1         <- nrow(obs_d1_w1)

pred_start <- sapply(seq_len(n_w1), function(w) {
  day_start <- 7 * (w - 1) + 1
  day_end   <- min(7 * w, nrow(check))
  if (day_start > nrow(check)) return(0)
  sum(daily_adm_d1[day_start:day_end])
})

dev.new()
plot(seq_len(n_w1), obs_d1_w1$obs_admissions,
     type = "l", col = "black",
     xlab = "Week from model day 0 (2020-07-27)",
     ylab = "Weekly admissions",
     main = "Decile 1, wave 1: observed (black) vs predicted at starting values (red)")
lines(seq_len(n_w1), pred_start, col = "red")
legend("topright",
       legend = c("observed", "predicted"),
       col = c("black", "red"), lty = 1)
cat("Inspect pred vs obs plot before proceeding to MCMC.\n\n")

# ------------------------------------------------------------
# 7. fit_decile: monty MCMC
#
#    Fits log_beta and log_size only.
#    I0_frac fixed at 1e-4 (not identifiable from beta in a
#    single-wave deterministic model -- see header note).
# ------------------------------------------------------------
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
  
  pop_size <- rural_age %>%
    filter(IMD == imd_decile) %>%
    summarise(pop = sum(Population)) %>%
    pull(pop)
  
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
  
  vcv     <- diag(2) * c(0.02, 0.04)
  sampler <- monty::monty_sampler_random_walk(vcv)
  
  samples <- monty::monty_sample(
    posterior, sampler, n_samples,
    initial  = c(log(0.06), log(10)),
    n_chains = n_chains
  )
  
  cat("  Done:", n_samples, "samples x", n_chains, "chains\n")
  return(samples)
}

# ------------------------------------------------------------
# 8. Formal run: decile 1, 2000 samples x 3 chains
#
#    Convergence criteria before scaling to all 10 deciles:
#    (a) beta: all 3 chains stable and overlapping after burn-in
#    (b) size: not trending toward 0, chains overlapping
#    (c) size median > 1 (model is genuinely fitting, not
#        absorbing residuals via overdispersion)
# ------------------------------------------------------------
dir.create("output/fitting", recursive = TRUE, showWarnings = FALSE)

cat("Fitting decile 1: 2000 samples x 3 chains...\n")
fit_d1 <- fit_decile(imd_decile = 1, n_samples = 2000, n_chains = 3)

if (!is.null(fit_d1)) {
  saveRDS(fit_d1, "output/fitting/fitted_samples_imd1.rds")
  cat("Saved: output/fitting/fitted_samples_imd1.rds\n")
  
  # Trace plots: all 3 chains for each parameter
  chain_cols <- c("black", "steelblue", "firebrick")
  dev.new()
  par(mfrow = c(2, 3))
  
  for (ch in 1:3) {
    pars <- fit_d1$pars[,,ch]
    plot(exp(pars[1,]), type = "l", col = chain_cols[ch],
         main  = paste("beta — chain", ch),
         ylab  = "beta", xlab = "Index",
         ylim  = range(sapply(1:3, function(c)
           range(exp(fit_d1$pars[1,,c])))))
  }
  for (ch in 1:3) {
    pars <- fit_d1$pars[,,ch]
    plot(exp(pars[2,]), type = "l", col = chain_cols[ch],
         main  = paste("size — chain", ch),
         ylab  = "size", xlab = "Index",
         ylim  = range(sapply(1:3, function(c)
           range(exp(fit_d1$pars[2,,c])))))
  }
  
  # Posterior summary across all chains (discard first 500 as burn-in)
  burnin <- 500
  cat("\n--- Posterior summary, decile 1 (chains 1-3, after burn-in) ---\n")
  for (ch in 1:3) {
    pars <- fit_d1$pars[, (burnin + 1):2000, ch]
    cat(sprintf("Chain %d:  beta median = %.5f   size median = %.3f\n",
                ch,
                median(exp(pars[1,])),
                median(exp(pars[2,]))))
  }
  cat("\nIf all 3 chains agree on beta and size medians,\n")
  cat("proceed to section 9 (full 10-decile run).\n")
}

# ------------------------------------------------------------
# 9. Full 10-decile run (uncomment after section 8 passes)
# ------------------------------------------------------------
# cat("\nFitting all 10 IMD deciles...\n")
# all_fits <- lapply(1:10, function(d) {
#   fit <- fit_decile(imd_decile = d, n_samples = 2000, n_chains = 3)
#   if (!is.null(fit)) {
#     saveRDS(fit, paste0("output/fitting/fitted_samples_imd", d, ".rds"))
#   }
#   fit
# })
# cat("All deciles done. Results saved to output/fitting/\n")