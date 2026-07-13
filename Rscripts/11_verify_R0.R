# ==============================================================
# Script: 11_verify_R0.R
#
# Purpose: Verify that contact matrices and clinical fraction parameters
#          produce plausible R0 values, using the next-generation matrix
#          method from Goodfellow et al. (2024):
#
#            N[i,j] = p * M[i,j] * (pi_j*(gamma + r_c) + xi*(1-pi_j)*r_s)
#            R0     = dominant eigenvalue of N
#
#          Expected range from Goodfellow et al. (2024):
#          2.09 (rural, decile 7) to 2.71 (urban, most deprived).
#
# Note: This is a parameter validation script. The fitted beta and
#       R0 values used in the main analysis are derived in script 13
#       (MCMC fitting) and script 17 (school closure scenario).
#
# Inputs:
#   data/parameters/clin_frac.csv
#   data/parameters/G.csv
#   data/parameters/rural_age.csv
# ==============================================================


library(dplyr)
library(readr)

age_levels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                "70 to 74","75+")

clin_frac <- read_csv("data/parameters/clin_frac.csv",
                      show_col_types = FALSE) %>%
  mutate(Ageband = factor(Ageband, levels = age_levels))

G <- as.matrix(read_csv("data/parameters/G.csv", show_col_types = FALSE))

rural_age <- read_csv("data/parameters/rural_age.csv",
                      show_col_types = FALSE) %>%
  mutate(Age = factor(Age, levels = age_levels))

# Parameters (Goodfellow et al. 2024)
p     <- 0.06
gamma <- 2.1
r_c   <- 2.9
r_s   <- 5
xi    <- 0.5

compute_contact_matrix <- function(imd_decile, urban = TRUE) {
  pop_vec <- rural_age %>%
    filter(IMD == imd_decile,
           rural == if (urban) "Urban" else "Rural") %>%
    arrange(Age) %>%
    pull(Population)
  M <- matrix(nrow = 17, ncol = 17)
  for (i in 1:17) for (j in 1:17) M[i,j] <- G[i,j] * pop_vec[j] / sum(pop_vec)
  M
}

compute_R0 <- function(imd_decile, urban = TRUE) {
  M    <- compute_contact_matrix(imd_decile, urban)
  pi_a <- clin_frac %>%
    filter(IMD == imd_decile) %>%
    arrange(Ageband) %>%
    pull(clin_frac)
  N <- matrix(nrow = 17, ncol = 17)
  for (i in 1:17) for (j in 1:17) {
    N[i,j] <- p * M[i,j] * (pi_a[j]*(gamma + r_c) + xi*(1-pi_a[j])*r_s)
  }
  max(Mod(eigen(N, only.values = TRUE)$values))
}

r0_results <- data.frame(
  imd_decile = 1:10,
  R0_urban   = sapply(1:10, compute_R0, urban = TRUE),
  R0_rural   = sapply(1:10, compute_R0, urban = FALSE)
)

print(r0_results)
cat("\nGoodfellow et al. (2024) reports R0: 2.09 (rural, decile 7) to 2.71 (urban, decile 1)\n")
cat("Urban decile 1:", round(r0_results$R0_urban[1], 3),
    "| Rural decile 7:", round(r0_results$R0_rural[7], 3), "\n")

