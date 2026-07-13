# ==============================================================
# Script: 11_verify_R0.R
#
# Purpose: Independently verify that the model's parameters
#          (contact matrices, clinical fraction) produce a
#          plausible R0, using the next-generation-matrix method
#          from Goodfellow et al. (2024):
#
#            N_ij = p * M_ij * [pi_j*(gamma + r_c) + xi*(1-pi_j)*r_s]
#            R0   = dominant eigenvalue of N
#
# Expected (from Goodfellow et al. 2024, Results section):
#   R0 ranges from 2.09 (rural, decile 7) to 2.71 (urban, most
#   deprived decile)
#
# Inputs: data/parameters/clin_frac.csv, G.csv, rural_age.csv
# ==============================================================

library(dplyr)
library(readr)

age_levels <- c("Under 1","1 to 4","5 to 9","10 to 14","15 to 19",
                "20 to 24","25 to 29","30 to 34","35 to 39","40 to 44",
                "45 to 49","50 to 54","55 to 59","60 to 64","65 to 69",
                "70 to 74","75+")

clin_frac <- read_csv("data/parameters/clin_frac.csv", show_col_types = FALSE) %>%
  mutate(Ageband = factor(Ageband, levels = age_levels))

G <- as.matrix(read_csv("data/parameters/G.csv", show_col_types = FALSE))

rural_age <- read_csv("data/parameters/rural_age.csv", show_col_types = FALSE) %>%
  mutate(Age = factor(Age, levels = age_levels))

compute_contact_matrix <- function(imd_decile, urban = TRUE) {
  setting <- if (urban) "Urban" else "Rural"
  pop_vec <- rural_age %>%
    filter(IMD == imd_decile, rural == setting) %>%
    arrange(Age) %>%
    pull(Population)
  M <- matrix(nrow = 17, ncol = 17)
  for (i in 1:17) {
    for (j in 1:17) {
      M[i, j] <- G[i, j] * pop_vec[j] / sum(pop_vec)
    }
  }
  return(M)
}

# Parameters (durations in days, matching docs/model_parameter_table.md)
p     <- 0.06
gamma <- 2.1
r_c   <- 2.9
r_s   <- 5
xi    <- 0.5

compute_R0 <- function(imd_decile, urban = TRUE) {
  M <- compute_contact_matrix(imd_decile, urban)
  pi_a <- clin_frac %>%
    filter(IMD == imd_decile) %>%
    arrange(Ageband) %>%
    pull(clin_frac)
  
  N <- matrix(nrow = 17, ncol = 17)
  for (i in 1:17) {
    for (j in 1:17) {
      d_j <- pi_a[j] * (gamma + r_c) + xi * (1 - pi_a[j]) * r_s
      N[i, j] <- p * M[i, j] * d_j
    }
  }
  ev <- eigen(N, only.values = TRUE)$values
  max(Mod(ev))
}

cat("R0 by IMD decile (urban):\n")
r0_results <- data.frame(
  imd_decile = 1:10,
  R0_urban   = sapply(1:10, compute_R0, urban = TRUE)
)
print(r0_results)

cat("\nGoodfellow et al. (2024) reports R0 range: 2.09 (rural, decile 7) to 2.71 (urban, most deprived)\n")
cat("Our urban decile 1 R0:", round(r0_results$R0_urban[1], 3),
    "-- compare to paper's reported 2.71\n")
cat("\nIf these numbers are close, the contact matrices and clinical fraction\n")
cat("data are loading and combining correctly -- independent confirmation\n")
cat("that the model's transmission parameters are working as intended.\n")
