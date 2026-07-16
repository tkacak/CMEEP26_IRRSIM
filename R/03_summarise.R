# =====================================================================
# 03_summarise.R
# Performance measures with Monte Carlo standard errors (MCSE),
# following Siepe et al. (2024, Psychological Methods) and
# Morris, White & Crowther (2019, Statistics in Medicine):
#
#   Bias        = mean(est) - truth          MCSE = sqrt(var(est) / S)
#   RelBias     = Bias / truth               MCSE = MCSE(Bias) / |truth|
#   EmpSE       = sd(est)                    MCSE = EmpSE / sqrt(2 (S - 1))
#   RMSE        = sqrt(mean((est - truth)^2))
#               MCSE = sqrt(var((est - truth)^2) / S) / (2 * RMSE)
#   Rejection   = mean(p < .05)              MCSE = sqrt(R (1 - R) / S)
#               (= Power when truth != 0, Type I error when truth == 0)
#
# S = number of converged (non-missing) replications; the convergence
# rate itself is reported per Siepe et al.'s reporting checklist.
# =====================================================================

source("R/01_functions.R")
library(dplyr)

out_dir  <- "results"
cond_dir <- file.path(out_dir, "conditions")
alpha_level <- 0.05

summarise_condition <- function(df) {
  df %>%
    group_by(condition, nLevels, k, agree, nEvents, prob_type, coefficient) %>%
    summarise(
      n_reps      = n(),
      S           = sum(!is.na(estimate)),
      convergence = S / n_reps,

      truth = true_value(first(coefficient), first(agree),
                         first(nLevels), first(prob_type)),

      mean_est  = mean(estimate, na.rm = TRUE),
      Bias      = mean_est - truth,
      Bias_MCSE = sqrt(var(estimate, na.rm = TRUE) / S),

      RelBias      = if (first(truth) != 0) Bias / truth else NA_real_,
      RelBias_MCSE = if (first(truth) != 0) Bias_MCSE / abs(truth) else NA_real_,

      EmpSE      = sd(estimate, na.rm = TRUE),
      EmpSE_MCSE = EmpSE / sqrt(2 * (S - 1)),

      MSE       = mean((estimate - truth)^2, na.rm = TRUE),
      RMSE      = sqrt(MSE),
      RMSE_MCSE = sqrt(var((estimate - truth)^2, na.rm = TRUE) / S) / (2 * RMSE),

      S_test         = sum(!is.na(p_value)),
      Rejection      = if (first(S_test) > 0)
                         mean(p_value < alpha_level, na.rm = TRUE) else NA_real_,
      Rejection_MCSE = if (first(S_test) > 0)
                         sqrt(Rejection * (1 - Rejection) / S_test) else NA_real_,

      .groups = "drop"
    ) %>%
    mutate(
      # Rejection of H0: coeff = 0 is Power where the population value
      # is non-zero and Type I error where it is exactly zero. Note that
      # at agree = 0 with skewed margins, AC1/BP/PA have truth != 0, so
      # their rejection rate there is still power, not size.
      rejection_type = case_when(
        is.na(Rejection)          ~ NA_character_,
        abs(truth) < .Machine$double.eps^0.5 ~ "TypeI",
        TRUE                      ~ "Power"
      )
    )
}

files <- list.files(cond_dir, pattern = "^cond_.*\\.rds$", full.names = TRUE)
stopifnot(length(files) > 0)

perf <- bind_rows(lapply(files, function(f) summarise_condition(readRDS(f))))

saveRDS(perf, file.path(out_dir, "performance.rds"))
write.csv(perf, file.path(out_dir, "performance.csv"), row.names = FALSE)

cat("Summarised", length(files), "conditions ->",
    nrow(perf), "coefficient x condition rows\n")

# Quick sanity views ---------------------------------------------------

# Worst relative bias per coefficient (should flag small-n problems)
perf %>%
  filter(agree > 0) %>%
  group_by(coefficient) %>%
  slice_max(abs(RelBias), n = 1) %>%
  select(coefficient, nLevels, k, agree, nEvents, prob_type,
         RelBias, RelBias_MCSE) %>%
  print(n = Inf)

# Type I error overview (only meaningful in the agree = 0 conditions)
perf %>%
  filter(rejection_type == "TypeI") %>%
  group_by(coefficient) %>%
  summarise(mean_TypeI = mean(Rejection, na.rm = TRUE),
            max_TypeI  = max(Rejection, na.rm = TRUE),
            .groups = "drop") %>%
  print(n = Inf)
