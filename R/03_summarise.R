# =====================================================================
# 03_summarise.R
# Post-processing of the SimDesign output: reshape the wide condition-
# level results (columns named <measure>.<coefficient>) into a long
# performance table (one row per condition x coefficient), label
# rejection rates as Power vs Type I error, and write results/
# performance.rds / .csv used by 04_plots.R.
# =====================================================================

source("R/01_functions.R")
library(dplyr)
library(tidyr)

out_dir <- "results"

res <- readRDS(file.path(out_dir, "irr_simdesign.rds"))

perf <- res %>%
  as.data.frame() %>%
  select(nLevels, k, agree, nEvents, prob_type, contains(".")) %>%
  pivot_longer(cols = contains("."),
               names_to = c("measure", "coefficient"),
               names_sep = "\\.") %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(
    # Rejection of H0: coeff = 0 is Power where the population value is
    # non-zero and Type I error where it is exactly zero. Note that at
    # agree = 0 with skewed margins, AC1/BP have truth != 0, so their
    # rejection rate there is still power, not size.
    rejection_type = case_when(
      is.na(Rejection)                      ~ NA_character_,
      abs(truth) < .Machine$double.eps^0.5  ~ "TypeI",
      TRUE                                  ~ "Power"
    )
  )

saveRDS(perf, file.path(out_dir, "performance.rds"))
write.csv(perf, file.path(out_dir, "performance.csv"), row.names = FALSE)

cat("Wrote", nrow(perf), "condition x coefficient rows\n")

# Quick sanity views ---------------------------------------------------

# Convergence problems, if any
perf %>%
  filter(convergence < 1) %>%
  count(coefficient, nLevels, k, nEvents, wt = NULL, name = "n_conditions") %>%
  print(n = 20)

# Worst relative bias per coefficient (flags small-n problems)
perf %>%
  filter(agree > 0) %>%
  group_by(coefficient) %>%
  slice_max(abs(RelBias), n = 1) %>%
  select(coefficient, nLevels, k, agree, nEvents, prob_type,
         RelBias, RelBias_MCSE) %>%
  print(n = Inf)

# Type I error overview (only where truth == 0)
perf %>%
  filter(rejection_type == "TypeI") %>%
  group_by(coefficient) %>%
  summarise(mean_TypeI = mean(Rejection, na.rm = TRUE),
            max_TypeI  = max(Rejection, na.rm = TRUE),
            .groups = "drop") %>%
  print(n = Inf)
