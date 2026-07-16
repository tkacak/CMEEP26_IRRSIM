# =====================================================================
# 03_summarise.R
# Post-processing of the SimDesign output: reshape the wide condition-
# level results (columns named <measure>.<coefficient>) into a long
# performance table (one row per condition x coefficient), label
# rejection rates as Power vs Type I error, attach the comparison
# factors (coefficient family, weighting scheme), and write results/
# performance.rds / .csv used by 04_plots.R.
# Run with the R/ folder as the working directory.
# =====================================================================

source("01_functions.R")
library(dplyr)
library(tidyr)

out_dir <- "results"

res <- readRDS(file.path(out_dir, "irr_simdesign.rds"))

perf <- res %>%
  as.data.frame() %>%
  select(nLevels, k, rho, nEvents, prob_type, contains(".")) %>%
  pivot_longer(cols = contains("."),
               names_to = c("measure", "coefficient"),
               names_sep = "\\.") %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(
    # Rejection of H0: coeff = 0 is Power where the population value is
    # non-zero and Type I error where it is exactly zero. Note that at
    # rho = 0 with skewed margins, AC1/AC2 have truth != 0, so their
    # rejection rate there is still power, not size.
    rejection_type = case_when(
      is.na(Rejection)                      ~ NA_character_,
      abs(truth) < .Machine$double.eps^0.5  ~ "TypeI",
      TRUE                                  ~ "Power"
    ),
    # Weighted variants duplicate their unweighted counterparts when
    # nLevels == 2 (see coef_applicable in 01_functions.R); figures and
    # reports keep applicable == TRUE rows only
    applicable = coef_applicable(coefficient, nLevels),
    coef_name  = coef_display(coefficient),
    # Comparison factors for the weighted-vs-unweighted contrast
    family     = coef_family(coefficient),
    weighting  = factor(coef_weighting(coefficient),
                        levels = c("unweighted", "linear", "quadratic"))
  )

saveRDS(perf, file.path(out_dir, "performance.rds"))
write.csv(perf, file.path(out_dir, "performance.csv"), row.names = FALSE)

cat("Wrote", nrow(perf), "condition x coefficient rows\n")

# Quick sanity views ---------------------------------------------------

# Convergence problems, if any
perf %>%
  filter(convergence < 1) %>%
  count(coefficient, nLevels, k, nEvents, name = "n_conditions") %>%
  print(n = 20)

# Worst relative bias per coefficient (flags small-n problems)
perf %>%
  filter(rho > 0, applicable) %>%
  group_by(coefficient) %>%
  slice_max(abs(RelBias), n = 1) %>%
  select(coefficient, nLevels, k, rho, nEvents, prob_type,
         RelBias, RelBias_MCSE) %>%
  print(n = Inf)

# Weighted vs unweighted at a glance: mean |RelBias| and mean power by
# family x weighting (collapsed over conditions with rho > 0)
perf %>%
  filter(rho > 0, applicable) %>%
  group_by(family, weighting) %>%
  summarise(mean_abs_RelBias = mean(abs(RelBias), na.rm = TRUE),
            mean_EmpSE       = mean(EmpSE, na.rm = TRUE),
            mean_Power       = mean(Rejection[rejection_type == "Power"],
                                    na.rm = TRUE),
            .groups = "drop") %>%
  print(n = Inf)

# Type I error overview (only where truth == 0)
perf %>%
  filter(rejection_type == "TypeI") %>%
  group_by(coefficient) %>%
  summarise(mean_TypeI = mean(Rejection, na.rm = TRUE),
            max_TypeI  = max(Rejection, na.rm = TRUE),
            .groups = "drop") %>%
  print(n = Inf)
