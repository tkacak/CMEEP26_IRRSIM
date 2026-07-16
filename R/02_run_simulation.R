# =====================================================================
# 02_run_simulation.R
# SimDesign (Chalmers & Adkins, 2020) implementation of the full
# factorial Monte Carlo simulation comparing weighted vs unweighted
# agreement coefficients (kappa family and Gwet AC family) under the
# Underlying Variable Approach (Muthen, 1984).
#
# Design (ADEMP; see README.md):
#   nLevels   : 3, 4, 5 ordinal categories
#   k         : 2 raters (both rate all events)
#   rho       : 0.30 ... 0.90 by 0.10 latent correlation
#               (+ 0.00 null condition for Type I error)
#   nEvents   : 20, 30, 40, 50
#   prob_type : uniform, skew (category prevalence via fixed thresholds)
#   -> 192 conditions
#
# Run this script with the R/ folder as the working directory.
# Performance measures with analytic Monte Carlo standard errors
# (Siepe et al., 2024; Morris, White & Crowther, 2019) are computed in
# Summarise() against numerically computed population values. Raw
# replication-level results are stored per condition via
# save_results = TRUE; an interrupted run resumes automatically from
# SimDesign's tempfile when this script is re-run.
# =====================================================================

library(SimDesign)
source("01_functions.R")

# ------------------------- configuration -----------------------------

n_reps      <- 100   # test runs; set to 1000 for the final run (MCSE targets)
alpha_level <- 0.05
out_dir     <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- design grid ------------------------------

Design <- createDesign(
  nLevels   = 3:5,
  k         = 2,
  rho       = c(0, seq(0.30, 0.90, by = 0.10)),  # 0 = null (Type I error)
  nEvents   = c(20, 30, 40, 50),
  prob_type = c("uniform", "skew")
)

# ----------------------- generate / analyse ---------------------------

Generate <- function(condition, fixed_objects) {
  uva_generate(nEvents   = condition$nEvents,
               rho       = condition$rho,
               nLevels   = condition$nLevels,
               prob_type = condition$prob_type)
}

Analyse <- function(condition, dat, fixed_objects) {
  # returns est.<coef> and p.<coef>; failures become NA inside
  # fit_all_coefficients() so no data redraws occur (see 01_functions.R)
  fit_all_coefficients(dat, condition$nLevels)
}

# ----------------------------- summarise -------------------------------
# Per coefficient: Bias, RelBias, EmpSE, RMSE, rejection rate of
# H0: coeff = 0, and convergence rate - each with its MCSE:
#   Bias      MCSE = sqrt(var(est)/S)
#   RelBias   MCSE = MCSE(Bias)/|truth|
#   EmpSE     MCSE = EmpSE/sqrt(2(S-1))
#   RMSE      MCSE = sqrt(var((est-truth)^2)/S) / (2*RMSE)  [delta method]
#   Rejection MCSE = sqrt(R(1-R)/S)
# S = converged (non-missing) replications. truth = the coefficient's
# own population value under the UVA process (see 01_functions.R).

Summarise <- function(condition, results, fixed_objects) {

  alpha <- fixed_objects$alpha
  out <- c()

  for (cf in fixed_objects$coefs) {

    est <- results[[paste0("est.", cf)]]
    pv  <- results[[paste0("p.",  cf)]]

    R <- length(est)
    S <- sum(!is.na(est))

    truth <- true_value(cf, condition$rho,
                        condition$nLevels, condition$prob_type)

    mean_est  <- mean(est, na.rm = TRUE)
    Bias      <- mean_est - truth
    Bias_MCSE <- sqrt(stats::var(est, na.rm = TRUE) / S)

    RelBias      <- if (truth != 0) Bias / truth else NA_real_
    RelBias_MCSE <- if (truth != 0) Bias_MCSE / abs(truth) else NA_real_

    EmpSE      <- stats::sd(est, na.rm = TRUE)
    EmpSE_MCSE <- EmpSE / sqrt(2 * (S - 1))

    sq_err    <- (est - truth)^2
    RMSE      <- sqrt(mean(sq_err, na.rm = TRUE))
    RMSE_MCSE <- sqrt(stats::var(sq_err, na.rm = TRUE) / S) / (2 * RMSE)

    S_test <- sum(!is.na(pv))
    Rejection      <- if (S_test > 0) mean(pv < alpha, na.rm = TRUE) else NA_real_
    Rejection_MCSE <- if (S_test > 0) sqrt(Rejection * (1 - Rejection) / S_test) else NA_real_

    v <- c(truth = truth, convergence = S / R,
           Bias = Bias, Bias_MCSE = Bias_MCSE,
           RelBias = RelBias, RelBias_MCSE = RelBias_MCSE,
           EmpSE = EmpSE, EmpSE_MCSE = EmpSE_MCSE,
           RMSE = RMSE, RMSE_MCSE = RMSE_MCSE,
           Rejection = Rejection, Rejection_MCSE = Rejection_MCSE)

    out <- c(out, stats::setNames(v, paste(names(v), cf, sep = ".")))
  }

  out
}

# ------------------------------- run -----------------------------------

res <- runSimulation(
  design        = Design,
  replications  = n_reps,
  generate      = Generate,
  analyse       = Analyse,
  summarise     = Summarise,
  fixed_objects = list(coefs = coef_labels(), alpha = alpha_level),
  packages      = c("irrCAC"),
  seed          = genSeeds(design = Design, iseed = 22),
  parallel      = TRUE,
  filename      = file.path(out_dir, "irr_simdesign"),
  save_results  = TRUE,
  save_details  = list(save_results_dirname = file.path(out_dir, "raw-results")),
  progress      = TRUE,
  control       = list(allow_na = TRUE, allow_nan = TRUE)
)

# ------------------- reproducibility bookkeeping ----------------------

print(res)
print(SimExtract(res, what = "errors"))

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))

cat("Done. Summarised results: results/irr_simdesign.rds\n",
    "Raw per-condition results: results/raw-results/\n", sep = "")
