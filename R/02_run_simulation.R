# =====================================================================
# 02_run_simulation.R
# Full factorial Monte Carlo simulation.
#
# Design (ADEMP; see README.md):
#   nLevels   : 2, 3, 4, 5
#   k         : 2, 3, 4, 5 raters (all raters rate all events)
#   agree     : 0.30 ... 0.90 by 0.10 (+ 0.00 null condition for Type I)
#   nEvents   : 20, 30, 40, 50
#   prob_type : uniform, skew
#   n_reps    : 1000 replications per condition
#
# Each condition's raw replication-level results are written to
# results/conditions/cond_XXXX.rds so an interrupted run can resume.
# =====================================================================

source("R/01_functions.R")
library(future.apply)

# ------------------------- configuration -----------------------------

n_reps       <- 1000
include_null <- TRUE     # add agree = 0 to estimate Type I error
out_dir      <- "results"
cond_dir     <- file.path(out_dir, "conditions")
n_workers    <- max(1L, parallel::detectCores() - 1L)

dir.create(cond_dir, recursive = TRUE, showWarnings = FALSE)

# --------------------------- design grid ------------------------------

agree_vec <- seq(0.30, 0.90, by = 0.10)
if (include_null) agree_vec <- c(0, agree_vec)

param_grid <- expand.grid(
  nLevels   = 2:5,
  k         = 2:5,
  agree     = agree_vec,
  nEvents   = c(20, 30, 40, 50),
  prob_type = c("uniform", "skew"),
  stringsAsFactors = FALSE
)
param_grid$condition <- seq_len(nrow(param_grid))

cat("Conditions:", nrow(param_grid),
    "| replications per condition:", n_reps,
    "| total replications:", nrow(param_grid) * n_reps, "\n")

# --------------------------- run --------------------------------------

plan(multisession, workers = n_workers)
set.seed(20260716)  # future.seed = TRUE derives one L'Ecuyer stream per condition

invisible(future_lapply(seq_len(nrow(param_grid)), function(i) {

  res_file <- file.path(cond_dir, sprintf("cond_%04d.rds", i))
  if (file.exists(res_file)) return(NULL)  # resume support

  cond <- param_grid[i, ]

  reps <- lapply(seq_len(n_reps), function(r) {
    out <- run_one_rep(cond$nLevels, cond$k, cond$agree,
                       cond$nEvents, cond$prob_type)
    out$rep <- r
    out
  })

  out <- cbind(cond, do.call(rbind, reps), row.names = NULL)
  saveRDS(out, res_file)
  NULL
}, future.seed = TRUE))

plan(sequential)

# ------------------- reproducibility bookkeeping ----------------------

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo.txt"))
saveRDS(param_grid, file.path(out_dir, "param_grid.rds"))

cat("Done. Per-condition files in", cond_dir, "\n",
    "Monitor progress with: length(list.files('", cond_dir, "'))\n", sep = "")
