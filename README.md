# CMEEP26_IRRSIM — Comparing Inter-Rater Agreement Coefficients

Bu klasörde Arş. Gör. Tugay KAÇAK ve Doç. Dr. Abdullah Faruk KILIÇ'ın CMEEP26'da sunduğu puanlayıcılar arası uyum katsayılarının karşılaştırıldığı çalışmanın kodları bulunmaktadır.

A Monte Carlo simulation comparing the performance of chance-corrected
inter-rater agreement coefficients, implemented in the
[`SimDesign`](https://cran.r-project.org/package=SimDesign) framework
(Chalmers & Adkins, 2020). Ordinal ratings are generated with the
Underlying Variable Approach (UVA; Muthén, 1984, as used by Almehrizi,
2025). Reporting follows the simulation template of Siepe et al. (2024)
with performance measures and Monte Carlo standard errors (MCSE) from
Morris, White & Crowther (2019).

## Design (ADEMP)

**Aims.** Test whether *weighted* agreement coefficients outperform their
unweighted ("default") counterparts on ordinal rating scales: Cohen's kappa
vs quadratic-weighted kappa (QWK), and Gwet's AC1 vs AC2 with linear and
quadratic weights — compared on bias, relative bias, empirical SE, and
power (plus Type I error) across small rating designs.

**Data-generating mechanism.** Underlying Variable Approach (UVA; Muthén,
1984). Two standard-normal latent variables with correlation ρ (the latent
"true agreement" of the design) are discretised into L ordinal categories
using a common set of thresholds, so disagreements naturally concentrate on
adjacent categories — the scenario weighted coefficients were designed for.
Thresholds are the normal quantiles of the target category prevalences:
uniform (*p_i = 1/L*) or skewed (*p_i ∝ i*). Full factorial design:

| Factor      | Levels                                    |
|-------------|-------------------------------------------|
| `nLevels`   | 3, 4, 5                                    |
| `k` (raters)| 2 (both raters rate all events)            |
| `rho`       | 0.30–0.90 by 0.10, plus 0.00 (null)        |
| `nEvents`   | 20, 30, 40, 50                             |
| `prob_type` | uniform; skew (*p_i ∝ i*, via thresholds)  |

8 × 3 × 1 × 4 × 2 = **192 conditions**, `n_reps = 1000` replications each.

**Deviation from Almehrizi (2025).** In that flow the thresholds are drawn
at random per data set. Random thresholds make every replication estimate a
*different* population value, which rules out bias/RMSE/power against a
known truth (accordingly, that study reports only descriptive summaries:
mean estimates, exceedance percentages, correlations). Here the thresholds
are **fixed design factors**, so every coefficient has a known population
value in every condition and the full Siepe et al. performance framework
applies.

**Estimands / true values.** With ρ and the thresholds fixed, the
population L × L joint table of the two raters' categories is fully
determined (bivariate-normal rectangle probabilities, computed by 1-D
Gaussian quadrature in `uva_joint()` — base R only). Each coefficient's
population value follows by plugging this table and the margins into its
defining formula (`true_value()` in `R/01_functions.R`). Under UVA the
weighted and unweighted coefficients have genuinely **different estimands**
(unlike mixture-type DGPs where they coincide), so each estimator is
evaluated against its own true value; note also that the kappa-family true
values are substantially smaller than ρ itself — ρ lives on the latent
correlation scale, not the chance-corrected agreement scale.

**Methods (coefficients compared).** Two families × weighting schemes,
crossed so the weighted-vs-unweighted contrast is made *within* family:

| Label        | Family | Weighting  | Implementation |
|--------------|--------|------------|----------------|
| `Kappa`      | Kappa  | unweighted | `irrCAC::conger.kappa.raw` (= Cohen's kappa at k = 2) |
| `Kappa_quad` | Kappa  | quadratic  | `irrCAC::conger.kappa.raw(weights = "quadratic")` (= QWK) |
| `AC1`        | Gwet   | unweighted | `irrCAC::gwet.ac1.raw` |
| `AC2_linear` | Gwet   | linear     | `irrCAC::gwet.ac1.raw(weights = "linear")` |
| `AC2_quad`   | Gwet   | quadratic  | `irrCAC::gwet.ac1.raw(weights = "quadratic")` |

All `irrCAC` estimators receive `categ.labels = 1:nLevels` so unobserved
categories in small samples cannot shrink the category space, and they all
report a subject-level-variance test of H0: coefficient = 0 (used for
power). `Metrics::ScoreQuadraticWeightedKappa` from the v1 script was
dropped: it estimates the same quantity as `Kappa_quad` but without a test.

Estimation failures are caught *inside* `Analyse()` and returned as `NA`
rather than thrown, so SimDesign does **not** redraw the data on failure —
redraws would condition results on estimability and bias the small-sample
conditions. The NA rate is reported as the convergence rate.

**Condition-dependent coefficient set.** With `k = 2` Conger's kappa
reduces exactly to **Cohen's kappa** in every condition, and with
`nLevels ≥ 3` the weighted variants are genuinely distinct estimators, so
the full set is applicable everywhere. `coef_applicable()` still guards
the degenerate case should the grid change (`nLevels = 2`: all weighted
variants are mathematically identical to their unweighted counterparts
and would be dropped from figures).

**Performance measures** (per condition × coefficient, with MCSE; S =
converged replications):

- Bias = mean(est) − truth; MCSE = √(var(est)/S)
- Relative bias = Bias/truth; MCSE = MCSE(Bias)/|truth|
- Empirical SE; MCSE = EmpSE/√(2(S−1))
- RMSE; MCSE via delta method = MCSE(MSE)/(2·RMSE)
- Rejection rate of H0: coefficient = 0 at α = .05; MCSE = √(R(1−R)/S).
  This is **power** where the truth ≠ 0 and **Type I error** where the truth
  = 0 (the `rho = 0` conditions — but note under skew the AC family has
  non-zero truth even at `rho = 0`, which `rejection_type` accounts for).
- Convergence rate (proportion of non-missing estimates), per the Siepe et
  al. reporting checklist.

With S = 1000, the worst-case MCSE of a rejection rate is
√(.5 × .5 / 1000) ≈ .016, and the MCSE of bias is EmpSE/√1000.

## Repository layout

```
R/01_functions.R      Analytic true values + coefficient estimators
R/02_run_simulation.R SimDesign Generate/Analyse/Summarise + runSimulation
R/03_summarise.R      Reshape SimDesign output to a long performance table
R/04_plots.R          Relative bias / bias / power / Type I figures
archive/IRR_sim_v1.R  Original draft script (kept for provenance)
results/, figs/       Generated output (git-ignored)
```

## Running

Scripts are written to run with `R/` as the working directory:

```r
install.packages(c("SimDesign", "irrCAC", "dplyr", "tidyr", "ggplot2"))
setwd("R")
source("02_run_simulation.R")  # writes results/irr_simdesign.rds (+ raw-results/)
source("03_summarise.R")       # writes results/performance.rds / .csv
source("04_plots.R")           # writes figs/*.png
```

`n_reps` is set to 100 for test runs; use 1000 for the final run (MCSE
targets in the Performance measures section assume S = 1000).

`runSimulation()` handles parallelisation (`parallel = TRUE`), per-condition
seeds (`seed = genSeeds(Design, iseed = 22)`, so the run is fully
reproducible), and progress/ETA reporting. If the run is interrupted,
re-sourcing the script resumes automatically from SimDesign's tempfile.
Replication-level raw results are stored in `results/raw-results/`
(`save_results = TRUE`), the condition-level summary in
`results/irr_simdesign.rds`, and error/warning bookkeeping is available via
`SimExtract(res, "errors")`. `sessionInfo()` is saved alongside the
results, per the Siepe et al. checklist.

## References

- Almehrizi, R. S. (2025). Simulation study on weighted inter-rater
  agreement coefficients (weighted Lambda and alternatives).
- Chalmers, R. P., & Adkins, M. C. (2020). Writing effective and reliable
  Monte Carlo simulations with the SimDesign package. *The Quantitative
  Methods for Psychology, 16*(4), 248–280.
- Gwet, K. L. (2014). *Handbook of Inter-Rater Reliability* (4th ed.).
- Muthén, B. (1984). A general structural equation model with dichotomous,
  ordered categorical, and continuous latent variable indicators.
  *Psychometrika, 49*(1), 115–132.
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
  studies to evaluate statistical methods. *Statistics in Medicine, 38*(11),
  2074–2102.
- Siepe, B. S., Bartoš, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., &
  Pawel, S. (2024). Simulation studies for methodological research in
  psychology: A standardized template for planning, preregistration, and
  reporting. *Psychological Methods*.
- Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
  assessing rater reliability. *Psychological Bulletin, 86*(2), 420–428.
