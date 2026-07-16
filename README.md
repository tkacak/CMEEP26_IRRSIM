# CMEEP26_IRRSIM — Comparing Inter-Rater Agreement Coefficients

Bu klasörde Arş. Gör. Tugay KAÇAK ve Doç. Dr. Abdullah Faruk KILIÇ'ın CMEEP26'da sunduğu puanlayıcılar arası uyum katsayılarının karşılaştırıldığı çalışmanın kodları bulunmaktadır.

A Monte Carlo simulation comparing the performance of chance-corrected
inter-rater agreement coefficients, implemented in the
[`SimDesign`](https://cran.r-project.org/package=SimDesign) framework
(Chalmers & Adkins, 2020) with the data-generating process of the
[`IRRsim`](https://irrsim.bryer.org) package (Bryer). Reporting follows
the simulation template of Siepe et al. (2024) with performance measures and
Monte Carlo standard errors (MCSE) from Morris, White & Crowther (2019).

## Design (ADEMP)

**Aims.** Compare bias, relative bias, and power (plus Type I error) of
commonly used agreement coefficients for ordinal rating scales across small
rating designs.

**Data-generating mechanism.** `IRRsim::simulateRatingMatrix()`. For each
event a seed score is drawn from the response distribution *p*; with
probability `agree` **all** raters give the seed score, otherwise all raters
score independently from *p*. Full factorial design:

| Factor      | Levels                                    |
|-------------|-------------------------------------------|
| `nLevels`   | 2, 3, 4, 5                                 |
| `k` (raters)| 2, 3, 4, 5 (complete design, `k_per_event = k`) |
| `agree`     | 0.30–0.90 by 0.10, plus 0.00 (null)        |
| `nEvents`   | 20, 30, 40, 50                             |
| `prob_type` | uniform; skew (*p_i ∝ i*)                  |

8 × 4 × 4 × 4 × 2 = **1024 conditions**, `n_reps = 1000` replications each.

**Estimands / true values.** Under this DGP the pairwise joint distribution
of any two raters' scores is
`P(i, j) = agree · p_i · 1{i=j} + (1 − agree) · p_i · p_j`,
i.e., a mixture of perfect agreement (prob. `agree`) and independence.
Because every disagreement arises only from the independence component, for
*any* agreement weights the expected weighted observed disagreement equals
`(1 − agree)` times the expected chance disagreement computed from the
(identical) marginals. Hence the population values are:

| Coefficient                                   | True value |
|-----------------------------------------------|------------|
| Cohen/Conger kappa, Fleiss kappa (any weights)| `agree`    |
| Krippendorff's alpha (any metric)             | `agree` (asymptotically) |
| ICC(2,1) agreement                            | `agree` (pairwise product-moment correlation = `agree`) |
| Percent agreement                             | `agree + (1 − agree) Σp²` |
| Gwet AC1 / AC2, Brennan–Prediger              | closed form; = `agree` under uniform *p*, ≠ `agree` under skew (see `true_value()` in `R/01_functions.R`) |

This is why performance must be evaluated against **each coefficient's own
estimand**, not against `agree` for all of them.

**Methods (coefficients compared).**

| Label        | Coefficient                              | Implementation |
|--------------|-------------------------------------------|----------------|
| `PA`         | Percent agreement (descriptive baseline) | `irrCAC::pa.coeff.raw` |
| `Kappa`      | Conger's kappa (= Cohen's kappa at k = 2)| `irrCAC::conger.kappa.raw` |
| `Kappa_quad` | Quadratic-weighted kappa (= QWK)         | `irrCAC::conger.kappa.raw(weights = "quadratic")` |
| `Fleiss`     | Fleiss' kappa                            | `irrCAC::fleiss.kappa.raw` |
| `AC1`        | Gwet's AC1                               | `irrCAC::gwet.ac1.raw` |
| `AC2_quad`   | Gwet's AC2, quadratic weights            | `irrCAC::gwet.ac1.raw(weights = "quadratic")` |
| `Alpha_ord`  | Krippendorff's alpha, ordinal weights    | `irrCAC::krippen.alpha.raw(weights = "ordinal")` |
| `BP`         | Brennan–Prediger                         | `irrCAC::bp.coeff.raw` |
| `ICC21`      | ICC(2,1), two-way random, agreement      | `irr::icc` |

All `irrCAC` estimators receive `categ.labels = 1:nLevels` so unobserved
categories in small samples cannot shrink the category space, and they all
report a subject-level-variance test of H0: coefficient = 0 (used for
power). `Metrics::ScoreQuadraticWeightedKappa` from the v1 script was
dropped: it estimates the same quantity as `Kappa_quad` but without a test.

Estimation failures are caught *inside* `Analyse()` and returned as `NA`
rather than thrown, so SimDesign does **not** redraw the data on failure —
redraws would condition results on estimability and bias the small-sample
conditions. The NA rate is reported as the convergence rate.

**Performance measures** (per condition × coefficient, with MCSE; S =
converged replications):

- Bias = mean(est) − truth; MCSE = √(var(est)/S)
- Relative bias = Bias/truth; MCSE = MCSE(Bias)/|truth|
- Empirical SE; MCSE = EmpSE/√(2(S−1))
- RMSE; MCSE via delta method = MCSE(MSE)/(2·RMSE)
- Rejection rate of H0: coefficient = 0 at α = .05; MCSE = √(R(1−R)/S).
  This is **power** where the truth ≠ 0 and **Type I error** where the truth
  = 0 (the `agree = 0` conditions — but note under skew AC1/BP/PA have
  non-zero truth even at `agree = 0`, which `rejection_type` accounts for).
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

```r
install.packages(c("SimDesign", "IRRsim", "irr", "irrCAC",
                   "dplyr", "tidyr", "ggplot2"))
source("R/02_run_simulation.R")  # writes results/irr_simdesign.rds (+ raw-results/)
source("R/03_summarise.R")       # writes results/performance.rds / .csv
source("R/04_plots.R")           # writes figs/*.png
```

`runSimulation()` handles parallelisation (`parallel = TRUE`), per-condition
seeds (`seed = 20260716 + 1:1024`, so the run is fully reproducible), and
progress/ETA reporting. If the run is interrupted, re-sourcing the script
resumes automatically from SimDesign's tempfile. Replication-level raw
results are stored in `results/raw-results/` (`save_results = TRUE`), the
condition-level summary in `results/irr_simdesign.rds`, and error/warning
bookkeeping is available via `SimExtract(res, "errors")`. `sessionInfo()`
is saved alongside the results, per the Siepe et al. checklist.

## References

- Bryer, J. *IRRsim: Simulating Inter-Rater Reliability*. https://irrsim.bryer.org
- Chalmers, R. P., & Adkins, M. C. (2020). Writing effective and reliable
  Monte Carlo simulations with the SimDesign package. *The Quantitative
  Methods for Psychology, 16*(4), 248–280.
- Gwet, K. L. (2014). *Handbook of Inter-Rater Reliability* (4th ed.).
- Morris, T. P., White, I. R., & Crowther, M. J. (2019). Using simulation
  studies to evaluate statistical methods. *Statistics in Medicine, 38*(11),
  2074–2102.
- Siepe, B. S., Bartoš, F., Morris, T. P., Boulesteix, A.-L., Heck, D. W., &
  Pawel, S. (2024). Simulation studies for methodological research in
  psychology: A standardized template for planning, preregistration, and
  reporting. *Psychological Methods*.
- Shrout, P. E., & Fleiss, J. L. (1979). Intraclass correlations: Uses in
  assessing rater reliability. *Psychological Bulletin, 86*(2), 420–428.
