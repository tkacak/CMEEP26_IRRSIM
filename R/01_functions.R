# =====================================================================
# 01_functions.R
# Helpers for the SimDesign simulation: UVA data generation, numeric
# true (population) values, and coefficient estimators.
#
# Research question: do weighted agreement coefficients outperform
# their unweighted ("default") counterparts on ordinal scales?
# Two coefficient families x weighting schemes:
#   Kappa family : Cohen's kappa (unweighted) | quadratic (QWK)
#   Gwet family  : AC1 (unweighted) | AC2-linear | AC2-quadratic
#
# Data-generating mechanism: Underlying Variable Approach (UVA;
# Muthen, 1984, as used by Almehrizi, 2025). Two standard-normal latent
# variables with correlation rho are discretised into L ordinal
# categories with a common set of thresholds. Unlike Almehrizi's flow,
# thresholds are FIXED design factors (balanced vs skewed prevalence)
# rather than randomly drawn per data set: random thresholds would make
# every replication estimate a different population value, which is
# incompatible with bias/power evaluation against known truths.
#
# All estimator calls use :: so that SimDesign's `packages` argument
# can load dependencies on parallel workers. All shared constants are
# functions (not global data objects) because SimDesign exports
# workspace *functions* to parallel workers but not plain data objects.
# =====================================================================

# Coefficient labels used throughout Analyse/Summarise/post-processing
coef_labels <- function() {
  c("Kappa", "Kappa_quad", "AC1", "AC2_linear", "AC2_quad")
}

# Display names for tables/figures
coef_display <- function(coefficient) {
  c(Kappa      = "Cohen's kappa",
    Kappa_quad = "Quadratic-weighted kappa (QWK)",
    AC1        = "Gwet AC1",
    AC2_linear = "Gwet AC2 (linear)",
    AC2_quad   = "Gwet AC2 (quadratic)")[coefficient]
}

# Comparison factors: coefficient family and weighting scheme, so that
# post-processing can contrast weighted vs unweighted within family
coef_family <- function(coefficient) {
  ifelse(startsWith(coefficient, "Kappa"), "Kappa (Cohen)", "Gwet AC")
}

coef_weighting <- function(coefficient) {
  c(Kappa      = "unweighted",
    Kappa_quad = "quadratic",
    AC1        = "unweighted",
    AC2_linear = "linear",
    AC2_quad   = "quadratic")[coefficient]
}

# Weighted variants are mathematically identical to their unweighted
# counterparts when nLevels == 2 (all disagreement weights hit the same
# single off-diagonal band), so they would duplicate the same estimator.
# With the current design (nLevels 3-5) everything is applicable.
coef_applicable <- function(coefficient, nLevels) {
  ifelse(coefficient %in% c("Kappa_quad", "AC2_linear", "AC2_quad"),
         nLevels > 2, TRUE)
}

# ---------------------------------------------------------------------
# Category prevalence targets and UVA thresholds
#
# prob_type controls the population category prevalences via the
# discretisation thresholds (identical for both raters, so marginal
# distributions are symmetric):
#   uniform : all categories equally prevalent, p_i = 1/L
#   skew    : p_i proportional to i (upper categories more prevalent)
# Thresholds are the standard-normal quantiles of the cumulative
# prevalences, so the realised margins match get_probs() exactly.
# ---------------------------------------------------------------------

skew_probs <- function(nLevels) {
  w <- seq_len(nLevels)
  w / sum(w)
}

get_probs <- function(prob_type, nLevels) {
  switch(prob_type,
         uniform = rep(1 / nLevels, nLevels),
         skew    = skew_probs(nLevels),
         stop("Unknown prob_type: ", prob_type))
}

uva_thresholds <- function(prob_type, nLevels) {
  p <- get_probs(prob_type, nLevels)
  stats::qnorm(cumsum(p[-nLevels]))
}

# ---------------------------------------------------------------------
# UVA data generation for one condition (k = 2 raters)
#
# (Z1, Z2) ~ standard bivariate normal with correlation rho; both are
# discretised with the same thresholds. rho is the latent ("true")
# agreement parameter of the design.
# ---------------------------------------------------------------------

uva_generate <- function(nEvents, rho, nLevels, prob_type) {
  th <- uva_thresholds(prob_type, nLevels)
  z1 <- stats::rnorm(nEvents)
  z2 <- rho * z1 + sqrt(1 - rho^2) * stats::rnorm(nEvents)
  cbind(R1 = findInterval(z1, th) + 1L,
        R2 = findInterval(z2, th) + 1L)
}

# ---------------------------------------------------------------------
# True (population) values under the UVA data-generating process
#
# With rho and the thresholds fixed, the population L x L joint table
# of the two raters' categories is
#   P(i, j) = P(t_{i-1} < Z1 <= t_i, t_{j-1} < Z2 <= t_j)
# computed here by 1-D Gaussian quadrature (base R integrate(); no
# extra package needed):
#   P(i, j) = int_{t_{i-1}}^{t_i} phi(z) *
#             [Phi((t_j - rho z)/s) - Phi((t_{j-1} - rho z)/s)] dz,
#   s = sqrt(1 - rho^2).
# Every coefficient's population value follows by plugging P and the
# margins p into its defining formula. Note that under UVA the weighted
# and unweighted coefficients have genuinely different estimands, so
# each is evaluated against its own true value.
# ---------------------------------------------------------------------

agreement_weights <- function(L, weights = c("unweighted", "linear", "quadratic")) {
  weights <- match.arg(weights)
  i <- matrix(seq_len(L), L, L)
  j <- t(i)
  switch(weights,
         unweighted = (i == j) * 1,
         linear     = 1 - abs(i - j) / (L - 1),
         quadratic  = 1 - (i - j)^2 / (L - 1)^2)
}

uva_joint <- function(rho, thresholds) {
  L <- length(thresholds) + 1
  b <- c(-Inf, thresholds, Inf)
  p <- diff(stats::pnorm(b))
  if (abs(rho) < 1e-12) return(outer(p, p))
  s <- sqrt(1 - rho^2)
  P <- matrix(0, L, L)
  for (i in seq_len(L)) {
    for (j in seq_len(L)) {
      P[i, j] <- stats::integrate(
        function(z) stats::dnorm(z) *
          (stats::pnorm((b[j + 1] - rho * z) / s) -
           stats::pnorm((b[j]     - rho * z) / s)),
        lower = b[i], upper = b[i + 1],
        rel.tol = 1e-10, abs.tol = 1e-12
      )$value
    }
  }
  P / sum(P)  # remove tiny numerical drift
}

# Population Cohen's (weighted) kappa: chance = product of the margins
true_kappa_w <- function(P, p, W) {
  pa <- sum(W * P)
  pe <- sum(W * outer(p, p))
  (pa - pe) / (1 - pe)
}

# Population Gwet AC1/AC2: Gwet (2014) chance agreement
true_gwet_w <- function(P, p, W) {
  L  <- length(p)
  pa <- sum(W * P)
  pe <- sum(W) * sum(p * (1 - p)) / (L * (L - 1))
  (pa - pe) / (1 - pe)
}

# True value dispatcher, keyed by coef_labels()
true_value <- function(coefficient, rho, nLevels, prob_type) {
  p  <- get_probs(prob_type, nLevels)
  th <- uva_thresholds(prob_type, nLevels)
  P  <- uva_joint(rho, th)
  W  <- switch(coefficient,
               Kappa      = agreement_weights(nLevels, "unweighted"),
               Kappa_quad = agreement_weights(nLevels, "quadratic"),
               AC1        = agreement_weights(nLevels, "unweighted"),
               AC2_linear = agreement_weights(nLevels, "linear"),
               AC2_quad   = agreement_weights(nLevels, "quadratic"))
  if (startsWith(coefficient, "Kappa")) true_kappa_w(P, p, W)
  else true_gwet_w(P, p, W)
}

# ---------------------------------------------------------------------
# Estimation helpers
# ---------------------------------------------------------------------

safe_fit <- function(expr) {
  tryCatch(suppressWarnings(expr), error = function(e) NULL)
}

num1 <- function(v) {
  out <- suppressWarnings(as.numeric(as.character(v)))
  if (length(out) == 0) NA_real_ else out[1]
}

cac_est <- function(fit) {
  if (is.null(fit)) return(c(est = NA_real_, p = NA_real_))
  c(est = num1(fit$est$coeff.val), p = num1(fit$est$p.value))
}

# ---------------------------------------------------------------------
# Fit all coefficients on one rating matrix
#
# Returns a named numeric vector est.<coef> / p.<coef> (p = p-value of
# H0: coefficient = 0). Failed fits yield NA rather than throwing, so
# SimDesign does NOT redraw the data on estimation failure - redraws
# would condition results on estimability and bias the small-sample
# conditions. The NA rate is reported as the convergence rate.
# categ.labels = 1:nLevels keeps unobserved categories from shrinking
# the category space (this matters for the weighted coefficients).
# conger.kappa.raw() at k = 2 is exactly Cohen's (weighted) kappa.
# ---------------------------------------------------------------------

fit_all_coefficients <- function(x, nLevels) {

  xdf  <- as.data.frame(x)
  labs <- seq_len(nLevels)

  fits <- list(
    Kappa      = safe_fit(irrCAC::conger.kappa.raw(xdf, weights = "unweighted",
                                                   categ.labels = labs)),
    Kappa_quad = safe_fit(irrCAC::conger.kappa.raw(xdf, weights = "quadratic",
                                                   categ.labels = labs)),
    AC1        = safe_fit(irrCAC::gwet.ac1.raw(xdf, weights = "unweighted",
                                               categ.labels = labs)),
    AC2_linear = safe_fit(irrCAC::gwet.ac1.raw(xdf, weights = "linear",
                                               categ.labels = labs)),
    AC2_quad   = safe_fit(irrCAC::gwet.ac1.raw(xdf, weights = "quadratic",
                                               categ.labels = labs))
  )

  res <- vapply(fits, cac_est, numeric(2))
  ests  <- res["est", ]
  pvals <- res["p", ]

  stopifnot(identical(names(ests), coef_labels()))

  c(stats::setNames(ests,  paste0("est.", names(ests))),
    stats::setNames(pvals, paste0("p.",  names(pvals))))
}
