# =====================================================================
# 01_functions.R
# Helpers for the SimDesign simulation: response distributions,
# analytic true (population) values, and coefficient estimators.
#
# Research question: do weighted agreement coefficients outperform
# their unweighted ("default") counterparts on ordinal scales?
# Two coefficient families x weighting schemes:
#   Kappa family : Cohen's kappa (unweighted) | quadratic (QWK)
#   Gwet family  : AC1 (unweighted) | AC2-linear | AC2-quadratic
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
# Response probability distributions
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

# ---------------------------------------------------------------------
# True (population) values under the IRRsim data-generating process
#
# simulateRatingMatrix() works per event as: draw a seed score X ~ p;
# with probability `agree` ALL raters score X, otherwise every rater
# scores independently from p. Hence for any pair of raters the joint
# distribution of scores is
#     P(i, j) = agree * p_i * 1{i == j} + (1 - agree) * p_i * p_j
# From this all population coefficient values follow analytically
# (derivations in README.md):
#   * Cohen's kappa with ANY weights equals `agree` (every disagreement
#     comes from the independence component, so weighted observed
#     disagreement = (1 - agree) x weighted chance disagreement).
#   * Gwet's AC1/AC2 use a different chance correction; their true
#     values are computed below and differ from `agree` under skew.
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

# Population pairwise joint distribution of two raters' scores
pair_joint <- function(agree, p) {
  P <- (1 - agree) * outer(p, p)
  diag(P) <- diag(P) + agree * p
  P
}

# Expected raw percent agreement (reference only; not an estimand here)
true_pa <- function(agree, p) {
  agree + (1 - agree) * sum(p^2)
}

# Gwet's AC1 (unweighted) / AC2 (weighted); Gwet (2014) chance agreement
true_gwet <- function(agree, p, weights = "unweighted") {
  L <- length(p)
  W <- agreement_weights(L, weights)
  pa <- sum(W * pair_joint(agree, p))
  pe <- sum(W) * sum(p * (1 - p)) / (L * (L - 1))
  (pa - pe) / (1 - pe)
}

# True value dispatcher, keyed by coef_labels(). The kappa family
# equals `agree` for every weighting scheme.
true_value <- function(coefficient, agree, nLevels, prob_type) {
  p <- get_probs(prob_type, nLevels)
  switch(coefficient,
         AC1        = true_gwet(agree, p, "unweighted"),
         AC2_linear = true_gwet(agree, p, "linear"),
         AC2_quad   = true_gwet(agree, p, "quadratic"),
         agree)
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
