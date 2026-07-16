# =====================================================================
# 01_functions.R
# Helpers for the SimDesign simulation: response distributions,
# analytic true (population) values, and coefficient estimators.
# All estimator calls use :: so that SimDesign's `packages` argument
# can load dependencies on parallel workers.
# =====================================================================

# Coefficient labels used throughout Analyse/Summarise/post-processing.
# Defined as a function (not a global constant) because SimDesign exports
# workspace *functions* to parallel workers but not plain data objects —
# a global vector here raises "object not found" on the workers.
coef_labels <- function() {
  c("PA", "Kappa", "Kappa_quad", "Fleiss",
    "AC1", "AC2_quad", "Alpha_ord", "BP", "ICC21")
}

# Display names for tables/figures
coef_display <- function(coefficient) {
  c(PA         = "Percent agreement",
    Kappa      = "Cohen/Conger kappa",
    Kappa_quad = "Quadratic-weighted kappa (QWK)",
    Fleiss     = "Fleiss kappa",
    AC1        = "Gwet AC1",
    AC2_quad   = "Gwet AC2 (quadratic)",
    Alpha_ord  = "Krippendorff alpha (ordinal)",
    BP         = "Brennan-Prediger",
    ICC21      = "ICC(2,1)")[coefficient]
}

# Condition-dependent applicability of the coefficient set:
# with nLevels == 2 the weighted variants are mathematically identical to
# their unweighted counterparts (quadratic-weighted kappa == kappa,
# AC2 == AC1), so reporting them would duplicate the same estimator.
# Post-processing and figures keep only applicable == TRUE rows.
# (Conger's kappa reduces exactly to Cohen's kappa at k == 2, and
# Krippendorff's ordinal alpha reduces to nominal alpha at nLevels == 2,
# so those stay applicable everywhere.)
coef_applicable <- function(coefficient, nLevels) {
  ifelse(coefficient %in% c("Kappa_quad", "AC2_quad"), nLevels > 2, TRUE)
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
#   * Cohen/Conger/Fleiss kappa (any weights), Krippendorff's alpha
#     (any metric), and the intraclass correlation all equal `agree`.
#   * Percent agreement and Gwet/Brennan-Prediger coefficients do NOT
#     equal `agree` in general; they are computed below.
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

# Brennan-Prediger; chance agreement sum(W)/L^2
true_bp <- function(agree, p, weights = "unweighted") {
  L <- length(p)
  W <- agreement_weights(L, weights)
  pa <- sum(W * pair_joint(agree, p))
  pe <- sum(W) / L^2
  (pa - pe) / (1 - pe)
}

# True value dispatcher, keyed by COEF_LABELS. Kappa family, alpha,
# and ICC all equal `agree`.
true_value <- function(coefficient, agree, nLevels, prob_type) {
  p <- get_probs(prob_type, nLevels)
  switch(coefficient,
         PA       = true_pa(agree, p),
         AC1      = true_gwet(agree, p, "unweighted"),
         AC2_quad = true_gwet(agree, p, "quadratic"),
         BP       = true_bp(agree, p, "unweighted"),
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
# H0: coefficient = 0; NA where no test exists). Failed fits yield NA
# rather than throwing, so SimDesign does NOT redraw the data on
# estimation failure — redraws would condition results on estimability
# and bias the small-sample conditions. The NA rate is reported as the
# convergence rate in Summarise(). categ.labels = 1:nLevels keeps
# unobserved categories from shrinking the category space (this matters
# for the weighted coefficients).
# ---------------------------------------------------------------------

fit_all_coefficients <- function(x, nLevels) {

  xdf  <- as.data.frame(x)
  labs <- seq_len(nLevels)

  fits <- list(
    PA        = safe_fit(irrCAC::pa.coeff.raw(xdf, categ.labels = labs)),
    Kappa     = safe_fit(irrCAC::conger.kappa.raw(xdf, weights = "unweighted",
                                                  categ.labels = labs)),
    Kappa_quad = safe_fit(irrCAC::conger.kappa.raw(xdf, weights = "quadratic",
                                                   categ.labels = labs)),
    Fleiss    = safe_fit(irrCAC::fleiss.kappa.raw(xdf, weights = "unweighted",
                                                  categ.labels = labs)),
    AC1       = safe_fit(irrCAC::gwet.ac1.raw(xdf, weights = "unweighted",
                                              categ.labels = labs)),
    AC2_quad  = safe_fit(irrCAC::gwet.ac1.raw(xdf, weights = "quadratic",
                                              categ.labels = labs)),
    Alpha_ord = safe_fit(irrCAC::krippen.alpha.raw(xdf, weights = "ordinal",
                                                   categ.labels = labs)),
    BP        = safe_fit(irrCAC::bp.coeff.raw(xdf, weights = "unweighted",
                                              categ.labels = labs))
  )

  res <- vapply(fits, cac_est, numeric(2))
  ests  <- res["est", ]
  pvals <- res["p", ]
  pvals["PA"] <- NA_real_  # H0: PA = 0 is not a meaningful test

  # ICC(2,1): two-way random effects, absolute agreement, single rater
  icc_fit <- safe_fit(irr::icc(as.matrix(x), model = "twoway",
                               type = "agreement", unit = "single"))
  ests  <- c(ests,  ICC21 = if (is.null(icc_fit)) NA_real_ else num1(icc_fit$value))
  pvals <- c(pvals, ICC21 = if (is.null(icc_fit)) NA_real_ else num1(icc_fit$p.value))

  stopifnot(identical(names(ests), coef_labels()))

  c(stats::setNames(ests,  paste0("est.", names(ests))),
    stats::setNames(pvals, paste0("p.",  names(pvals))))
}
