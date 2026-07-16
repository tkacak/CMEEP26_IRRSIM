library(irr)
library(IRRsim)
library(irrCAC)
library(dplyr)
library(tidyr)
library(Metrics)

#--------------------------------------------------
# Helper functions
#--------------------------------------------------

safe_num <- function(x) {
  out <- tryCatch(suppressWarnings(x), error = function(e) NA_real_)
  if (is.null(out) || length(out) == 0) return(NA_real_)
  as.numeric(out)[1]
}

safe_prop <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# Skewed probability (y??ksek kategori bias)
skew_probs <- function(nLevels) {
  w <- 1:nLevels
  w / sum(w)
}

#--------------------------------------------------
# Simulation function
#--------------------------------------------------

run_one_sim <- function(nLevels, k, agree, nEvents, prob_type) {
  
  probs <- switch(prob_type,
                  "uniform" = rep(1/nLevels, nLevels),
                  "skew"    = skew_probs(nLevels))
  
  x <- simulateRatingMatrix(
    nLevels = nLevels,
    k = k,
    k_per_event = k,
    agree = agree,
    nEvents = nEvents,
    response.probs = probs
  )
  
  #---------------- ICC ----------------
  icc_out <- tryCatch(
    irr::icc(x, model="twoway", type="agreement"),
    error = function(e) NULL
  )
  
  ICC <- if (!is.null(icc_out)) safe_num(icc_out$value) else NA
  ICC_p <- if (!is.null(icc_out)) safe_num(icc_out$p.value) else NA
  
  #---------------- Krippendorff ----------------
  kripp_out <- tryCatch(
    irr::kripp.alpha(t(x), method="ordinal"),
    error = function(e) NULL
  )
  
  Kripp <- if (!is.null(kripp_out)) safe_num(kripp_out$value) else NA
  
  #---------------- Gwet AC1 Linear ----------------
  gwet_lin <- tryCatch(
    irrCAC::gwet.ac1.raw(x, weights="linear"),
    error = function(e) NULL
  )
  
  Gwet_L <- if (!is.null(gwet_lin)) safe_num(gwet_lin$est$coeff.val) else NA
  Gwet_L_p <- if (!is.null(gwet_lin)) safe_num(gwet_lin$est$p.value) else NA
  
  #---------------- Gwet AC1 Quadratic ----------------
  gwet_quad <- tryCatch(
    irrCAC::gwet.ac1.raw(x, weights="quadratic"),
    error = function(e) NULL
  )
  
  Gwet_Q <- if (!is.null(gwet_quad)) safe_num(gwet_quad$est$coeff.val) else NA
  Gwet_Q_p <- if (!is.null(gwet_quad)) safe_num(gwet_quad$est$p.value) else NA
  
  #---------------- Cohen Kappa ----------------
  coh <- tryCatch(
    irr::kappa2(x[,1:2], weight="unweighted"),
    error = function(e) NULL
  )
  
  Cohen <- if (!is.null(coh)) safe_num(coh$value) else NA
  Cohen_p <- if (!is.null(coh)) safe_num(coh$p.value) else NA
  
  #---------------- Weighted Cohen Kappa ----------------
  coh_s <- tryCatch(
    irr::kappa2(x[,1:2], weight="equal"),
    error = function(e) NULL
  )
  
  Cohen_s <- if (!is.null(coh_s)) safe_num(coh_s$value) else NA
  Cohen_s_p <- if (!is.null(coh_s)) safe_num(coh_s$p.value) else NA
  
  #---------------- QWK ----------------
  QWK <- tryCatch(
    Metrics::ScoreQuadraticWeightedKappa(
      x[,1], x[,2],
      min.rating = min(x[,1],x[,2]),
      max.rating = max(x[,1],x[,2])
    ),
    error = function(e) NA
  )
  
  data.frame(
    ICC, ICC_p,
    Kripp,
    Gwet_L, Gwet_L_p,
    Gwet_Q, Gwet_Q_p,
    Cohen, Cohen_p,
    Cohen_s, Cohen_s_p,
    QWK
  )
}

#--------------------------------------------------
# Simulation design
#--------------------------------------------------

nLevels_vec <- 2:5
k_vec <- 2
agree_vec <- seq(0.30, 0.90, 0.10)
nEvents_vec <- seq(20, 50, 10)

# ??? NEW: distribution types
prob_type_vec <- c("uniform", "skew")

n_reps <- 1000

param_grid <- expand.grid(
  nLevels = nLevels_vec,
  k = k_vec,
  agree = agree_vec,
  nEvents = nEvents_vec,
  prob_type = prob_type_vec
)

set.seed(123)

#--------------------------------------------------
# Run simulation
#--------------------------------------------------

sim_list <- lapply(seq_len(nrow(param_grid)), function(i) {
  
  if (i %% 10 == 0) {
    cat("Progress:", i, "/", nrow(param_grid), "\n")
  }
  
  one_param <- lapply(seq_len(n_reps), function(j) {
    run_one_sim(
      param_grid$nLevels[i],
      param_grid$k[i],
      param_grid$agree[i],
      param_grid$nEvents[i],
      param_grid$prob_type[i]
    )
  })
  
  df <- bind_rows(one_param)
  
  df$nLevels <- param_grid$nLevels[i]
  df$agree <- param_grid$agree[i]
  df$nEvents <- param_grid$nEvents[i]
  df$prob_type <- param_grid$prob_type[i]
  
  df
})

sim_df <- bind_rows(sim_list)

#--------------------------------------------------
# Summary function
#--------------------------------------------------

summarise_metric <- function(data, name, est, p=NULL) {
  
  data %>%
    group_by(nLevels, agree, nEvents, prob_type) %>%
    summarise(
      Metric = name,
      Target = first(agree),
      
      Mean = mean(.data[[est]], na.rm=TRUE),
      Bias = Mean - Target,
      Variance = var(.data[[est]], na.rm=TRUE),
      RMSE = sqrt(mean((.data[[est]] - Target)^2, na.rm=TRUE)),
      Relative_Bias = Bias / Target,
      
      TypeI = if (!is.null(p)) safe_prop(.data[[p]] < 0.05) else NA_real_,
      
      .groups="drop"
    )
}

#--------------------------------------------------
# Final report
#--------------------------------------------------

metric_report <- bind_rows(
  
  summarise_metric(sim_df, "ICC", "ICC", "ICC_p"),
  summarise_metric(sim_df, "Gwet AC1 Linear", "Gwet_L", "Gwet_L_p"),
  summarise_metric(sim_df, "Gwet AC1 Quadratic", "Gwet_Q", "Gwet_Q_p"),
  summarise_metric(sim_df, "Cohen Kappa", "Cohen", "Cohen_p"),
  summarise_metric(sim_df, "Weighted Cohen Kappa", "Cohen_s", "Cohen_s_p"),
  summarise_metric(sim_df, "Krippendorff Alpha", "Kripp", NULL),
  summarise_metric(sim_df, "QWK", "QWK", NULL)
  
)

#--------------------------------------------------
# Example outputs
#--------------------------------------------------

head(metric_report)

# Overall summary
overall_summary <- metric_report %>%
  group_by(Metric, agree,nEvents,prob_type) %>%
  summarise(
    Mean_Bias = mean(Bias, na.rm=TRUE),
    Mean_RMSE = mean(RMSE, na.rm=TRUE),
    Mean_TypeI = mean(TypeI, na.rm=TRUE),
    Mean_RB = mean(Relative_Bias,na.rm = T),
    .groups="drop"
  )

print(overall_summary)

ggplot(overall_summary,
       aes(x = as.factor(agree),
           y = Mean_RB,
           color = Metric,
           group = Metric,
           shape = Metric)) +
  
  #geom_hline(yintercept = c(0.05,-0.05), linetype = "dashed", color = "black") +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  facet_grid(nEvents~prob_type) +
  
  theme_bw()
