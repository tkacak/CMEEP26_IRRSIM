# =====================================================================
# 04_plots.R
# Figures: relative bias, bias (with MCSE bands), and power curves.
# One PNG per (nEvents x prob_type) cell; panels are nLevels x k.
# =====================================================================

library(dplyr)
library(ggplot2)

out_dir <- "results"
fig_dir <- "figs"
dir.create(fig_dir, showWarnings = FALSE)

perf <- readRDS(file.path(out_dir, "performance.rds"))

plot_measure <- function(data, y, y_mcse, ylab, ref = 0,
                         band = c(-0.05, 0.05)) {
  ggplot(data,
         aes(x = agree, y = .data[[y]],
             colour = coefficient, shape = coefficient,
             group = coefficient)) +
    { if (!is.null(band)) geom_hline(yintercept = band, linetype = "dotted") } +
    geom_hline(yintercept = ref, linetype = "dashed") +
    geom_errorbar(aes(ymin = .data[[y]] - .data[[y_mcse]],
                      ymax = .data[[y]] + .data[[y_mcse]]),
                  width = 0.015, alpha = 0.6) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    facet_grid(nLevels ~ k,
               labeller = labeller(nLevels = function(x) paste0("L = ", x),
                                   k = function(x) paste0("k = ", x))) +
    scale_x_continuous(breaks = seq(0.3, 0.9, 0.2)) +
    labs(x = "True agreement (agree)", y = ylab,
         colour = NULL, shape = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
}

for (pt in unique(perf$prob_type)) {
  for (ne in unique(perf$nEvents)) {

    d <- perf %>% filter(prob_type == pt, nEvents == ne, agree > 0)
    tag <- sprintf("%s_n%02d", pt, ne)

    ggsave(file.path(fig_dir, paste0("relbias_", tag, ".png")),
           plot_measure(d, "RelBias", "RelBias_MCSE",
                        "Relative bias", band = c(-0.05, 0.05)),
           width = 10, height = 8, dpi = 300)

    ggsave(file.path(fig_dir, paste0("bias_", tag, ".png")),
           plot_measure(d, "Bias", "Bias_MCSE", "Bias", band = NULL),
           width = 10, height = 8, dpi = 300)

    ggsave(file.path(fig_dir, paste0("power_", tag, ".png")),
           plot_measure(d %>% filter(rejection_type == "Power"),
                        "Rejection", "Rejection_MCSE",
                        "Power (H0: coefficient = 0)",
                        ref = 0.80, band = NULL),
           width = 10, height = 8, dpi = 300)
  }
}

# Type I error figure (agree = 0 conditions, where truth == 0)
d0 <- perf %>% filter(rejection_type == "TypeI")
if (nrow(d0) > 0) {
  p0 <- ggplot(d0, aes(x = factor(nEvents), y = Rejection,
                       colour = coefficient, shape = coefficient,
                       group = coefficient)) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_errorbar(aes(ymin = Rejection - Rejection_MCSE,
                      ymax = Rejection + Rejection_MCSE),
                  width = 0.15, alpha = 0.6) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    facet_grid(nLevels ~ k + prob_type) +
    labs(x = "nEvents", y = "Type I error (nominal .05)",
         colour = NULL, shape = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "type1_error.png"), p0,
         width = 12, height = 8, dpi = 300)
}

cat("Figures written to", fig_dir, "\n")
