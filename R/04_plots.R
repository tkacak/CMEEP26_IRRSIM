# =====================================================================
# 04_plots.R
# Figures for the weighted-vs-unweighted comparison: relative bias,
# bias, empirical SE, and power curves. Panels are nLevels x coefficient
# family; colour/shape encode the weighting scheme, so each panel
# directly contrasts unweighted vs linear vs quadratic within a family.
# One PNG per (nEvents x prob_type) cell.
# Run with the R/ folder as the working directory.
# =====================================================================

library(dplyr)
library(ggplot2)

out_dir <- "results"
fig_dir <- "figs"
dir.create(fig_dir, showWarnings = FALSE)

# applicable == FALSE rows are weighted variants at nLevels == 2, where
# they are mathematically identical to their unweighted counterparts
perf <- readRDS(file.path(out_dir, "performance.rds")) %>%
  filter(applicable)

plot_measure <- function(data, y, ylab, ref = 0,
                         band = c(-0.05, 0.05)) {
  ggplot(data,
         aes(x = agree, y = .data[[y]],
             colour = weighting, shape = weighting,
             group = weighting)) +
    { if (!is.null(band)) geom_hline(yintercept = band, linetype = "dotted") } +
    geom_hline(yintercept = ref, linetype = "dashed") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    facet_grid(nLevels ~ family,
               labeller = labeller(nLevels = function(x) paste0("L = ", x))) +
    scale_x_continuous(breaks = seq(0.3, 0.9, 0.2)) +
    labs(x = "True agreement (agree)", y = ylab,
         colour = "Weighting", shape = "Weighting") +
    theme_bw() +
    theme(legend.position = "bottom")
}

for (pt in unique(perf$prob_type)) {
  for (ne in unique(perf$nEvents)) {

    d <- perf %>% filter(prob_type == pt, nEvents == ne, agree > 0)
    tag <- sprintf("%s_n%02d", pt, ne)

    ggsave(file.path(fig_dir, paste0("relbias_", tag, ".png")),
           plot_measure(d, "RelBias",
                        "Relative bias", band = c(-0.05, 0.05)),
           width = 8, height = 8, dpi = 300)

    ggsave(file.path(fig_dir, paste0("bias_", tag, ".png")),
           plot_measure(d, "Bias", "Bias", band = NULL),
           width = 8, height = 8, dpi = 300)

    ggsave(file.path(fig_dir, paste0("empse_", tag, ".png")),
           plot_measure(d, "EmpSE", "Empirical SE", band = NULL),
           width = 8, height = 8, dpi = 300)

    ggsave(file.path(fig_dir, paste0("power_", tag, ".png")),
           plot_measure(d %>% filter(rejection_type == "Power"),
                        "Rejection",
                        "Power (H0: coefficient = 0)",
                        ref = 0.80, band = NULL),
           width = 8, height = 8, dpi = 300)
  }
}

# Type I error figure (agree = 0 conditions, where truth == 0).
# Under skewed margins the AC family has truth != 0 at agree = 0, so
# those cells appear only in the uniform columns.
d0 <- perf %>% filter(rejection_type == "TypeI")
if (nrow(d0) > 0) {
  p0 <- ggplot(d0, aes(x = factor(nEvents), y = Rejection,
                       colour = weighting, shape = weighting,
                       group = weighting)) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.8) +
    facet_grid(nLevels ~ family + prob_type,
               labeller = labeller(nLevels = function(x) paste0("L = ", x))) +
    labs(x = "nEvents", y = "Type I error (nominal .05)",
         colour = "Weighting", shape = "Weighting") +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(file.path(fig_dir, "type1_error.png"), p0,
         width = 12, height = 8, dpi = 300)
}

cat("Figures written to", fig_dir, "\n")
