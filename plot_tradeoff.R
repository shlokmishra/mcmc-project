# plot_tradeoff.R
library(ggplot2)
library(gridExtra)

load("loglinear_datageneratingparameters.RData")

analyze_results <- function() {
  files <- list.files(pattern = "tradeoff_eps.*\\.rds")
  
  metrics <- lapply(files, function(f) {
    res <- readRDS(f)
    data.frame(
      epsilon = res$epsilon,
      method = ifelse(grepl("adaptive", f), "adaptive", res$method),
      bias = mean(abs(colMeans(res$p_samples) - true_params$p)),
      coverage = mean(sapply(1:5, function(i) {
        quantile(res$p_samples[,i], c(0.05, 0.95))[1] <= true_params$p[i] &&
          quantile(res$p_samples[,i], c(0.05, 0.95))[2] >= true_params$p[i]
      })),
      ess = mean(effectiveSize(res$p_samples)),
      runtime = res$runtime
    )
  }) %>% do.call(rbind, .)
  
  # Generate plots
  p1 <- ggplot(metrics, aes(epsilon, bias, color=method)) + 
    geom_line() + scale_x_log10() + labs(title="Bias vs Privacy Level")
  
  p2 <- ggplot(metrics, aes(epsilon, coverage, color=method)) +
    geom_line() + geom_hline(yintercept=0.9, linetype=2) +
    scale_x_log10() + labs(title="Coverage vs Privacy Level")
  
  p3 <- ggplot(metrics, aes(epsilon, ess/runtime, color=method)) +
    geom_line() + scale_x_log10() + labs(title="ESS/Second Efficiency")
  
  grid.arrange(p1, p2, p3, ncol=1)
}

analyze_results()
ggsave("privacy_tradeoff_analysis.png", width=10, height=8)
