library(ggplot2)
library(gridExtra)


plot_comparison <- function(eps) {
  cat(sprintf("Loading RDS for epsilon = %.1f...\n", eps))
  
  # Check files exist
  orig_file <- sprintf("results_eps%.1f_original.rds", eps)
  smart_file <- sprintf("results_eps%.1f_smart.rds", eps)
  
  if (!file.exists(orig_file)) stop(paste("Missing:", orig_file))
  if (!file.exists(smart_file)) stop(paste("Missing:", smart_file))
  
  orig <- readRDS(orig_file)
  smart <- readRDS(smart_file)
  
  # Validate structure
  if (is.null(orig$p_samples) || is.null(smart$p_samples)) stop("Missing p_samples in RDS")
  
  ess_orig <- tryCatch(effectiveSize(orig$p_samples), error = function(e) { print(e); return(rep(NA, 5)) })
  ess_smart <- tryCatch(effectiveSize(smart$p_samples), error = function(e) { print(e); return(rep(NA, 5)) })
  
  # Create ESS comparison data
  ess_df <- data.frame(
    Method = rep(c("Original", "Smart"), each = 5),
    Parameter = rep(1:5, 2),
    ESS = c(ess_orig, ess_smart)
  )
  
  # Create acceptance vs ESS data
  accept_df <- data.frame(
    Method = c("Original", "Smart"),
    Acceptance = c(mean(orig$acc_rates), mean(smart$acc_rates)),
    ESS = c(mean(ess_orig, na.rm = TRUE), mean(ess_smart, na.rm = TRUE))
  )
  
  # ESS Comparison Plot
  p1 <- ggplot(ess_df, aes(x = factor(Parameter), y = ESS, fill = Method)) +
    geom_col(position = position_dodge()) +
    scale_fill_brewer(palette = "Set1") +
    labs(title = sprintf("ESS Comparison (ε=%.1f)", eps), x = "Parameter")
  
  # Acceptance vs ESS Plot
  p2 <- ggplot(accept_df, aes(x = Acceptance, y = ESS, color = Method)) +
    geom_point(size = 4) +
    scale_color_brewer(palette = "Set1") +
    labs(title = "Acceptance vs ESS Tradeoff")
  
  list(p1 = p1, p2 = p2)
}



for (eps in c(0.1, 0.3, 1, 3, 10)) {
  plots <- plot_comparison(eps)
  
  if (inherits(plots$p1, "gg")) {
    ggsave(sprintf("ess_comparison_eps%.1f.png", eps), plot = plots$p1, width = 6, height = 4)
  } else {
    cat(sprintf("Skipping ESS plot for eps=%.1f\n", eps))
  }
  
  if (inherits(plots$p2, "gg")) {
    ggsave(sprintf("acceptance_vs_ess_eps%.1f.png", eps), plot = plots$p2, width = 6, height = 4)
  } else {
    cat(sprintf("Skipping acceptance vs ESS plot for eps=%.1f\n", eps))
  }
}
