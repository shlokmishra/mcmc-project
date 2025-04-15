library(ggplot2)

analyze_hybrid <- function(eps) {
  orig <- readRDS(sprintf("results_eps%.1f_original.rds", eps))
  smart <- readRDS(sprintf("results_eps%.1f_smart.rds", eps))
  hybrid <- readRDS(sprintf("results_eps%.1f_hybrid.rds", eps))
  
  data.frame(
    Method = rep(c("Original","Smart","Hybrid"), each=5),
    Parameter = rep(1:5,3),
    ESS = c(effectiveSize(orig$p_samples),
            effectiveSize(smart$p_samples),
            effectiveSize(hybrid$p_samples))
  ) %>% 
    ggplot(aes(factor(Parameter), ESS, fill=Method)) +
    geom_col(position=position_dodge()) +
    labs(title=sprintf("ESS Comparison (ε=%.1f)", eps), x="Parameter")
}

analyze_hybrid(0.1) # Repeat for 1.0 and 10.0
analyze_hybrid(1.0)
analyze_hybrid(10)
