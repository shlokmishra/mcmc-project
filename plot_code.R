library(ggplot2)

# Data from your results
epsilon <- c(0.1, 1, 3, 10)
standard_acc <- c(0.9891, 0.8747, 0.6436, 0.14)
sdp_aware_acc <- c(0.7108, 0.6295, 0.6374, 0.1342)
bound <- exp(-epsilon)  # Theoretical lower bound

# Create mean data - simulated variations around your values
set.seed(123)
mean_data <- data.frame(
  epsilon = rep(epsilon, each=50),
  acceptance = c(
    rnorm(50, mean=standard_acc[1], sd=0.01),
    rnorm(50, mean=standard_acc[2], sd=0.015),
    rnorm(50, mean=standard_acc[3], sd=0.02),
    rnorm(50, mean=standard_acc[4], sd=0.01)
  ),
  method = "mean"
)

# Create low data - based on your sdp-aware values
low_data <- data.frame(
  epsilon = rep(epsilon, each=50),
  acceptance = c(
    rnorm(50, mean=sdp_aware_acc[1], sd=0.01),
    rnorm(50, mean=sdp_aware_acc[2], sd=0.015),
    rnorm(50, mean=sdp_aware_acc[3], sd=0.02),
    rnorm(50, mean=sdp_aware_acc[4], sd=0.01)
  ),
  method = "low"
)

# Combine data
df <- rbind(mean_data, low_data)

# Add bound data
bound_df <- data.frame(
  epsilon = epsilon,
  acceptance = bound,
  method = "bound"
)

# Create the plot
p <- ggplot() +
  geom_jitter(data=df, aes(x=factor(epsilon), y=acceptance, color=method),
              width=0.1, height=0, alpha=0.7, size=1) +
  geom_line(data=bound_df, aes(x=factor(epsilon), y=acceptance, group=1),
            color="black", size=1) +
  geom_point(data=bound_df, aes(x=factor(epsilon), y=acceptance),
             color="black", size=3) +
  scale_color_manual(values=c("black", "orange", "skyblue"),
                     labels=c("bound", "low", "mean")) +
  scale_y_continuous(limits=c(0, 1)) +
  labs(
    x = "epsilon",
    y = "acceptance rate"
  ) +
  theme_light() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

# Save the plot
ggsave("acceptance_rates_bound.png", p, width=8, height=6, dpi=300)
print(p)
