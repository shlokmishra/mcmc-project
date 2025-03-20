rm(list = ls());
library("extraDistr")
library("VGAM")
library("coda")
library("ggplot2")
library("reshape2")
library("gridExtra")
library("grid")

# Run parameters
niter <- 5000  # Number of Gibbs iterations per method
beta_param <- 0.1  # Tuning parameter for s_dp-aware proposals
epsvec <- c(0.1, 1, 3, 10)  # Different privacy levels to test
set.seed(123)  # For reproducibility

# Load data
load(file = "loglinear_datageneratingparameters.RData")  # Contains K, I, Jk, alpha, alp
load(file = "loglinear_notdpdata_N100.RData")  # Contains true data nijk, N

# Helper functions
get_n <- function(g){
  return(as.vector(table(factor(g, levels = 1 : I))))
}

get_nijk <- function(g, ans){
  nijk <- list()
  for(k in 1 : K){
    nijk[[k]] <- as.matrix(table(factor(ans[,k], levels = 1 : Jk[k]), factor(g, levels = 1 : I)))
  }
  return(nijk)
}

get_sse <- function(nijk, mijk){
  dis <- 0
  for(k in 1 : K){
    dis <- dis + sum((nijk[[k]] - mijk[[k]])^2)
  }
  return(dis)
}

# Create initial state
rgibbs_init <- function(){
  p <- rgamma(n = I, shape = alp, rate = 1)
  p <- p / sum(p)
  pijk <- list()
  for(k in 1: K){
    pijk_ <- matrix(0, nrow = Jk[k], ncol= I)
    for(i in 1 : I){
      pijk_[,i] <- rgamma(n = Jk[k], shape = alpha[[k]][,i], rate = 1)
      pijk_[,i] <- pijk_[,i] / sum(pijk_[,i])
    }
    pijk[[k]] <- pijk_
  }
  state <- list()
  state$p <- p
  state$pijk <- pijk
  # Sample individual group identities
  g <- sample.int(I, N, TRUE, p)
  ans <- matrix(nrow = N, ncol = K)
  for(k in 1 : K){
    for(id in 1 : N){
      ans[id,k] <- extraDistr::rcat(1, pijk[[k]][,g[id]])
    }
  }
  state$g <-  g
  state$ans <- ans
  state$n <- get_n(g)
  state$nijk <- get_nijk(g, ans)
  return(state)
}

# Standard Gibbs sampler (original approach)
gibbs_iter_standard <- function(state, mijk, epsilon){
  # Sample p given n from Dirichlet
  state$p <- extraDistr::rdirichlet(1, as.vector(state$n) + alp)
  for(k in 1 : K){
    for(i in 1 : I){
      state$pijk[[k]][,i] <- extraDistr::rdirichlet(1, state$nijk[[k]][,i] + alpha[[k]][,i])
    }
  }
  
  # Independent MH update for each g[id] and ans[id,]
  acnt <- 0
  for(id in 1 : N){  # For each person
    lratio <- 0
    # Standard blind proposal from current parameters
    newgid <- extraDistr::rcat(1, state$p)
    newans <- rep(NA, K)
    
    for(k in 1 : K){  # For each feature
      newans[k] <- extraDistr::rcat(1, state$pijk[[k]][,newgid])  # Sample new answer xk
      
      # Update ratio of eta's in log scale
      lratio <- lratio + extraDistr::dlaplace(
        state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid] + 1, 
        0, 2 * K / epsilon, TRUE
      )
      lratio <- lratio + extraDistr::dlaplace(
        state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]] - 1, 
        0, 2 * K / epsilon, TRUE
      )
      lratio <- lratio - extraDistr::dlaplace(
        state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid], 
        0, 2 * K / epsilon, TRUE
      )
      lratio <- lratio - extraDistr::dlaplace(
        state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]], 
        0, 2 * K / epsilon, TRUE
      )
    }
    
    if(log(runif(1)) < lratio){  # If accept then book keeping
      for(k in 1 : K){
        state$nijk[[k]][state$ans[id,k], state$g[id]] <- state$nijk[[k]][state$ans[id,k], state$g[id]] - 1
        state$nijk[[k]][newans[k], newgid] <- state$nijk[[k]][newans[k], newgid] + 1
      }
      state$n[state$g[id]] <- state$n[state$g[id]] - 1
      state$n[newgid] <- state$n[newgid] + 1
      state$g[id] <- newgid
      state$ans[id,] <- newans
      acnt <- acnt + 1  # Increment acceptance count
    }
    lratio <- 0
  }
  
  state$accept <- acnt / N
  state$dis <- get_sse(state$nijk, mijk)
  return(state)
}

# Modified Gibbs sampler with s_dp-aware proposals
gibbs_iter_sdp_aware <- function(state, mijk, epsilon, beta){
  # Sample p given n from Dirichlet
  state$p <- extraDistr::rdirichlet(1, as.vector(state$n) + alp)
  for(k in 1 : K){
    for(i in 1 : I){
      state$pijk[[k]][,i] <- extraDistr::rdirichlet(1, state$nijk[[k]][,i] + alpha[[k]][,i])
    }
  }
  
  # Independent MH update for each g[id] and ans[id,]
  acnt <- 0
  for(id in 1 : N){  # For each person
    lratio <- 0
    
    # Create s_dp-aware proposal distributions
    # Calculate weights based on discrepancy between mijk and nijk
    p_weights <- rep(0, I)
    for(i in 1:I) {
      weight_sum <- 0
      for(k in 1:K) {
        for(j in 1:Jk[k]) {
          diff <- mijk[[k]][j,i] - state$nijk[[k]][j,i]
          weight_sum <- weight_sum + exp(beta * diff)
        }
      }
      p_weights[i] <- state$p[i] * (weight_sum/(K*sum(Jk)))
    }
    p_weights <- p_weights / sum(p_weights)
    
    # Modified proposal for group
    newgid <- extraDistr::rcat(1, p_weights)
    newans <- rep(NA, K)
    
    # Store original proposal probs for Metropolis-Hastings correction
    orig_p_prob <- state$p[newgid]
    new_p_prob <- p_weights[newgid]
    
    # Log of proposal ratio correction
    proposal_correction <- log(orig_p_prob) - log(new_p_prob)
    
    for(k in 1 : K){  # For each feature
      # Create biased proposal for answers
      ans_weights <- state$pijk[[k]][,newgid]
      for(j in 1:Jk[k]) {
        diff <- mijk[[k]][j,newgid] - state$nijk[[k]][j,newgid]
        ans_weights[j] <- ans_weights[j] * exp(beta * diff)
      }
      ans_weights <- ans_weights / sum(ans_weights)
      
      # Sample from biased proposal
      newans[k] <- extraDistr::rcat(1, ans_weights)
      
      # Correct for modified proposal in MH ratio
      orig_ans_prob <- state$pijk[[k]][newans[k],newgid]
      new_ans_prob <- ans_weights[newans[k]]
      proposal_correction <- proposal_correction + log(orig_ans_prob) - log(new_ans_prob)
      
      # Standard likelihood ratio calculation
      lratio <- lratio + extraDistr::dlaplace(
        state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid] + 1, 
        0, 2 * K / epsilon, TRUE
      )
      lratio <- lratio + extraDistr::dlaplace(
        state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]] - 1, 
        0, 2 * K / epsilon, TRUE
      )
      lratio <- lratio - extraDistr::dlaplace(
        state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid], 
        0, 2 * K / epsilon, TRUE
      )
      lratio <- lratio - extraDistr::dlaplace(
        state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]], 
        0, 2 * K / epsilon, TRUE
      )
    }
    
    # Add proposal ratio correction to maintain detailed balance
    lratio <- lratio + proposal_correction
    
    if(log(runif(1)) < lratio){  # If accept then book keeping
      for(k in 1 : K){
        state$nijk[[k]][state$ans[id,k], state$g[id]] <- state$nijk[[k]][state$ans[id,k], state$g[id]] - 1
        state$nijk[[k]][newans[k], newgid] <- state$nijk[[k]][newans[k], newgid] + 1
      }
      state$n[state$g[id]] <- state$n[state$g[id]] - 1
      state$n[newgid] <- state$n[newgid] + 1
      state$g[id] <- newgid
      state$ans[id,] <- newans
      acnt <- acnt + 1  # Increment acceptance count
    }
    lratio <- 0
  }
  
  state$accept <- acnt / N
  state$dis <- get_sse(state$nijk, mijk)
  return(state)
}

# Function to run experiments across epsilon values
run_comparison <- function(epsvec, niter, beta) {
  results <- list()
  
  for(eps_idx in 1:length(epsvec)) {
    epsilon <- epsvec[eps_idx]
    cat("Running with epsilon =", epsilon, "\n")
    
    # Add noise to create privatized data
    mijk <- list()
    for(k in 1:K) {
      mijk[[k]] <- nijk[[k]] + VGAM::rlaplace(n = Jk[k] * I, location = 0, scale = 2 * K / epsilon)
    }
    
    # Initialize both samplers from the same starting point
    init_state <- rgibbs_init()
    state_std <- init_state
    state_sdp <- init_state
    
    # Storage for results
    p_gibbs_std <- matrix(NA, ncol = I, nrow = niter)
    p_gibbs_sdp <- matrix(NA, ncol = I, nrow = niter)
    acc_gibbs_std <- rep(NA, niter)
    acc_gibbs_sdp <- rep(NA, niter)
    dis_gibbs_std <- rep(NA, niter)
    dis_gibbs_sdp <- rep(NA, niter)
    
    # Run both samplers
    for(iter in 1:niter) {
      # Standard sampler
      state_std <- gibbs_iter_standard(state_std, mijk, epsilon)
      p_gibbs_std[iter,] <- state_std$p
      acc_gibbs_std[iter] <- state_std$accept
      dis_gibbs_std[iter] <- state_std$dis
      
      # s_dp-aware sampler
      state_sdp <- gibbs_iter_sdp_aware(state_sdp, mijk, epsilon, beta)
      p_gibbs_sdp[iter,] <- state_sdp$p
      acc_gibbs_sdp[iter] <- state_sdp$accept
      dis_gibbs_sdp[iter] <- state_sdp$dis
      
      if(iter %% 500 == 0) {
        cat("  Iteration", iter, "of", niter, "completed\n")
        cat("  Standard acceptance rate:", mean(acc_gibbs_std[max(1, iter-100):iter]), "\n")
        cat("  s_dp-aware acceptance rate:", mean(acc_gibbs_sdp[max(1, iter-100):iter]), "\n")
      }
    }
    
    # Calculate effective sample sizes
    ess_std <- effectiveSize(mcmc(p_gibbs_std[-(1:500),]))
    ess_sdp <- effectiveSize(mcmc(p_gibbs_sdp[-(1:500),]))
    
    # Store results
    results[[eps_idx]] <- list(
      epsilon = epsilon,
      p_std = p_gibbs_std,
      p_sdp = p_gibbs_sdp,
      acc_std = acc_gibbs_std,
      acc_sdp = acc_gibbs_sdp,
      dis_std = dis_gibbs_std,
      dis_sdp = dis_gibbs_sdp,
      ess_std = ess_std,
      ess_sdp = ess_sdp,
      mean_acc_std = mean(acc_gibbs_std[-(1:500)]),
      mean_acc_sdp = mean(acc_gibbs_sdp[-(1:500)])
    )
  }
  
  return(results)
}

# Run the comparison
results <- run_comparison(epsvec, niter, beta_param)

# Extract acceptance rates
acc_data <- data.frame(
  epsilon = numeric(),
  method = character(),
  acceptance_rate = numeric()
)

for(i in 1:length(results)) {
  acc_data <- rbind(acc_data, 
                    data.frame(
                      epsilon = results[[i]]$epsilon,
                      method = "Standard",
                      acceptance_rate = results[[i]]$mean_acc_std
                    ))
  acc_data <- rbind(acc_data, 
                    data.frame(
                      epsilon = results[[i]]$epsilon,
                      method = "s_dp-aware",
                      acceptance_rate = results[[i]]$mean_acc_sdp
                    ))
}

# Extract ESS improvement
ess_improvement <- data.frame(
  epsilon = sapply(results, function(r) r$epsilon),
  improvement = sapply(results, function(r) mean(r$ess_sdp)/mean(r$ess_std))
)

# Generate acceptance rate plot
acc_plot <- ggplot(acc_data, aes(x=factor(epsilon), y=acceptance_rate, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) +
  labs(x="Privacy budget (ε)", y="Acceptance rate", 
       title="Acceptance rates comparison across privacy levels") +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  theme(legend.title = element_blank())

# Save p lot
ggsave("acceptance_rates_comparison.png", acc_plot, width=10, height=6)

# Generate trace plots for epsilon = 1
eps1_idx <- which(epsvec == 1)
if(length(eps1_idx) > 0) {
  trace_data_std <- data.frame(
    iteration = 1:niter,
    p1 = results[[eps1_idx]]$p_std[,1]
  )
  
  trace_data_sdp <- data.frame(
    iteration = 1:niter,
    p1 = results[[eps1_idx]]$p_sdp[,1]
  )
  
  trace_plot_std <- ggplot(trace_data_std, aes(x=iteration, y=p1)) +
    geom_line() +
    ylim(0, 0.5) +
    labs(x="Iteration", y="p1 value", 
         title="Standard proposals") +
    theme_minimal()
  
  trace_plot_sdp <- ggplot(trace_data_sdp, aes(x=iteration, y=p1)) +
    geom_line() +
    ylim(0, 0.5) +
    labs(x="Iteration", y="p1 value", 
         title="s_dp-aware proposals") +
    theme_minimal()
  
  # Combine plots
  combined_trace <- grid.arrange(trace_plot_std, trace_plot_sdp, ncol=2)
  
  # Save combined plot
  ggsave("trace_plot_comparison.png", combined_trace, width=12, height=5)
}

# Print summary
cat("\nResults Summary:\n")
for(i in 1:length(results)) {
  cat("Epsilon =", results[[i]]$epsilon, "\n")
  cat("  Standard acceptance rate:", round(results[[i]]$mean_acc_std, 4), "\n")
  cat("  s_dp-aware acceptance rate:", round(results[[i]]$mean_acc_sdp, 4), "\n")
  cat("  Improvement factor:", round(results[[i]]$mean_acc_sdp/results[[i]]$mean_acc_std, 2), "×\n")
  cat("  ESS improvement:", round(mean(results[[i]]$ess_sdp)/mean(results[[i]]$ess_std), 2), "×\n\n")
}

# Save results
save(results, file="sdp_aware_proposal_results.RData")
write.csv(acc_data, "acceptance_rates.csv", row.names=FALSE)
write.csv(ess_improvement, "ess_improvement.csv", row.names=FALSE)


# 1. Create an improved version of the acceptance rate plot (more compact)
acc_plot <- ggplot(acc_data, aes(x=factor(epsilon), y=acceptance_rate, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.6) +
  labs(x="Privacy budget (ε)", y="Acceptance rate") +
  theme_minimal() +
  scale_fill_brewer(palette="Set1") +
  theme(legend.position="top", 
        legend.title = element_blank(),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

# 2. Create ESS improvement plot
ess_plot <- ggplot(ess_improvement, aes(x=factor(epsilon), y=improvement)) +
  geom_bar(stat="identity", fill="#4575b4") +
  geom_text(aes(label=round(improvement, 2)), vjust=-0.5, size=3.5) +
  labs(x="Privacy budget (ε)", y="ESS improvement factor") +
  theme_minimal() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

# Save the acceptance rate plot
ggsave("acceptance_rate_plot.png", acc_plot, width=6, height=4, dpi=300)

# Save the ESS improvement plot
ggsave("ess_improvement_plot.png", ess_plot, width=6, height=4, dpi=300)
