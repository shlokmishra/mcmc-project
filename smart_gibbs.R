rm(list = ls())
library(extraDistr)
library(coda)
library(ggplot2)
library(gridExtra)

# ----------------------------
# Configuration
# ----------------------------
iexp <- 1
noiselevel <- 1
niter <- 10000
beta <- 0.5  # Weight parameter
method <- "smart"  # "original" or "smart"
args <- commandArgs(trailingOnly = TRUE)

# Default values
method <- "smart"
noiselevel <- 1

# Custom parsing
for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--method") {
    method <- args[i + 1]
  } else if (args[i] == "--noiselevel") {
    noiselevel <- as.numeric(args[i + 1])
  } else if (args[i] == "--beta") {
    beta <- as.numeric(args[i + 1])
  }
}

# ----------------------------
# Data Loading
# ----------------------------
set.seed(iexp)
epsvec <- rev(c(0.1, 0.3, 1, 3, 10))
epsilon <- epsvec[noiselevel]
load("loglinear_datageneratingparameters.RData")
load("loglinear_notdpdata_N100.RData")

# ----------------------------
# DP Noise Addition
# ----------------------------
mijk <- list()
for(k in 1:K) {
  mijk[[k]] <- nijk[[k]] + VGAM::rlaplace(Jk[k] * I, 0, 2*K/epsilon)
}

# ----------------------------
# Helper Functions
# ----------------------------
get_n <- function(g) as.vector(table(factor(g, levels=1:I)))

get_nijk <- function(g, ans) {
  lapply(1:K, function(k) 
    as.matrix(table(factor(ans[,k], levels=1:Jk[k]), factor(g, levels=1:I))))
}

get_sse <- function(nijk) sum(sapply(1:K, function(k) sum((nijk[[k]] - mijk[[k]])^2)))

compute_weights <- function(nijk, beta) {
  lapply(1:K, function(k) {
    w <- exp(beta * (mijk[[k]] - nijk[[k]]))
    t(apply(w, 2, function(x) x/sum(x)))
  })
}

# ----------------------------
# Gibbs Initialization
# ----------------------------
rgibbs_init <- function() {
  state <- list()
  state$p <- rdirichlet(1, rep(2, I))[1,]
  state$pijk <- lapply(1:K, function(k) 
    t(sapply(1:I, function(i) rdirichlet(1, rep(2, Jk[k])))))
  state$g <- sample(1:I, N, TRUE, state$p)
  state$ans <- sapply(1:K, function(k) 
    sapply(state$g, function(i) rcat(1, state$pijk[[k]][i,])))
  state$n <- get_n(state$g)
  state$nijk <- get_nijk(state$g, state$ans)
  state
}

# ----------------------------
# Gibbs Iteration (Modified)
# ----------------------------
gibbs_iter <- function(state) {
  # Parameter updates
  state$p <- rdirichlet(1, state$n + 2)[1,]
  for(k in 1:K) {
    state$pijk[[k]] <- t(sapply(1:I, function(i) 
      rdirichlet(1, state$nijk[[k]][,i] + 2)))
  }
  
  # Data augmentation with smart proposals
  acnt <- 0
  if(method == "smart") {
    w <- compute_weights(state$nijk, beta)
  }
  
  for(id in 1:N) {
    if(method == "smart") {
      # Weighted proposal
      g_weights <- state$p * sapply(1:I, function(i) 
        mean(sapply(w, function(m) m[i, state$ans[id,]])))
      newgid <- rcat(1, g_weights)
      newans <- sapply(1:K, function(k) rcat(1, w[[k]][newgid,]))
    } else {
      # Original proposal
      newgid <- rcat(1, state$p)
      newans <- sapply(1:K, function(k) rcat(1, state$pijk[[k]][newgid,]))
    }
    
    # Acceptance ratio
    lratio <- sum(sapply(1:K, function(k) {
      dlaplace(state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid] + 1, 0, 2*K/epsilon, TRUE) +
        dlaplace(state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]] - 1, 0, 2*K/epsilon, TRUE) -
        dlaplace(state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid], 0, 2*K/epsilon, TRUE) -
        dlaplace(state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]], 0, 2*K/epsilon, TRUE)
    }))
    
    if(log(runif(1)) < lratio) {
      for(k in 1:K) {
        state$nijk[[k]][state$ans[id,k], state$g[id]] <- state$nijk[[k]][state$ans[id,k], state$g[id]] - 1
        state$nijk[[k]][newans[k], newgid] <- state$nijk[[k]][newans[k], newgid] + 1
      }
      state$n[state$g[id]] <- state$n[state$g[id]] - 1
      state$n[newgid] <- state$n[newgid] + 1
      state$g[id] <- newgid
      state$ans[id,] <- newans
      acnt <- acnt + 1
    }
  }
  
  state$accept <- acnt/N
  state$dis <- get_sse(state$nijk)
  state
}

# ----------------------------
# Run Sampler
# ----------------------------
state <- rgibbs_init()
p_gibbs <- matrix(NA, niter, I)
acc_gibbs <- numeric(niter)
jump_sizes <- numeric(niter)

for(iter in 1:niter) {
  state <- gibbs_iter(state)
  p_gibbs[iter,] <- state$p
  acc_gibbs[iter] <- state$accept
  jump_sizes[iter] <- if(iter > 1) sqrt(sum((p_gibbs[iter,] - p_gibbs[iter-1,])^2)) else 0
  
  if(iter %% 100 == 0) {
    cat(sprintf("[%s] Iter %d/%d | Accept: %.2f | Jump: %.3f\n",
                method, iter, niter, state$accept, jump_sizes[iter]))
  }
}

# ----------------------------
# Save Results
# ----------------------------
results <- list(
  p_samples = p_gibbs,
  acc_rates = acc_gibbs,
  jump_sizes = jump_sizes,
  epsilon = epsilon,
  method = method,
  beta = beta
)

save_filename <- sprintf("results_eps%s_%s.rds",
                         formatC(epsilon, format="f", digits=1),
                         method)
saveRDS(results, save_filename)
# ----------------------------
# Generate Figures
# ----------------------------
if(interactive()) {
  # Trace Plots
  df_trace <- data.frame(
    Iteration = 1:niter,
    p1 = p_gibbs[,1],
    p2 = p_gibbs[,2],
    p3 = p_gibbs[,3],
    p4 = p_gibbs[,4],
    p5 = p_gibbs[,5]
  )
  
  p_trace <- ggplot(df_trace, aes(x=Iteration)) +
    geom_line(aes(y=p1, color="p1")) +
    geom_line(aes(y=p2, color="p2")) +
    geom_line(aes(y=p3, color="p3")) +
    geom_line(aes(y=p4, color="p4")) +
    geom_line(aes(y=p5, color="p5")) +
    labs(title=sprintf("Trace Plot (ε=%.1f, %s)", epsilon, method),
         y="Parameter Value", color="Parameter")
  
  # Acceptance Rate Plot
  df_acc <- data.frame(
    Iteration = 1:niter,
    Acceptance = acc_gibbs
  )
  
  p_acc <- ggplot(df_acc, aes(x=Iteration, y=Acceptance)) +
    geom_line() +
    geom_hline(yintercept=exp(-epsilon), linetype="dashed", color="red") +
    labs(title=sprintf("Acceptance Rate (ε=%.1f, %s)", epsilon, method))
  
  # Combine Plots
  grid.arrange(p_trace, p_acc, ncol=1)
  
  # Effective Sample Size
  ess <- effectiveSize(p_gibbs)
  cat(sprintf("\nEffective Sample Sizes (ε=%.1f, %s):\n", epsilon, method))
  print(ess)
}
