# adaptive_gibbs.R
rm(list=ls())
library(extraDistr)
library(coda)

# ----------------------------
# Configuration
# ----------------------------
niter <- 2000
epsilon <- 1.0
beta <- 0.5  # Initial beta
method <- "adaptive"

args <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--epsilon") {
    epsilon <- as.numeric(args[i + 1])
  } else if (args[i] == "--beta") {
    beta <- as.numeric(args[i + 1])
  }
}
cat(sprintf("Running adaptive Gibbs with epsilon = %.2f and initial beta = %.2f\n", epsilon, beta))

# ----------------------------
# Data Loading
# ----------------------------
load("loglinear_datageneratingparameters.RData")
load("loglinear_notdpdata_N100.RData")
set.seed(1)

# DP noise addition
mijk <- list()
for(k in 1:K) {
  mijk[[k]] <- nijk[[k]] + VGAM::rlaplace(Jk[k] * I, 0, 2*K/epsilon)
}

# ----------------------------
# Helper Functions
# ----------------------------
clamp <- function(x, min, max) pmin(pmax(x, min), max)

get_n <- function(g) as.vector(table(factor(g, levels=1:I)))

get_nijk <- function(g, ans) {
  lapply(1:K, function(k) 
    as.matrix(table(factor(ans[,k], levels=1:Jk[k]), factor(g, levels=1:I))))
}

get_sse <- function(nijk) sum(sapply(1:K, function(k) sum((nijk[[k]] - mijk[[k]])^2)))

compute_weights <- function(nijk, beta) {
  lapply(1:K, function(k) {
    w <- exp(beta * (mijk[[k]] - nijk[[k]]))
    t(apply(w, 2, function(x) x / sum(x)))
  })
}

# ----------------------------
# Initialization
# ----------------------------
rgibbs_init <- function() {
  start_time <- Sys.time()
  state <- list()
  state$p <- rdirichlet(1, rep(2, I))[1,]
  state$pijk <- lapply(1:K, function(k) 
    t(sapply(1:I, function(i) rdirichlet(1, rep(2, Jk[k])))))
  state$g <- sample(1:I, N, TRUE, state$p)
  state$ans <- sapply(1:K, function(k) 
    sapply(state$g, function(i) rcat(1, state$pijk[[k]][i,])))
  state$n <- get_n(state$g)
  state$nijk <- get_nijk(state$g, state$ans)
  state$beta <- beta
  state$acc_rates <- numeric()
  state$start_time <- start_time
  state
}

# ----------------------------
# Adaptive Gibbs Iteration
# ----------------------------
adaptive_gibbs_iter <- function(state, iter) {
  # Update parameters
  state$p <- rdirichlet(1, state$n + 2)[1,]
  for(k in 1:K) {
    state$pijk[[k]] <- t(sapply(1:I, function(i) 
      rdirichlet(1, state$nijk[[k]][,i] + 2)))
  }
  
  # Modified beta adaptation with safeguards
  if (iter > 100 && iter %% 50 == 0) {
    recent_accept <- mean(tail(na.omit(state$acc_rates), 50))
    state$beta <- clamp(
      ifelse(recent_accept < 0.25, state$beta * 0.95, state$beta * 1.05),
      0.01, 2.0
    )
  }
  
  # Compute and normalize weights
  w_raw <- compute_weights(state$nijk, state$beta)
  w <- lapply(w_raw, function(mat) {
    mat <- mat - max(mat)  # Stabilize
    exp_mat <- exp(mat)
    exp_mat / sum(exp_mat)
  })
  
  # Data augmentation step
  acnt <- 0
  for(id in 1:N) {
    g_weights <- state$p * sapply(1:I, function(i) 
      mean(sapply(w, function(m) m[i, state$ans[id,]])))
    newgid <- rcat(1, g_weights)
    newans <- sapply(1:K, function(k) rcat(1, w[[k]][newgid,]))
    
    # MH acceptance
    lratio <- sum(sapply(1:K, function(k) {
      dlaplace(state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid] + 1, 0, 2*K/epsilon, TRUE) +
        dlaplace(state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]] - 1, 0, 2*K/epsilon, TRUE) -
        dlaplace(state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid], 0, 2*K/epsilon, TRUE) -
        dlaplace(state$nijk[[k]][state$ans[id,k], state$g[id]] - mijk[[k]][state$ans[id,k], state$g[id]], 0, 2*K/epsilon, TRUE)
    }))
    
    if (log(runif(1)) < lratio) {
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
  
  # Record acceptance
  accept_rate <- acnt / N
  state$acc_rates[iter] <- accept_rate
  state$accept <- accept_rate
  state
}

# ----------------------------
# Run Adaptive Gibbs Sampler
# ----------------------------
state <- rgibbs_init()
p_gibbs <- matrix(NA, niter, I)
acc_gibbs <- numeric(niter)
beta_trace <- numeric(niter)

for(iter in 1:niter) {
  state <- adaptive_gibbs_iter(state, iter)
  p_gibbs[iter,] <- state$p
  acc_gibbs[iter] <- state$accept
  beta_trace[iter] <- state$beta
  
  if(iter %% 100 == 0) {
    cat(sprintf("[adaptive] Iter %d/%d | Accept: %.2f | Beta: %.3f\n",
                iter, niter, state$accept, state$beta))
  }
}

# ----------------------------
# Save Results
# ----------------------------
results <- list(
  p_samples = p_gibbs,
  acc_rates = acc_gibbs,
  beta_trace = beta_trace,
  epsilon = epsilon,
  method = "adaptive",
  runtime = as.numeric(difftime(Sys.time(), state$start_time, units = "secs"))
)

filename <- sprintf("adaptive_eps%.2f.rds", epsilon)
saveRDS(results, filename)
cat(sprintf("Saved results to %s\n", filename))
