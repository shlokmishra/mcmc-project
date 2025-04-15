# hybrid_gibbs.R
rm(list=ls())
library(extraDistr)
library(coda)

# ----------------------------
# Configuration & Arguments
# ----------------------------
niter <- 10000
epsilon <- 1.0
method <- "hybrid"
beta <- 0.5  # used in smart phase

args <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--epsilon") {
    epsilon <- as.numeric(args[i + 1])
  } else if (args[i] == "--beta") {
    beta <- as.numeric(args[i + 1])
  }
}
cat(sprintf("Running hybrid Gibbs with epsilon = %.2f and beta = %.2f\n", epsilon, beta))

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
  state$start_time <- start_time
  state
}

# ----------------------------
# Hybrid Gibbs Iteration
# ----------------------------
hybrid_gibbs_iter <- function(state, iter) {
  # Parameter updates
  state$p <- rdirichlet(1, state$n + 2)[1,]
  for(k in 1:K) {
    state$pijk[[k]] <- t(sapply(1:I, function(i) 
      rdirichlet(1, state$nijk[[k]][,i] + 2)))
  }
  
  # Proposal phase selection
  acnt <- 0
  use_smart <- iter >= 0.3 * niter
  if (use_smart) {
    w <- compute_weights(state$nijk, beta)
  }
  
  for(id in 1:N) {
    if (use_smart) {
      g_weights <- state$p * sapply(1:I, function(i) 
        mean(sapply(w, function(m) m[i, state$ans[id,]])))
      newgid <- rcat(1, g_weights)
      newans <- sapply(1:K, function(k) rcat(1, w[[k]][newgid,]))
    } else {
      newgid <- rcat(1, state$p)
      newans <- sapply(1:K, function(k) rcat(1, state$pijk[[k]][newgid,]))
    }
    
    # MH acceptance ratio
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
  
  state$accept <- acnt / N
  state
}

# ----------------------------
# Run Sampler
# ----------------------------
state <- rgibbs_init()
p_gibbs <- matrix(NA, niter, I)
acc_gibbs <- numeric(niter)

for(iter in 1:niter) {
  state <- hybrid_gibbs_iter(state, iter)
  p_gibbs[iter,] <- state$p
  acc_gibbs[iter] <- state$accept
  
  if(iter %% 100 == 0) {
    cat(sprintf("[hybrid] Iter %d/%d | Accept: %.2f\n", iter, niter, state$accept))
  }
}

# ----------------------------
# Save Results
# ----------------------------
results <- list(
  p_samples = p_gibbs,
  acc_rates = acc_gibbs,
  epsilon = epsilon,
  method = "hybrid",
  beta = beta,
  runtime = as.numeric(difftime(Sys.time(), state$start_time, units = "secs"))
)

filename <- sprintf("hybrid_eps%.2f.rds", epsilon)
saveRDS(results, filename)
cat(sprintf("Saved results to %s\n", filename))
