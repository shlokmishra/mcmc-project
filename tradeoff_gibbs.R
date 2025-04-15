# tradeoff_gibbs.R
rm(list=ls())
library(extraDistr)
library(coda)

# ----------------------------
# Configuration and Argument Parsing
# ----------------------------
niter <- 2000
method <- "original"
epsilon <- 0.1
beta <- 0.5

args <- commandArgs(trailingOnly = TRUE)
for (i in seq(1, length(args), by = 2)) {
  if (args[i] == "--method") {
    method <- args[i + 1]
  } else if (args[i] == "--epsilon") {
    epsilon <- as.numeric(args[i + 1])
  } else if (args[i] == "--beta") {
    beta <- as.numeric(args[i + 1])
  }
}
cat(sprintf("Running with method = %s, epsilon = %.2f, beta = %.2f\n", method, epsilon, beta))

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
# Gibbs Initialization with Timing
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



gibbs_iter <- function(state) {
  # ----------------------------
  # Runtime Safeguard
  # ----------------------------
  if (any(sapply(state$nijk, function(x) any(x < 0)))) {
    warning("Negative counts detected, resetting...")
    state$nijk <- lapply(state$nijk, function(x) pmax(x, 0))
  }
  
  # ----------------------------
  # Parameter updates
  # ----------------------------
  state$p <- rdirichlet(1, state$n + 2)[1,]
  for (k in 1:K) {
    state$pijk[[k]] <- t(sapply(1:I, function(i)
      rdirichlet(1, state$nijk[[k]][, i] + 2)))
  }
  
  # ----------------------------
  # Proposal weight computation
  # ----------------------------
  acnt <- 0
  if (method == "smart") {
    w <- compute_weights(state$nijk, beta)
  }
  
  # ----------------------------
  # Data augmentation loop
  # ----------------------------
  for (id in 1:N) {
    if (method == "smart") {
      g_weights <- state$p * sapply(1:I, function(i)
        mean(sapply(w, function(m) m[i, state$ans[id, ]])))
      newgid <- rcat(1, g_weights)
      newans <- sapply(1:K, function(k) rcat(1, w[[k]][newgid, ]))
    } else {
      newgid <- rcat(1, state$p)
      newans <- sapply(1:K, function(k) rcat(1, state$pijk[[k]][newgid, ]))
    }
    
    # ----------------------------
    # Acceptance ratio with clamping
    # ----------------------------
    lratio <- sum(sapply(1:K, function(k) {
      term1 <- dlaplace(
        pmax(state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid] + 1, 0),
        0, 2*K/epsilon, TRUE
      )
      term2 <- dlaplace(
        pmax(state$nijk[[k]][state$ans[id, k], state$g[id]] - mijk[[k]][state$ans[id, k], state$g[id]] - 1, 0),
        0, 2*K/epsilon, TRUE
      )
      term3 <- dlaplace(
        pmax(state$nijk[[k]][newans[k], newgid] - mijk[[k]][newans[k], newgid], 0),
        0, 2*K/epsilon, TRUE
      )
      term4 <- dlaplace(
        pmax(state$nijk[[k]][state$ans[id, k], state$g[id]] - mijk[[k]][state$ans[id, k], state$g[id]], 0),
        0, 2*K/epsilon, TRUE
      )
      term1 + term2 - term3 - term4
    }))
    
    if (log(runif(1)) < lratio) {
      for (k in 1:K) {
        state$nijk[[k]][state$ans[id, k], state$g[id]] <- state$nijk[[k]][state$ans[id, k], state$g[id]] - 1
        state$nijk[[k]][newans[k], newgid] <- state$nijk[[k]][newans[k], newgid] + 1
      }
      state$n[state$g[id]] <- state$n[state$g[id]] - 1
      state$n[newgid] <- state$n[newgid] + 1
      state$g[id] <- newgid
      state$ans[id, ] <- newans
      acnt <- acnt + 1
    }
  }
  
  # ----------------------------
  # Save iteration statistics
  # ----------------------------
  state$accept <- acnt / N
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
# Save Results with Runtime
# ----------------------------
results <- list(
  p_samples = p_gibbs,
  acc_rates = acc_gibbs,
  jump_sizes = jump_sizes,
  epsilon = epsilon,
  method = method,
  beta = beta,
  runtime = as.numeric(difftime(Sys.time(), state$start_time, units = "secs"))
)

filename <- sprintf("tradeoff_eps%.2f_%s.rds", epsilon, method)
saveRDS(results, filename)
cat(sprintf("Saved results to %s\n", filename))
