#!/bin/bash
EPS_VALUES=(0.05 0.1 0.3 0.5 1 2 3 5 10 20)
NITER=2000  # Reduced iterations

for eps in "${EPS_VALUES[@]}"; do
  # Original method
  Rscript tradeoff_gibbs.R --method original --epsilon $eps --niter $NITER
  
  # Smart method
  Rscript tradeoff_gibbs.R --method smart --epsilon $eps --niter $NITER --beta 0.5
  
  # Adaptive method
  Rscript adaptive_gibbs.R --epsilon $eps --niter $NITER
done

