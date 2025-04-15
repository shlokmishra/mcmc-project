# Run key comparisons (ε=0.1,1,10)
for eps in 0.1 1 10; do
  # Original
  Rscript smart_gibbs.R --method original --epsilon $eps
  # Smart
  Rscript smart_gibbs.R --method smart --epsilon $eps
  # Hybrid 
  Rscript hybrid_gibbs.R --epsilon $eps
done

