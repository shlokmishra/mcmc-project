#!/bin/bash

for noiselevel in 1 2 3 4 5; do
  echo "Running original method with noiselevel $noiselevel"
  Rscript smart_gibbs.R --method original --noiselevel $noiselevel

  echo "Running smart method with noiselevel $noiselevel"
  Rscript smart_gibbs.R --method smart --noiselevel $noiselevel --beta 0.5
done

