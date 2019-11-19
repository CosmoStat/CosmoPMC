
Examples to test and run tempered PMC.

1. `original/1_mvnorm_centered_2D`

   No tempering, to test Bayesian evidence and approximations.
   Run

   ``` $COSMOPMC/bin/cosmo_pmc.pl -P $COSMOPMC -f "0.5 0.5"```

   The likelihood is a 2D multi-variate normal. The prior volume is
   the unit square. Since the likelihood is normalised, the Bayesian evidence
   is the inverse prior volume, which is 1, computed and written to the
   file `evidence_analytic`..

   The file `evidence` should be consistent with 1.

   The Laplace approximation with the Fisher matrix replacing the input likelihood
   uses as maximum-likelihood a point very close to the
   true value (given as starting point for the maximum-search by `-f "0.5 0.5")`.
   The file `evidence_fisher` should thus also be consistent with 1.

   Another Laplace approximation is computed using a the sample covariance
   instead of the true input covariance matrix. The results are performed for each
   iteration, and written to `iter_?/evidence_covariance`. For the final iteration
   the result should be close to 1.

