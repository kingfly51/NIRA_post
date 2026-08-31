# NIRApost 1.2.0 (2026-08-31)

Adds an exact analytical workflow alongside the original Monte-Carlo
pipeline, improves how uncertainty is quantified, and accelerates the
moderation-effect prerequisite test.

### Highlights

- **Analytical NIRA (`analyticalNIRAtest()`)** — computes the expected number
  of activated nodes *exactly* by enumerating all 2^p Ising microstates with a
  numerically stable log-sum-exp scheme. This eliminates Monte Carlo error
  entirely and runs in milliseconds for small-to-moderate networks (p ≲ 20),
  reproducing Isinglandr's exact results to machine precision.
- **Bridge / two-community NIRA (`analyticalBridgeNIRAtest()`)** — new exact
  analysis for networks partitioned into two communities: intervene on each
  node of a *source* community and obtain the exact mean activation of the
  *outcome* community, with the same resampling-based inference.
- **Uncertainty quantified by resampling the *original* data** — following Cui
  et al. (2026), the remaining (sampling) error is characterized by
  bootstrapping participants and re-estimating the network on each resample,
  instead of relying on simulated samples. Each function returns a bootstrap
  test (`$stat`: SE, CI, p for H0: Δ = 0, adjusted p, `percenttop_1` ranking
  stability) and a random-target permutation test (`$random`: is node *i* a
  better target than a randomly chosen node?).
- **Consistent resampling estimator (`resample_default`)** — the bootstrap
  re-estimation now uses the same estimator as the point estimate by default
  (bootnet when network parameters are supplied, the fast pseudo-likelihood
  estimator otherwise), so confidence intervals are correctly centred on the
  point estimate. Available in both `analyticalNIRAtest()` and
  `analyticalBridgeNIRAtest()`.
- **`permutationNIRAtest()` returns `$random`** — a new element comparing each
  node's simulated effect with a randomly chosen target node (rank-based
  p-value, uncorrected), mirroring the analytical workflow.
- **Faster moderation-effect testing (`runMgmmAnalysis`)** — new `ncores`
  argument parallelises the candidate moderator nodes, giving an
  approximately `ncores`-fold speed-up on multi-core machines while producing
  results identical to the serial run. The bridge (two-community) NIRA
  workflow now runs its whole-network moderation prerequisite test in
  `Bridge_network_NIRA.R`.

### Other improvements

- Updated README and package vignette documenting all workflows (Monte Carlo,
  analytical, and bridge NIRA), including how to interpret and report
  bridge-NIRA results.
- New tests: exact-value verification against independent brute-force
  enumeration, consistency between the two analytical functions, and
  regression tests for the resampling estimator.
- No changes to the original V1.1.0 Monte-Carlo functions — all existing code
  remains backward compatible.

> **Note:** the analytical workflow is intended for small-to-moderate networks
> (p ≲ 20). For larger networks where 2^p enumeration is infeasible, the
> V1.1.0 Monte-Carlo workflow remains the recommended choice.
