# NIRApost 1.2.2 (2026-08-31)

Faster moderation-effect testing, and the bridge (two-community) NIRA workflow
now includes its prerequisite moderation test.

### runMgmmAnalysis: parallel resampling (`ncores`)

- New `ncores` argument parallelises the `p` candidate moderator nodes with
  `parallel::parLapply`. Because `mgm::resample()` has no internal
  parallelisation, this gives an approximately `ncores`-fold speed-up on
  multi-core machines (about 4x on 6 cores for the bridge dataset), while
  producing results **identical** to the serial run.
- Serial execution (`ncores = 1`) remains the default and is unchanged.

### Bridge NIRA workflow: moderation test added as prerequisite

- `Bridge_network_NIRA.R` now runs `runMgmmAnalysis()` (Step 3) on the whole
  network before the bridge intervention, consistent with the other
  workflows. An empty `significant_moderators` confirms the fixed-edge-weight
  assumption holds.
- README and the package tutorial (Workflow 3) document the interpretation
  and reporting of bridge-NIRA results.

---

# NIRApost 1.2.1 (2026-08-29)

API refinement of the analytical workflow:

- `analyticalNIRAtest` no longer has a `resample_type` argument; it always
  returns both the bootstrap test (`$stat`) and the random-target test
  (`$random`).
- `permutationNIRAtest` gains a `$random` element (rank-based target
  comparison, uncorrected); its `$stat` and `$plot_data` are unchanged.
- Fixed an inconsistency where the bootstrap resampling estimator could
  differ from the point-estimate estimator. New `resample_default` argument
  in both `analyticalNIRAtest()` and `analyticalBridgeNIRAtest()`:
  by default (`NULL`) the resampling uses the same estimator as the point
  estimate (bootnet when parameters are supplied, the fast pseudo-likelihood
  estimator otherwise), so confidence intervals are correctly centred on the
  point estimate. Set `resample_default = "fast"` to force the fast
  estimator.

---

# NIRApost 1.2.0 (2026-08-29)

Adds an exact analytical workflow alongside the original Monte-Carlo
pipeline, and improves how uncertainty is quantified.

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
- **`permutationNIRAtest()` returns `$random`** — a new element comparing each
  node's simulated effect with a randomly chosen target node (rank-based
  p-value, uncorrected), mirroring the analytical workflow.

### Other improvements

- Updated package vignette documenting all three workflows and when to choose
  each.
- New tests: exact-value verification against independent brute-force
  enumeration, consistency between the two analytical functions, and
  regression tests for the resampling estimator.
- No changes to the original V1.1.0 Monte-Carlo functions — all existing code
  remains backward compatible.

> **Note:** the analytical workflow is intended for small-to-moderate networks
> (p ≲ 20). For larger networks where 2^p enumeration is infeasible, the
> V1.1.0 Monte-Carlo workflow remains the recommended choice.
