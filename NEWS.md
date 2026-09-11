# What’s New in v1.2.0

NIRApost v1.2.0 introduces exact analytical intervention analysis alongside the Monte Carlo workflow described in the original AMPPS tutorial. The update supports both whole-network and community-specific outcomes, quantifies sampling uncertainty through bootstrap resampling, and adds parallel moderation screening.

## Highlights

### Exact analytical NIRA (`analyticalNIRAtest()`)

Computes the expected number of activated nodes by enumerating all 2^p Ising microstates. The calculation eliminates Monte Carlo approximation error in intervention expectations and uses stable log-weight normalization. Exact enumeration is particularly useful for smaller networks, while the original Monte Carlo workflow remains available when enumeration becomes computationally impractical.

### Bridge and community-specific analysis (`analyticalBridgeNIRAtest()`)

Extends the analytical workflow to networks with defined communities. The function perturbs each node in a source community and calculates the exact expected activation count in an outcome community, allowing researchers to examine model-implied intervention effects across communities.

### Sampling uncertainty assessed from the original data

Both analytical functions bootstrap participants, re-estimate the network, and recompute intervention outcomes. Their `$stat` output includes intervention effects, Glass’s Δ, bootstrap standard errors, 95% percentile confidence intervals, approximate centered-bootstrap p-values for H₀: Δᵢ = 0, and multiplicity-adjusted p-values. The `percenttop_1` statistic summarizes how frequently each node has the largest absolute effect across bootstrap samples, providing a measure of ranking stability. Tied first places receive equal shares of the first-place credit.

### Configurable estimation and parallel resampling

The `resample_default` argument controls the network estimator used during bootstrap resampling. The default uses IsingFit when network parameters are supplied and the fast pseudo-likelihood estimator otherwise. Users should ensure that bootstrap estimation matches their original fitting procedure. The `parallel`, `ncores`, and `seed` arguments support parallel execution and reproducibility.

### Integrated analytical visualization

Each analytical function returns a customizable ggplot object displaying exact mean activation counts, bootstrap 95% confidence intervals, and significance annotations. Additional outputs provide exact expectations, bootstrap effect differences, and analysis settings for further examination.

### Parallel moderation screening (`runMgmmAnalysis()`)

A new `ncores` argument distributes candidate moderator analyses across CPU workers, helping reduce computation time. Moderation screening remains available as a model diagnostic for both the Monte Carlo and analytical workflows.

### Corrected moderation intervals

Moderation screening explicitly uses the 2.5th and 97.5th percentiles, consistent with the documented 95% interval. Previous moderation results should be rerun because the earlier implementation used default 90% intervals.

## Interpretation

The analytical workflows distinguish uncertainty from sampling participants from variability introduced by Monte Carlo simulation. Bootstrap intervals describe a fixed absolute perturbation determined from the original fitted network. Ranking stability and significance against zero should not be interpreted as evidence that one intervention target significantly outperforms another.
