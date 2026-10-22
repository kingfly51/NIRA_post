# NIRA_post

NIRApost is an R package for post-processing of the NodeIdentifyR Algorithm
(NIRA) applied to Ising network models in psychological research. It provides
**two complementary workflows** for single-network and dual-network (bridge)
intervention analysis:

- **Workflow 1 — Monte Carlo simulation** (the original NIRA pipeline):
  network estimation → simulation of original and post-intervention samples →
  permutation testing → repeated-simulation stability assessment.
- **Workflow 2 — Analytical solution** (broadened in v1.3.0): the expected
  number of activated nodes is computed *exactly* by enumerating all 2^p Ising
  microstates, eliminating Monte Carlo error. The exact functions report
  deterministic quantities and draw exact plots **without confidence
  intervals**; sampling uncertainty and ranking stability are handled by
  separate bootstrap functions that plot a stability graph analogous to the
  Monte Carlo stability test.

The `single_gds_NIRA.R` file is the complete set of code appearing in the
Simulation Intervention for Cross-Sectional Network Models article
(Wang et al., 2026), and `Supplementary Material Code.R` is the code used in
the article's supplementary materials. `Bridge_network_NIRA.R` and
`single_gds_analyticalNIRA.R` demonstrate the two workflows.

---

## 1. Why use NIRA?

1. Traditional centrality indices (e.g., strength centrality) only reflect
   structural information and may not fully capture directional influence
   patterns between nodes (e.g., in-influence vs. out-influence). NIRA
   addresses this by simulating interventions to estimate projected impacts
   on the network.
2. By simulating node removal/adjustment, NIRA quantifies each node's
   projected global effect on the network, complementing centrality-based
   approaches (e.g., high-centrality nodes may not always be projected
   optimal targets).
3. Cross-sectional models lack explicit causal pathways. NIRA infers
   potential relationships through simulated interventions, which may aid in
   identifying core or bridging symptoms.

---

## 2. The two workflows

### Workflow 1 — Monte Carlo simulation

Functions: `permutationNIRAtest()`, `stabilityNIRAtest()`, `findMaxN()`,
`plotMaxN()`, `plotMeanNIRA()`, preceded by the moderation-effect test
`runMgmmAnalysis()`. This workflow quantifies Monte Carlo error and is
recommended for larger networks where 2^p enumeration is infeasible.

### Workflow 2 — Analytical solution

The analytical workflow is split into two functions so that exact results and
sampling uncertainty are never mixed in one plot.

**Exact solution — `analyticalNIRAtest()`**

Computes the exact mean activation under any one-node intervention by
enumerating all 2^p microstates (Cui et al., 2026):

- `$stat` — exact mean activation, intervention effect (Δ) and Glass's Δ
  (`glass_d`) per node.
- `$exact` — exact mean activation for the original network and each
  intervention.
- `$plot` — ggplot of exact mean activation. It deliberately contains **no
  confidence interval and no significance stars**.

**Bootstrap stability — `analyticalStabilityNIRAtest()`**

Quantifies the only remaining uncertainty (**sampling error**) by resampling
the **original dataset** (not simulated samples), re-estimating the network
and recomputing the exact solution on every resample:

- `$stat` — exact effects plus bootstrap SE, 95% CI of the mean activation
  and of Δ, *p*-value for H0: Δ = 0, adjusted *p*-value, and
  `percenttop_1 … percenttop_p` (ranking stability).
- `$resample_effects` — the resampled effect matrix.
- `$maxn` — rank-frequency table ready for `plotMaxN()`.
- `$plot` — stacked-bar **ranking-stability plot** across bootstrap
  resamples, analogous to the Monte Carlo stability test
  (`stabilityNIRAtest()` + `plotMaxN()`).

### Dual-network (bridge) NIRA

For a network partitioned into two communities (e.g., health behaviors vs.
mental-health symptoms), `analyticalBridgeNIRAtest()` intervenes on each node
of a **source** community and computes the exact mean activation of the
**outcome** community; `analyticalBridgeStabilityNIRAtest()` performs the
corresponding bootstrap stability analysis. See **Section 6** for a worked
interpretation.

As in the other workflows, the moderation-effect test
(`runMgmmAnalysis()`) is a **prerequisite** and must be run on the **whole**
network (both communities together) before the bridge intervention: the
threshold-shift intervention model applies to the full joint distribution, and
an intervention on a source node propagates through the entire network, so a
moderator anywhere would invalidate the fixed-edge-weight assumption. The
`Bridge_network_NIRA.R` script includes this step (Step 3).

### Consistent resampling estimator

When network parameters are supplied (typically from
`bootnet::estimateNetwork`), the bootstrap re-estimates the network with the
**same** estimator by default (`resample_default = NULL`), so confidence
intervals are correctly centred on the point estimate. Use
`resample_default = "fast"` only when the point estimate also used the fast
estimator.

---

## 3. How to use NIRApost

For the complete tutorial, please check vignettes/NIRApost-tutorial.Rmd

---

## 4. How to install NIRApost

Users can install the NIRApost package using the following command. Of course,
the prerequisite is that the user needs to install the `devtools` package
first. In addition, users can also download and install locally.

```r
install.packages("devtools")

devtools::install_github("kingfly51/NIRA_post")
```

---

## 5. Introduction to each file in NIRApost

The `data/` directory contains five built-in datasets distributed with the
NIRApost package, stored in `.rda` format for immediate accessibility upon
package loading. Users can directly load these datasets into their R
environment—for example, the `single_gds` dataset can be accessed via:

```r
data("single_gds")
```

The `data_raw/` directory stores the original `.xlsx` files for all five
datasets.

The `doc/` directory contains four key files:
1. `NIRApost-tutorial.html`: A self-contained webpage providing a
   comprehensive tutorial on using the NIRApost package.
2. `NIRApost-tutorial.pdf`: The PDF document used to generate the tutorial,
   containing all code, text, and formatting instructions.
3. `NIRApost.R`: An R script version of the tutorial code (without Markdown
   text), suitable for direct execution or adaptation.
4. `index.html`: Main interface, when the browser is opened, it can directly
   jump to three other files.

The `man/` directory contains the help documentation files for NIRApost's
core functions. Taking the `permutationNIRAtest()` function as an example,
users can access its documentation in R using the following command:

```r
help(permutationNIRAtest, package = "NIRApost")
```

The `R/` directory contains all core functions of the NIRApost package,
including the exact analytical functions `analyticalNIRAtest()` and
`analyticalBridgeNIRAtest()` and their bootstrap stability counterparts
`analyticalStabilityNIRAtest()` and `analyticalBridgeStabilityNIRAtest()`.

---

## 6. Interpreting bridge-NIRA results (worked example)

The following example intervenes on each node of a *source* community
(JKXW, 8 health-behavior nodes) and reports the exact mean activation of the
*outcome* community (MH, 3 mental-health nodes: depression, anxiety, stress),
using an aggravating intervention (threshold raised, node more likely to be
activated).

```r
# 1) Exact solution: deterministic, no confidence intervals
res2 <- analyticalBridgeNIRAtest(
  data               = GS_non_null,
  groups             = groups,                # rep(c("JKXW", "MH"), c(8, 3))
  intervention_group = "JKXW",
  outcome_group      = "MH",
  perturbation_type  = "aggravating",
  edge_weights       = edgeWeightMatrix,      # fit$graph from bootnet
  thresholds         = thresholdVector,       # fit$intercepts from bootnet
  amount_of_SDs_perturbation = 2
)

# 2) Bootstrap stability: sampling uncertainty + ranking-stability plot
res2_stab <- analyticalBridgeStabilityNIRAtest(
  data               = GS_non_null,
  groups             = groups,
  intervention_group = "JKXW",
  outcome_group      = "MH",
  perturbation_type  = "aggravating",
  edge_weights       = edgeWeightMatrix,
  thresholds         = thresholdVector,
  amount_of_SDs_perturbation = 2,
  nReps              = 1000,
  method             = "bonferroni",
  seed               = 2025,
  parallel           = TRUE,
  ncores             = 6
)
```

### 6.1 How to read `$exact` and `$stat`

`$exact` gives the exact mean number of activated outcome nodes under each
condition:

| condition | mean_activation |
|---|---|
| original | 0.646 |
| 喝酒 / 吸烟 / 久坐 / 熬夜 / 睡太多 / 独居 / 不吃早餐 | 0.642–0.648 |
| 睡太少 | 0.685 |

**Interpretation.** In this sample the mental-health outcome community
(depression, anxiety, stress) is on average active at only ~0.65 nodes out of
3, i.e. most participants report no more than one symptom. Intervening on
most health-behavior nodes leaves the mental-health activation essentially
unchanged (Δ ≈ 0), because the estimated bridge edges from the JKXW community
to the MH community are weak and mostly pruned by the regularized estimator.
The exception is **睡太少 (too little sleep)**: raising this node's activation
increases the expected number of activated mental-health nodes by
Δ = 0.685 − 0.646 = 0.039 (Glass's Δ = 0.041).

The exact `$stat` gives `effect` / `glass_d` **per source node** — the
intervention effect and its standardised size (Glass's Δ: Δ divided by the
exact SD of the *original* activation distribution, **not** a pooled-SD
Cohen's *d*). Both carry the sign of the intervention: positive for
aggravating, negative for alleviating.

`analyticalBridgeStabilityNIRAtest()` adds the resampling-based uncertainty:

- `se_effect`, `ci_effect_lower/upper` — SE and 95% bootstrap CI of Δ.
- `ci_lower`, `ci_upper` — 95% bootstrap CI of the mean activation.
- `p`, `p.adjust` — bootstrap *p*-value for H0: Δ = 0, and the
  Bonferroni-corrected value.
- `percenttop_1 … percenttop_p` — the proportion of resamples in which the
  node occupied each rank by absolute effect on the outcome community. This
  is the analytical analogue of `findMaxN()`, and `$plot` draws the
  corresponding stacked-bar stability graph.
  Note that, like `findMaxN()`, `percenttop_1` scores by *absolute* effect
  (mirroring the classical stability metric). Under an alleviating
  intervention a node whose effect runs counter to the intended direction
  (positive Δ, e.g. a weak negative bridge edge) is still counted by its
  magnitude; in practice such counter-directional effects are typically
  negligible, but check the sign of `effect` when a node with a
  counter-directional Δ has a high `percenttop_1`.

**Reading the example output.** 睡太少 has by far the largest `percenttop_1`
(0.514, i.e. it ranked first in 51.4% of the 1000 resamples), followed by
不吃早餐 (0.174) and 独居 (0.129). All bootstrap *p*-values are non-significant
(`p.adjust = 1` for every node), so in this sample no single health-behavior
node produces a statistically reliable increase in mental-health activation.
`percenttop_1` is the metric that discriminates here: it ranks the source
nodes by how consistently each is the strongest bridge target, which is the
question of interest for target selection given that every Δ is
deterministically non-zero.

### 6.2 How to report it in a paper

A concise results paragraph could be written as follows:

> Using the analytical bridge-NIRA workflow, we intervened on each health-
> behavior node and computed the exact expected number of activated
> mental-health symptoms (depression, anxiety, stress). In the original
> network the outcome community was active at a mean of 0.65 nodes (SD-based
> scale 0–3). Aggravating intervention on 7 of the 8 behavior nodes produced
> negligible changes in mental-health activation (Δ < 0.004, all adjusted
> *p* = 1). Only intervening on too-little-sleep increased the expected
> number of activated symptoms, by Δ = 0.039 (Glass's Δ = 0.041), but this
> effect was not statistically reliable (95% CI [−0.001, 0.158]; bootstrap
> *p* = 0.489, Bonferroni-adjusted *p* = 1). Too-little-sleep nevertheless
> ranked first in 51.4% of the 1000 resamples (`percenttop_1`), suggesting it
> is the most consistently promising bridge target.

**Reporting checklist.** (1) state the outcome community and the scale of
`mean_activation`; (2) report the original (pre-intervention) value; (3)
report Δ, Glass's Δ, 95% CI, and the adjusted *p*-value for each target, or
at least for the top-ranked nodes; (4) interpret `percenttop_1` as ranking
stability rather than significance; (5) note the model assumption that the
moderation-effect test (`runMgmmAnalysis()`) found no moderators, so
intervening on one node leaves all edge weights unchanged.
