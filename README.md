# NIRA_post

NIRApost is an R package for post-processing of the NodeIdentifyR Algorithm
(NIRA) applied to Ising network models in psychological research. It provides
**two complementary workflows** for single-network and dual-network (bridge)
intervention analysis:

- **Workflow 1 — Monte Carlo simulation** (the original NIRA pipeline):
  network estimation → simulation of original and post-intervention samples →
  permutation testing → repeated-simulation stability assessment.
- **Workflow 2 — Analytical solution** (new in v1.2.0): the expected number of
  activated nodes is computed *exactly* by enumerating all 2^p Ising
  microstates, eliminating Monte Carlo error. Uncertainty is quantified by
  resampling the **original dataset** (bootstrap + random-target test).

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

### Workflow 2 — Analytical solution (`analyticalNIRAtest`)

A single function, `analyticalNIRAtest()`, computes the exact mean activation
under any one-node intervention by enumerating all 2^p microstates
(Cui et al., 2026). Because the result is exact given the network parameters,
the only remaining uncertainty is **sampling error**, which is quantified by
resampling the **original dataset** (not simulated samples):

- `$stat` — bootstrap inference per node: exact mean activation, intervention
  effect (Δ), Cohen's *d*, SE, 95% CI, *p*-value for H0: Δ = 0, adjusted
  *p*-value, and `percenttop_1` (ranking stability).
- `$random` — random-target test: is node *i* a better target than a randomly
  chosen node? (uncorrected).
- `$exact` — exact mean activation for the original network and each
  intervention.
- `$plot` — ggplot of exact mean activation with 95% CI and significance
  stars.

### Dual-network (bridge) NIRA — `analyticalBridgeNIRAtest`

For a network partitioned into two communities (e.g., health behaviors vs.
mental-health symptoms), `analyticalBridgeNIRAtest()` intervenes on each node
of a **source** community and computes the exact mean activation of the
**outcome** community. See **Section 6** for a worked interpretation.

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
including the analytical workflow functions `analyticalNIRAtest()` and
`analyticalBridgeNIRAtest()`.

---

## 6. Interpreting bridge-NIRA results (worked example)

The following example intervenes on each node of a *source* community
(JKXW, 8 health-behavior nodes) and reports the exact mean activation of the
*outcome* community (MH, 3 mental-health nodes: depression, anxiety, stress),
using an aggravating intervention (threshold raised, node more likely to be
activated).

```r
res2 <- analyticalBridgeNIRAtest(
  data               = GS_non_null,
  groups             = groups,                # rep(c("JKXW", "MH"), c(8, 3))
  intervention_group = "JKXW",
  outcome_group      = "MH",
  perturbation_type  = "aggravating",
  edge_weights       = edgeWeightMatrix,      # fit$graph from bootnet
  thresholds         = thresholdVector,       # fit$intercepts from bootnet
  amount_of_SDs_perturbation = 2,
  nResample          = 1000,
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
Δ = 0.685 − 0.646 = 0.039 (Cohen's *d* = 0.041).

`$stat` then summarises the same quantities **per source node** together with
their resampling-based uncertainty:

- `effect` / `cohen_d` — exact intervention effect and its standardised size.
- `se_effect`, `ci_effect_lower/upper` — SE and 95% bootstrap CI of Δ.
- `p`, `p.adjust` — bootstrap *p*-value for H0: Δ = 0, and the
  Bonferroni-corrected value.
- `percenttop_1` — the proportion of resamples in which the node had the
  largest absolute effect on the outcome community. This is the analytical
  analogue of `findMaxN()` and the *discriminating* metric here: because the
  analytical solution makes every Δ deterministically non-zero, the H0: Δ = 0
  test is trivially non-significant for all nodes, whereas `percenttop_1`
  ranks which source is most consistently the strongest bridge target.
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
`$random` approaches the same question from the target-comparison angle: for
each source node it gives the (uncorrected) probability that a *randomly
chosen source* node would produce an effect at least as large. In this
example all random-target *p*-values are large (well above 0.05), indicating
that none of the health-behavior nodes is demonstrably superior to a random
source target — consistent with the weak bridge edges recovered in the
estimated network.

### 6.2 How to report it in a paper

A concise results paragraph could be written as follows:

> Using the analytical bridge-NIRA workflow, we intervened on each health-
> behavior node and computed the exact expected number of activated
> mental-health symptoms (depression, anxiety, stress). In the original
> network the outcome community was active at a mean of 0.65 nodes (SD-based
> scale 0–3). Aggravating intervention on 7 of the 8 behavior nodes produced
> negligible changes in mental-health activation (Δ < 0.004, all adjusted
> *p* = 1). Only intervening on too-little-sleep increased the expected
> number of activated symptoms, by Δ = 0.039 (Cohen's *d* = 0.041), but this
> effect was not statistically reliable (95% CI [−0.001, 0.158]; bootstrap
> *p* = 0.489, Bonferroni-adjusted *p* = 1). Too-little-sleep nevertheless
> ranked first in 51.4% of the 1000 resamples (`percenttop_1`), suggesting it
> is the most consistently promising bridge target, although its superiority
> over a randomly chosen behavior node did not reach significance in the
> random-target test.

**Reporting checklist.** (1) state the outcome community and the scale of
`mean_activation`; (2) report the original (pre-intervention) value; (3)
report Δ, Cohen's *d*, 95% CI, and the adjusted *p*-value for each target, or
at least for the top-ranked nodes; (4) interpret `percenttop_1` as ranking
stability rather than significance; (5) note the model assumption that the
moderation-effect test (`runMgmmAnalysis()`) found no moderators, so
intervening on one node leaves all edge weights unchanged.
