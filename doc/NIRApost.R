#' ---
#' title: "NIRApost Tutorial"
#' author: "Wang Fei, Wu Yiming, Wu Yibo, Zhu Tingshao"
#' date: "2025-05-25"
#' output:
#'   rmarkdown::html_vignette:
#'     toc: true
#'     toc_depth: 3
#' vignette: >
#'   %\VignetteIndexEntry{NIRApost Tutorial}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#' 
## ----setup, include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, eval = FALSE)

#' 
#' ## Overview
#' 
#' NIRApost is an R package designed for post-processing of the NodeIdentifyR
#' Algorithm (NIRA) applied to Ising network models in psychological research.
#' It implements a three-step validation procedure:
#' 
#' 1. **Moderation effect testing** (`runMgmmAnalysis`) — a theoretical
#'    prerequisite that must be completed *before* running NIRA
#' 2. **Permutation testing** (`permutationNIRAtest`, `plotMeanNIRA`) — formal
#'    significance testing with multiple-comparison correction
#' 3. **Stability testing** (`stabilityNIRAtest`, `findMaxN`, `plotMaxN`) —
#'    verification that results are reproducible across independent simulation
#'    runs
#' 
#' The package is used in conjunction with **nodeIdentifyR**. For full details,
#' please refer to the companion tutorial article:
#' 
#' > *Simulation Intervention for Cross-Sectional Network Models: An R Tutorial*
#' 
#' ---
#' 
#' ## Package Load
#' 
## ----package_load--------------------------------------------------------------------
library(NIRApost)       # Load NIRApost Package
library(bootnet)        # Load bootnet Package
library(nodeIdentifyR)  # Load nodeIdentifyR Package
library(dplyr)          # Load dplyr Package

#' 
#' ---
#' 
#' ## Step 1: Import Data and Fit Ising Network
#' 
#' We use the built-in `single_gds` dataset (GAD-7 anxiety scale, N = 2,404).
#' An Ising network is estimated via `bootnet::estimateNetwork()` with the "AND"
#' rule. Two components are extracted:
#' 
#' - **Edge weight matrix** — pairwise connection weights between all nodes
#' - **Threshold parameter vector** — intercept terms from each node's logistic
#'   regression in the Ising model
#' 
## ----Ising_Network_Fit---------------------------------------------------------------
# Load the built-in GAD-7 anxiety dataset
data("single_gds")

# Fit the Ising model
gs_fit <- estimateNetwork(single_gds,
                           default = "IsingFit",
                           rule    = "AND")

# Extract edge weight matrix and threshold parameters
edgeWeightMatrix <- gs_fit$graph
thresholdVector  <- gs_fit$intercepts

#' 
#' ---
#' 
#' ## Step 2: Moderation Effect Test (Prerequisite for NIRA)
#' 
#' ### Why This Must Be Done First
#' 
#' NIRA perturbs one node's threshold while keeping all edge weights unchanged.
#' This is valid **only if** no node moderates other pairwise edges. If node A
#' moderates the B–C relationship, then changing A's threshold implicitly alters
#' the B–C edge weight, making the original edge weight matrix theoretically
#' inappropriate for post-intervention simulation.
#' 
#' **This test must be run before NIRA.** An empty `significant_moderators` list
#' confirms the assumption holds.
#' 
#' ### Implementation
#' 
#' `runMgmmAnalysis()` iteratively treats each node as a candidate moderator,
#' fits a Moderated Mixed Graphical Model (MGMM) via `mgm::mgm()`, and uses
#' resampling via `mgm::resample()` to retain only stable effects (95% CI
#' excluding zero).
#' 
#' | Parameter | Description |
#' |---|---|
#' | `data` | Dataset as a numeric matrix |
#' | `plotResults` | `TRUE` exports per-moderator network plots |
#' | `rule` | `"AND"` (stringent, default) or `"OR"` (inclusive) |
#' | `lambdaGam` | EBIC regularization parameter (default: 0.25) |
#' | `nB` | Resampling iterations (use ≥ 100 for publication) |
#' 
## ----Moderation_Effect_Test----------------------------------------------------------
# nB = 10 for demonstration only; use >= 100 for publication
sig_moder <- runMgmmAnalysis(
  data        = as.matrix(single_gds),
  plotResults = FALSE,
  rule        = "AND",
  lambdaGam   = 0.25,
  nB          = 10
)

# Empty list => no moderation effects => NIRA is appropriate
print(sig_moder$significant_moderators)

#' 
#' ---
#' 
#' ## Step 3: NIRA Simulation Intervention
#' 
#' With the moderation assumption confirmed, we proceed with NIRA using
#' functions from **nodeIdentifyR**.
#' 
#' `simulateResponses()` generates N + 1 datasets: for each of the N nodes,
#' one dataset with that node's threshold reduced by 2 SDs (alleviating
#' intervention), plus one unperturbed original dataset.
#' 
#' To perform aggravating interventions, change `"alleviating"` to
#' `"aggravating"` throughout.
#' 
## ----Simulated_alleviate_intervention------------------------------------------------
set.seed(2025)

# Generate N+1 simulated samples
gs_IsingSamples <- simulateResponses(
  edge_weight                = edgeWeightMatrix,
  thresholds                 = thresholdVector,
  perturbation_type          = "alleviating",
  amount_of_SDs_perturbation = 2
)

# Calculate sum scores for each participant in each simulated sample
gs_sumIsingSamples <- calculateSumScores(gs_IsingSamples)

# Convert to long format for plotting
gs_sumIsingSamplesLong <- prepareDFforPlottingAndANOVA(gs_sumIsingSamples)

# Plot mean total scores (ascending order for alleviating interventions)
plotSumScores(
  sum_scores_long   = gs_sumIsingSamplesLong,
  perturbation_type = "alleviating"
)

#' 
#' ---
#' 
#' ## Step 4: Permutation Test
#' 
#' ### Why Visual Inspection Is Insufficient
#' 
#' The `plotSumScores` output has four critical limitations:
#' 
#' 1. **No formal evidence** — visual comparison cannot confirm statistical
#'    significance
#' 2. **CI overlap is unreliable** — checking 95% error bar overlap is an
#'    inaccurate significance criterion
#' 3. **Magnitude sensitivity** — smaller perturbations make visual assessment
#'    highly unreliable
#' 4. **Multiple-comparison inflation** — N simultaneous comparisons with no
#'    correction inflate the Type I error rate
#' 
#' ### Implementation
#' 
#' `permutationNIRAtest()` calls `permutationTest()` for every node in batch,
#' applies multiple-comparison correction, and returns a unified statistics
#' table and plotting-ready data frame. `plotMeanNIRA()` visualises the result
#' as a line plot with 95% CIs and significance stars
#' (*** p < 0.001, ** p < 0.01, * p < 0.05).
#' 
#' | Parameter | Description |
#' |---|---|
#' | `gs_sumIsingSamplesLong` | Long-format output from `prepareDFforPlottingAndANOVA()` |
#' | `method` | Multiple-comparison correction: `"bonferroni"` (default), `"BH"`, `"holm"`, etc. |
#' 
## ----permutationTest-----------------------------------------------------------------
# Batch permutation tests with Bonferroni correction
stat_result <- permutationNIRAtest(gs_sumIsingSamplesLong,
                                    method = "bonferroni")

# View per-node statistics (mean, SD, SE, CI, Cohen's d, p, p.adjust)
head(stat_result$stat)

# Visualise means with 95% CIs and significance stars
plotMeanNIRA(stat_result, "alleviating")

#' 
#' ---
#' 
#' ## Step 5: Stability Test
#' 
#' ### Why Stability Testing Is Necessary
#' 
#' A single simulation run may give different results depending on the random
#' seed. Only results consistent across many independent simulation runs can be
#' trusted.
#' 
#' **Note**: this procedure performs *repeated simulation* — independently
#' re-running the full simulation from scratch using the original network
#' parameters. The observed data are never resampled.
#' 
#' ### Implementation
#' 
#' `stabilityNIRAtest()` repeats the full simulation `nReps` times independently.
#' `findMaxN()` aggregates how often each node ranked in the top-n positions
#' across all repetitions. `plotMaxN()` visualises the stability as a stacked
#' bar plot.
#' 
#' | Parameter | Description |
#' |---|---|
#' | `edge_weights` | Edge weight matrix from the Ising network |
#' | `thresholds` | Threshold parameter vector |
#' | `perturbation_type` | `"aggravating"` or `"alleviating"` |
#' | `amount_of_SDs_perturbation` | Perturbation magnitude in SDs |
#' | `nReps` | Number of independent simulation repetitions |
#' | `parallel` | `TRUE` enables parallel execution (recommended for nReps ≥ 100) |
#' | `ncores` | CPU cores (`NULL` = auto-detect) |
#' 
## ----stability_test------------------------------------------------------------------
# Parallel execution (recommended)
# nReps = 100 for demonstration; use >= 1000 for publication
simResult <- stabilityNIRAtest(
  edge_weights               = edgeWeightMatrix,
  thresholds                 = thresholdVector,
  perturbation_type          = "alleviating",
  amount_of_SDs_perturbation = 2,
  nReps                      = 100,
  parallel                   = TRUE,
  ncores                     = 6
)

# Serial execution (uncomment to use; displays a progress bar)
# simResult <- stabilityNIRAtest(
#   edgeWeightMatrix, thresholdVector,
#   "alleviating", 2,
#   nReps    = 100,
#   parallel = FALSE,
#   ncores   = NULL
# )

# Compute rank frequencies across all simulation repetitions
repeat_node <- findMaxN(simResult, n = 7)
print(repeat_node)

# Stacked bar plot of node ranking stability
plotMaxN(repeat_node, low = 1, high = 7)

#' 
#' ---
#' 
#' ## Summary: Three-Step Validation Procedure
#' 
#' | Step | Purpose | Primary function(s) | Pass criterion |
#' |---|---|---|---|
#' | **Step 2** Moderation test | Verify NIRA assumption | `runMgmmAnalysis()` | `significant_moderators` is empty |
#' | **Step 4** Permutation test | Formal significance + correction | `permutationNIRAtest()`, `plotMeanNIRA()` | Adjusted p < 0.05 |
#' | **Step 5** Stability test | Reproducibility across seeds | `stabilityNIRAtest()`, `findMaxN()`, `plotMaxN()` | Consistent top node |
#' 
#' ---
#' 
#' ## Session Info
#' 
## ----session_info, eval=TRUE, echo=FALSE---------------------------------------------
sessionInfo()

