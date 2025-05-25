#' Bootstrap Analysis for NodeIdentifyR Algorithm (NIRA)
#'
#' @description
#' Performs case-resampling bootstrap analysis to evaluate the stability of intervention effects
#' in Ising network models. Part of the NodeIdentifyR Algorithm (NIRA) framework for
#' psychological network analysis.
#'
#' @param edge_weights Numeric matrix (`p × p`) representing edge weights in the Ising network.
#'   - Should be symmetric with zero diagonal
#'   - Positive values indicate excitatory connections
#'   - Negative values indicate inhibitory connections
#' @param thresholds Numeric vector (length `p`) of threshold parameters (intercepts) for nodes
#' @param perturbation_type Intervention direction:
#'   - `"aggravating"` (+): Increases symptom severity
#'   - `"alleviating"` (-): Decreases symptom severity
#' @param amount_of_SDs_perturbation Numeric magnitude of perturbation (in SDs of threshold distribution)
#'   - Typical range: 1-3 SDs
#' @param nboots Integer number of bootstrap replicates (default: `1000`)
#'   - Recommendation: ≥1000 for publication-quality results
#' @param parallel Logical for parallel processing (default: `FALSE`)
#'   - Set `TRUE` when: nodes >10 *or* nboots >100
#' @param ncores Integer number of CPU cores (default: `NULL`)
#'   - If `NULL` and `parallel=TRUE`, uses `parallel::detectCores() - 1`
#'
#' @return
#' A list containing:
#' - `BootSamples`: List of length `nboots` containing bootstrap samples
#' - `mean`: Matrix (`nboots × k`) of category means
#'   - Columns represent network states (original + perturbed conditions)
#' - `sd`: Matrix (`nboots × k`) of category standard deviations
#'
#' @examples
#' ```r
#' # Basic usage
#' set.seed(123)
#' result <- BootNIRAresult(
#'   edge_weights = matrix(runif(25, -0.5, 0.5), 5, 5),
#'   thresholds = rnorm(5),
#'   perturbation_type = "aggravating",
#'   amount_of_SDs_perturbation = 2,
#'   nboots = 50  # Use >=1000 in practice
#' )
#'
#' # Parallel computation
#' par_result <- BootNIRAresult(
#'   edge_weights = my_weights,
#'   thresholds = my_thresholds,
#'   perturbation_type = "alleviating",
#'   amount_of_SDs_perturbation = 1.5,
#'   nboots = 1000,
#'   parallel = TRUE,
#'   ncores = 4
#' )
#' ```
#'
#' @author
#' Wang Fei, State Key Laboratory of Cognitive Neuroscience and Learning,
#' Beijing Normal University (2025-03-22)
#'
#' @importFrom parallel detectCores stopCluster makePSOCKcluster makeForkCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar% %do%
#' @importFrom progress progress_bar
#' @importFrom stats sd
#' @importFrom magrittr %>%
#' @importFrom nodeIdentifyR simulateResponses calculateSumScores prepareDFforPlottingAndANOVA
#' @export

BootNIRAresult <- function(edge_weights, thresholds, perturbation_type,
                           amount_of_SDs_perturbation, nboots,
                           parallel = FALSE, ncores = NULL) {

  Bootresult <- list(
    BootSamples = vector(mode = "list", length = nboots),
    mean = list(),
    sd = list()
  )


  if (parallel) {
    if (is.null(ncores)) {
      ncores <- parallel::detectCores() - 1
    }

    cl <- if (.Platform$OS.type == "unix") {
      makeForkCluster(ncores)
    } else {
      makePSOCKcluster(ncores)
    }

    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  }


  if (!parallel && interactive()) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent completed: :eta remaining time",
      total = nboots,
      clear = FALSE,
      width = 60
    )
  }


  combine_results <- function(list1, list2) {
    list(
      sample = c(list1$sample,list2$sample),
      means = rbind(list1$means, list2$means),
      sds = rbind(list1$sds, list2$sds)
    )
  }

  if (parallel){
    result <- foreach(i = 1:nboots, .combine = "combine_results") %dopar% {
      gs_IsingSamples <- nodeIdentifyR::simulateResponses(edge_weights, thresholds,
                                                          perturbation_type, amount_of_SDs_perturbation)
      gs_sumIsingSamples <- nodeIdentifyR::calculateSumScores(gs_IsingSamples)
      library(dplyr)
      gs_sumIsingSamplesLong <- nodeIdentifyR::prepareDFforPlottingAndANOVA(gs_sumIsingSamples)

      categories <- sort(unique(gs_sumIsingSamplesLong[[2]]))
      category_means <- sapply(categories, function(cat) {
        mean(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
      })
      category_sd <- sapply(categories, function(cat) {
        sd(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
      })


      if (!parallel && interactive()) {
        pb$tick()
      }

      list(
        sample = list(gs_sumIsingSamplesLong),
        means = category_means,
        sds = category_sd
      )
    }
    } else {
      result <- foreach(i = 1:nboots, .combine = "combine_results") %do% {
        gs_IsingSamples <- nodeIdentifyR::simulateResponses(edge_weights, thresholds,
                                                            perturbation_type, amount_of_SDs_perturbation)
        gs_sumIsingSamples <- nodeIdentifyR::calculateSumScores(gs_IsingSamples)
        library(dplyr)
        gs_sumIsingSamplesLong <- nodeIdentifyR::prepareDFforPlottingAndANOVA(gs_sumIsingSamples)

        categories <- sort(unique(gs_sumIsingSamplesLong[[2]]))
        category_means <- sapply(categories, function(cat) {
          mean(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
        })
        category_sd <- sapply(categories, function(cat) {
          sd(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
        })
        pb$tick()
        list(
          sample = list(gs_sumIsingSamplesLong),
          means = category_means,
          sds = category_sd
        )
      }
    }

  Bootresult$BootSamples <- result$sample
  Bootresult$mean <-result$means
  Bootresult$sd <-result$sd
  categories <- sort(unique(result[["sample"]][[1]]$sample))
  colnames(Bootresult$mean) <- categories
  colnames(Bootresult$sd) <- categories
  return(Bootresult)
}





