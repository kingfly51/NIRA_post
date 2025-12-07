#' Stability test for NodeIdentifyR Algorithm (NIRA)
#'
#' @description
#' Repeated the NIRA simulation processes to generate multiple samples for evaluating the stability of 
#' intervention effects in Ising network models. Part of the NodeIdentifyR Algorithm (NIRA) framework for
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
#'   - Typical range: 0-3 SDs
#' @param nboots Integer number of simulation replicates (default: `1000`)
#'   - Recommendation: ≥1000 for publication-quality results
#' @param parallel Logical for parallel processing (default: `FALSE`)
#'   - Set `TRUE` when: nodes >10 *or* nboots >100
#' @param ncores Integer number of CPU cores (default: `NULL`)
#'   - If `NULL` and `parallel=TRUE`, uses `parallel::detectCores() - 1`
#' @param seed Integer seed for random number generation (default: `2025`)
#'   - Default: `2025` for reproducible results
#'   - Always uses the specified seed value
#'
#' @return
#' A list containing:
#' - `BootSamples`: List of length `nboots` containing simulation replicates samples
#' - `mean`: Matrix (`nboots × k`) of category means
#'   - Columns represent network states (original + perturbed conditions)
#' - `sd`: Matrix (`nboots × k`) of category standard deviations
#'
#' @examples
#' ```r
#' # Basic usage with default seed (2025)
#' result1 <- BootNIRAresult(
#'   edge_weights = matrix(runif(25, -0.5, 0.5), 5, 5),
#'   thresholds = rnorm(5),
#'   perturbation_type = "aggravating",
#'   amount_of_SDs_perturbation = 2,
#'   nboots = 50
#'   # Uses default seed = 2025
#' )
#'
#' # Specify different seed
#' result2 <- BootNIRAresult(
#'   edge_weights = my_weights,
#'   thresholds = my_thresholds,
#'   perturbation_type = "alleviating",
#'   amount_of_SDs_perturbation = 1.5,
#'   nboots = 1000,
#'   seed = 1234  # Different seed
#' )
#' ```
#'
#' @author
#' Wang Fei, State Key Laboratory of Cognitive Science and Mental Health, 
#' Institute of Psychology, Chinese Academy of Sciences (2025-12-4)
#'
#' @importFrom parallel detectCores stopCluster makePSOCKcluster
#' @importFrom progress progress_bar
#' @importFrom stats sd
#' @importFrom nodeIdentifyR simulateResponses calculateSumScores prepareDFforPlottingAndANOVA
#' @export

BootNIRAresult <- function(edge_weights, thresholds, perturbation_type,
                  amount_of_SDs_perturbation, nboots =1000,
                  parallel = FALSE, ncores = NULL, 
                  seed = 2025) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  Bootresult <- list(
    BootSamples = vector(mode = "list", length = nboots),
    mean = list(),
    sd = list()
  )
  
  if (parallel) {
    if (is.null(ncores)) {
      ncores <- parallel::detectCores() - 1
    }
    
    cl <- parallel::makeCluster(ncores)
    
    parallel::clusterSetRNGStream(cl, seed)
    
    parallel::clusterEvalQ(cl, {
      library(nodeIdentifyR)
      library(dplyr)
    })
    
    parallel::clusterExport(cl, 
                           varlist = c("edge_weights", "thresholds", 
                                      "perturbation_type", "amount_of_SDs_perturbation"),
                           envir = environment())
    
    all_results <- parallel::parLapply(cl, 1:nboots, function(i) {
      gs_IsingSamples <- nodeIdentifyR::simulateResponses(edge_weights, thresholds,
                                                          perturbation_type, amount_of_SDs_perturbation)
      gs_sumIsingSamples <- nodeIdentifyR::calculateSumScores(gs_IsingSamples)
      gs_sumIsingSamplesLong <- nodeIdentifyR::prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
      
      categories <- sort(unique(gs_sumIsingSamplesLong[[2]]))
      category_means <- sapply(categories, function(cat) {
        mean(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
      })
      category_sd <- sapply(categories, function(cat) {
        sd(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
      })
      
      list(
        sample = gs_sumIsingSamplesLong,
        means = category_means,
        sds = category_sd
      )
    })
    
    parallel::stopCluster(cl)
    
    Bootresult$BootSamples <- lapply(all_results, function(x) x$sample)
    Bootresult$mean <- do.call(rbind, lapply(all_results, function(x) x$means))
    Bootresult$sd <- do.call(rbind, lapply(all_results, function(x) x$sds))
    
    if (length(Bootresult$BootSamples) > 0) {
      categories <- sort(unique(Bootresult$BootSamples[[1]][[2]]))
      colnames(Bootresult$mean) <- categories
      colnames(Bootresult$sd) <- categories
    }
    
  } else {
    if (interactive()) {
      pb <- progress::progress_bar$new(
        format = "[:bar] :percent completed: :eta remaining time",
        total = nboots,
        clear = FALSE,
        width = 60
      )
    }
    
    all_means <- list()
    all_sds <- list()
    
    for (i in 1:nboots) {
      gs_IsingSamples <- nodeIdentifyR::simulateResponses(edge_weights, thresholds,
                                                          perturbation_type, amount_of_SDs_perturbation)
      gs_sumIsingSamples <- nodeIdentifyR::calculateSumScores(gs_IsingSamples)
      gs_sumIsingSamplesLong <- nodeIdentifyR::prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
      
      Bootresult$BootSamples[[i]] <- gs_sumIsingSamplesLong
      
      categories <- sort(unique(gs_sumIsingSamplesLong[[2]]))
      category_means <- sapply(categories, function(cat) {
        mean(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
      })
      category_sd <- sapply(categories, function(cat) {
        sd(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
      })
      
      all_means[[i]] <- category_means
      all_sds[[i]] <- category_sd
      
      if (interactive()) {
        pb$tick()
      }
    }
    
    Bootresult$mean <- do.call(rbind, all_means)
    Bootresult$sd <- do.call(rbind, all_sds)
    
    if (length(Bootresult$BootSamples) > 0) {
      categories <- sort(unique(Bootresult$BootSamples[[1]][[2]]))
      colnames(Bootresult$mean) <- categories
      colnames(Bootresult$sd) <- categories
    }
  }
  
  return(Bootresult)
}