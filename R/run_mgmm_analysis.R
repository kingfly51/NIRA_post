#' Moderated Graphical Model Analysis
#'
#' @description
#' Performs comprehensive moderated network analysis using Mixed Graphical Models (MGMM),
#' identifying significant moderator variables and their interaction effects through
#' bootstrap validation.
#'
#' @param data Numeric matrix or data frame (n × p) where:
#'   - Rows: Observations
#'   - Columns: Variables (nodes)
#' @param plot_results Logical indicating whether to save diagnostic plots for each moderator (default: `FALSE`)
#' @param rule Regularization rule for edge selection:
#'   - `"AND"` (default): More conservative, both directions must be significant
#'   - `"OR"`: More inclusive, either direction significant
#' @param lambdaGam EBIC tuning parameter (default: 0.25). Values:
#'   - 0: BIC-like penalty
#'   - 1: EBIC with full penalty
#' @param nB Number of bootstrap samples (default: 100). For stable results:
#'   - ≥100 for exploratory analysis
#'   - ≥1000 for publication
#'
#' @return
#' A list containing:
#'
#' **1. Full Results (`all_results`)**
#' List of length p (number of nodes), each element contains:
#' - `model`: Fitted mgm object
#' - `interactions`: Detected interaction effects
#' - `bootstrap`: Bootstrap results object
#' - `table`: Summary table of edge weights and CIs
#'
#' **2. Significant Moderators (`significant_moderators`)**
#' Filtered list where each element contains:
#' - `node`: Moderator index
#' - `significant_edges`: Indices of significantly moderated edges
#' - `interaction_details`: Statistical details for significant interactions
#'
#' @examples
#' ```r
#' # Simulate binary data
#' set.seed(123)
#' dat <- matrix(rbinom(100*10, 1, 0.5), ncol = 10)
#'
#' # Run analysis (minimal example)
#' results <- run_mgmm_analysis(dat, nB = 50)
#'
#' # Extract significant moderators
#' sig_mods <- results$significant_moderators
#' if(length(sig_mods) > 0) print(sig_mods[[1]]$interaction_details)
#' ```
#'
#' @author
#' Wang Fei, State Key Laboratory of Cognitive Neuroscience and Learning,
#' Beijing Normal University (2025-03-22)
#'
#' @section Algorithm Workflow:
#' 1. **Model Fitting**:
#'    - Estimates separate MGMs with each node as moderator
#'    - Uses EBIC for sparse selection
#' 2. **Bootstrap Validation**:
#'    - Case resampling (nB iterations)
#'    - Computes 95% CIs for edge weights
#' 3. **Significance Testing**:
#'    - Identifies edges where CI doesn't span zero
#'
#' @seealso
#' * [`mgm::mgm()`] for core modeling function
#' * [`bootnet::bootnet()`] for alternative network bootstrap methods
#'
#' @importFrom mgm mgm resample plotRes
#' @importFrom grDevices png dev.off
#' @importFrom utils capture.output
#' @export

run_mgmm_analysis <- function(data, plot_results = FALSE, rule = "AND",
                             lambdaGam = 0.25, nB = 100) {

  # Input validation
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("data must be a matrix or data frame")
  }

  if (!rule %in% c("AND", "OR")) {
    stop('rule must be either "AND" or "OR"')
  }

  if (lambdaGam < 0 || lambdaGam > 1) {
    stop("lambdaGam must be between 0 and 1")
  }

  if (nB < 10) {
    warning("nB < 10 may yield unreliable results")
  }

  node_num = ncol(data)
  # Initialization result list
  results_list <- vector("list", node_num)
  names(results_list) <- paste0("Moderator_", 1:node_num)

  # Run analysis on each node as a moderating variable
  for (i in 1:node_num) {
    cat("\n------ Running analysis with Moderator", i, "------\n")

    # Fit MGMM model
    mgm_mod <- mgm::mgm(
      data = as.matrix(data),
      type = rep("c", node_num),
      level = rep(2, node_num),
      lambdaSel = "EBIC",
      lambdaGam = lambdaGam,
      ruleReg = rule,
      moderators = i,
      scale = FALSE
    )

    # The bootstrap resampling
    boot_res <- mgm::resample(
      object = mgm_mod,
      data = as.matrix(data),
      nB = nB
    )

    # Obtain result table
    tab <- mgm::plotRes(object = boot_res, table = TRUE)
    tab <- tab[,-c(7,8)]

    # Store results
    results_list[[i]] <- list(
      model = mgm_mod,
      interactions = mgm_mod$interactions$indicator[[2]],
      bootstrap = boot_res,
      table = tab
    )

    # If drawing is required
    if (plot_results) {
      png(paste0("Moderator_", i, "_Plot.png"), width = 800, height = 600)
      plotRes(object = boot_res)
      dev.off()
    }
  }

  # Identify significant moderating variables
  significant_moderators <- list()

  for (i in 1:node_num) {
    boot_summary <- results_list[[i]]$table
    sig_edges <- which(boot_summary[, "Mod_qtl_low"] > 0 | boot_summary[, "Mod_qtl_high"] < 0)

    if (length(sig_edges) > 0) {
      significant_moderators[[paste0("Moderator_", i)]] <- list(
        node = i,
        significant_edges = sig_edges,
        interaction_details = results_list[[i]]$table[sig_edges, ]
      )
    }
  }

  return(list(
    all_results = results_list,
    significant_moderators = significant_moderators
  ))
}
