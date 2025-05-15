run_mgmm_analysis <- function(data, node_num, plot_results = FALSE, lambdaGam = 0.25, nB = 100) {
  #'Perform MGMM analysis and detect significant moderating variables
  #'
  #' @param data data matrix
  #' @param node_num Number of network nodes
  #' @param plot_results Whether to draw the result graph, default is False
  #' @param lambdaGam EBIC parameter, default is 0.25
  #' @param nB The bootstrap resampling is conducted with 100 replicates by default.
  #'
  #' @return List containing significant moderating variables
  
  # Load necessary packages
  if (!require(mgm)) {
    stop("Please install the mgm package first: install.packages('mgm')")
  }
  
  # Initialization result list
  results_list <- vector("list", node_num)
  names(results_list) <- paste0("Moderator_", 1:node_num)
  
  # Run analysis on each node as a moderating variable
  for (i in 1:node_num) {
    cat("\n------ Running analysis with Moderator", i, "------\n")
    
    # Fit MGMM model
    mgm_mod <- mgm(
      data = as.matrix(data),
      type = rep("c", node_num),
      level = rep(2, node_num),
      lambdaSel = "EBIC", 
      lambdaGam = lambdaGam,
      ruleReg = "AND",
      moderators = i,  
      scale = FALSE  
    )
    
    # The bootstrap resampling
    boot_res <- resample(
      object = mgm_mod,
      data = as.matrix(data),
      nB = nB
    )
    
    # Obtain result table
    tab <- plotRes(object = boot_res, table = TRUE)
    
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