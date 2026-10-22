# Internal helpers shared by the analytical NIRA functions.
# Nothing in this file is exported.

# Validate and normalise binary input data.  Returns the complete-case matrix,
# the number of nodes and the node labels.
analytical_prepare_data <- function(data) {
  if (!is.matrix(data) && !is.data.frame(data))
    stop("'data' must be a binary matrix or data frame.")
  data <- as.matrix(data)
  if (anyNA(data)) {
    keep <- stats::complete.cases(data)
    warning(sum(!keep), " row(s) with missing values removed.")
    data <- data[keep, , drop = FALSE]
  }
  if (!all(data %in% c(0, 1)))
    stop("'data' must contain only binary entries (0/1).")
  if (nrow(data) < 3 || ncol(data) < 2)
    stop("'data' must contain at least 3 rows and 2 columns.")

  p <- ncol(data)
  if (is.null(colnames(data))) colnames(data) <- paste0("V", seq_len(p))
  node_names <- colnames(data)
  if (anyNA(node_names) || any(!nzchar(node_names)) ||
      anyDuplicated(node_names) || "original" %in% node_names)
    stop("Node names must be unique, nonempty, and different from 'original'.")

  if (2^p > 2^20)
    warning("2^p = ", 2^p, " microstates; exact enumeration may be slow for p > 20.")

  list(data = data, p = p, node_names = node_names)
}

# Saturated nodewise logistic regression (pseudo-likelihood) estimator,
# roughly 25x faster than bootnet::estimateNetwork.  Used for the point
# estimate when parameters are not supplied and inside the resampling loop.
estimate_ising_fast <- function(X) {
  pn   <- ncol(X)
  thr  <- numeric(pn)
  W    <- matrix(0, pn, pn)
  for (j in seq_len(pn)) {
    others <- seq_len(pn)[-j]
    g <- stats::glm.fit(x = cbind(1, X[, others, drop = FALSE]),
                        y = X[, j],
                        family = stats::binomial())
    thr[j]       <- g$coefficients[1]
    W[j, others] <- g$coefficients[-1]
  }
  W <- (W + t(W)) / 2
  list(graph = W, thresholds = thr)
}

# Fit the Ising network when parameters are not supplied, or validate supplied
# parameters.  Returns finite, symmetric edge weights and thresholds.
analytical_fit_network <- function(data, thresholds, edge_weights, p) {
  if (is.null(thresholds)) {
    est          <- estimate_ising_fast(data)
    thresholds   <- est$thresholds
    edge_weights <- est$graph
  } else {
    thresholds   <- as.numeric(thresholds)
    edge_weights <- as.matrix(edge_weights)
    if (length(thresholds) != p)
      stop("'thresholds' must have length p = ", p, ".")
    if (!all(dim(edge_weights) == c(p, p)))
      stop("'edge_weights' must be a ", p, " x ", p, " matrix.")
    if (any(diag(edge_weights) != 0))
      warning("Diagonal of 'edge_weights' set to 0.")
    diag(edge_weights) <- 0
  }

  if (length(thresholds) != p || !all(is.finite(thresholds)) ||
      !is.numeric(edge_weights) || !all(is.finite(edge_weights)))
    stop("Network estimation must produce finite thresholds and edge weights; check constant or collinear columns.")
  if (!isTRUE(all.equal(unname(edge_weights), unname(t(edge_weights)),
                        tolerance = 1e-10)))
    stop("'edge_weights' must be symmetric.")

  list(thresholds = thresholds, edge_weights = edge_weights)
}

# Resolve the absolute threshold shift and the spin-coded shift vector.
analytical_resolve_delta <- function(thresholds, amount_of_SDs_perturbation,
                                     perturbation_type) {
  if (!is.numeric(amount_of_SDs_perturbation) ||
      length(amount_of_SDs_perturbation) != 1L ||
      !is.finite(amount_of_SDs_perturbation) ||
      amount_of_SDs_perturbation <= 0)
    stop("'amount_of_SDs_perturbation' must be a finite positive number.")

  sdt <- stats::sd(thresholds)
  if (is.finite(sdt) && sdt > 0) {
    delta <- amount_of_SDs_perturbation * sdt
  } else {
    warning("sd(thresholds) == 0; a raw perturbation of ",
            amount_of_SDs_perturbation, " is applied instead.")
    delta <- amount_of_SDs_perturbation
  }
  sgn <- if (perturbation_type == "aggravating") 1 else -1

  # Spin-convention threshold change per node (0.5 * {0,1}-coding change).
  list(delta = delta,
       delta_spin = rep(sgn * 0.5 * delta, length(thresholds)))
}

# Exact mean (and SD) of the number of activated nodes within outcome_idx,
# for the original network and for each one-node perturbation of source_idx.
# The 2^p microstates are enumerated ONCE per network (spin coding, stable
# log-sum-exp).  A perturbation of node i changes only the single-node field,
# so the perturbed weights are obtained by exactly reweighting the microstates
# by exp(delta_spin[i] * s_i): no re-enumeration is needed for the
# length(source_idx) interventions.
ising_exact_expectation <- function(thr, W, delta_spin, source_idx, outcome_idx) {
  p <- length(thr)
  if (length(delta_spin) == 1L) delta_spin <- rep(delta_spin, p)

  m  <- 0.5 * thr + 0.25 * (W %*% rep(1, p))  # spin thresholds
  Ws <- 0.25 * W                              # spin weights
  S  <- as.matrix(expand.grid(rep(list(c(-1, 1)), p)))
  h  <- -0.5 * rowSums((S %*% Ws) * S) - S %*% m  # Hamiltonian H(v)
  lf <- -h                                         # log P(v) = -H(v)
  lf <- lf - max(lf)                               # log-sum-exp
  f  <- exp(lf)

  n_active <- (rowSums(S[, outcome_idx, drop = FALSE]) +
                 length(outcome_idx)) / 2
  Z   <- sum(f)
  ma0 <- sum(n_active * f) / Z
  m2  <- sum(n_active^2 * f) / Z
  sd0 <- sqrt(max(m2 - ma0^2, 0))

  ma_i <- numeric(length(source_idx))
  for (k in seq_along(source_idx)) {
    i       <- source_idx[k]
    lp      <- as.vector(lf) + delta_spin[i] * S[, i]
    w       <- exp(lp - max(lp))
    ma_i[k] <- sum(n_active * w) / sum(w)
  }

  list(mean_original = ma0, sd_original = sd0, mean_nodes = ma_i)
}

# Decide whether resampling uses bootnet or the fast estimator.  The default
# (NULL) matches the point estimate: supplied parameters -> bootnet; else fast.
analytical_resolve_resample <- function(resample_default, supplied_thresholds) {
  use_bootnet <- TRUE
  if (is.null(resample_default)) {
    use_bootnet <- !is.null(supplied_thresholds)
  } else if (identical(resample_default, "fast")) {
    use_bootnet <- FALSE
  } else {
    if (!is.character(resample_default) || length(resample_default) != 1)
      stop("'resample_default' must be 'fast', a bootnet 'default' name, or NULL.")
    if (!requireNamespace("bootnet", quietly = TRUE))
      stop("'resample_default' requires the 'bootnet' package to be installed.")
    use_bootnet <- TRUE
  }

  if (use_bootnet) {
    if (!requireNamespace("bootnet", quietly = TRUE))
      stop("Resampling with bootnet requires the 'bootnet' package to be installed. ",
           "Install it, or set resample_default = 'fast'.")
    # When the point-estimate parameters came from
    # bootnet::estimateNetwork(default = "IsingFit") (the typical workflow),
    # default to the same estimator.
    if (is.null(resample_default)) resample_default <- "IsingFit"
  }

  list(use_bootnet = use_bootnet, resample_default = resample_default)
}

# Bootstrap the exact solution over rows of the original data.  Every
# repetition resamples participants with replacement, re-estimates the Ising
# network and recomputes the exact expectation.  Returns a
# nReps x (1 + length(source_idx)) matrix whose first column is "original".
bootstrap_exact_means <- function(data, delta_spin, source_idx, outcome_idx,
                                  nReps, use_bootnet, resample_default,
                                  seed, parallel, ncores) {
  one_rep <- function(b) {
    attempt <- 0L
    ok <- FALSE
    repeat {
      attempt <- attempt + 1L
      Xr <- data[sample.int(nrow(data), replace = TRUE), , drop = FALSE]
      if (use_bootnet) {
        e <- tryCatch({
          bfit <- bootnet::estimateNetwork(Xr, default = resample_default)
          list(thresholds = as.numeric(bfit$intercepts),
               graph      = bfit$graph)
        }, error = function(err) NULL)
      } else {
        e <- tryCatch(estimate_ising_fast(Xr), error = function(err) NULL)
      }
      ok <- !is.null(e) &&
            all(is.finite(e$thresholds)) &&
            all(is.finite(e$graph))
      if (ok || attempt >= 20L) break
    }
    if (!ok) return(c(original = NA_real_, rep(NA_real_, length(source_idx))))
    oo <- ising_exact_expectation(e$thresholds, e$graph, delta_spin,
                                  source_idx, outcome_idx)
    c(original = oo$mean_original, oo$mean_nodes)
  }

  if (parallel) {
    if (is.null(ncores)) ncores <- max(1L, parallel::detectCores() - 1L)
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterSetRNGStream(cl, seed)
    parallel::clusterExport(
      cl,
      varlist = c("data", "delta_spin", "source_idx", "outcome_idx",
                  "use_bootnet", "resample_default",
                  "estimate_ising_fast", "ising_exact_expectation", "one_rep"),
      envir   = environment())
    if (use_bootnet)
      parallel::clusterEvalQ(cl, requireNamespace("bootnet", quietly = TRUE))
    out <- parallel::parLapply(cl, seq_len(nReps), one_rep)
  } else {
    if (!is.null(seed)) set.seed(seed)
    out <- lapply(seq_len(nReps), one_rep)
  }

  resample_means <- do.call(rbind, out)
  if (anyNA(resample_means)) {
    keep <- stats::complete.cases(resample_means)
    dropped <- sum(!keep)
    if (dropped == nReps)
      stop("All resamples failed network estimation.")
    if (dropped > 0)
      warning(dropped, " resample(s) failed network estimation and were dropped.")
    resample_means <- resample_means[keep, , drop = FALSE]
  }
  resample_means
}

# Ranking-stability table across bootstrap resamples.  For every resample the
# nodes are ranked by absolute intervention effect (largest first); tied nodes
# share the occupied positions equally, so each column of proportions sums to
# one.  The output matches the layout of findMaxN(): percenttop_1..nRanks
# followed by repeattop_1..nRanks, and can be passed directly to plotMaxN().
rank_stability <- function(effects, nRanks) {
  B          <- nrow(effects)
  p          <- ncol(effects)
  node_names <- colnames(effects)
  nRanks     <- min(as.integer(nRanks), p)

  percent <- matrix(
    0, nrow = p, ncol = nRanks,
    dimnames = list(node_names, paste0("percenttop_", seq_len(nRanks))))
  magnitude <- abs(effects)

  for (b in seq_len(B)) {
    row    <- magnitude[b, ]
    ord    <- order(-row)
    sorted <- row[ord]
    tol    <- 1e-12 * pmax(1, sorted[-length(sorted)], sorted[-1])
    gap    <- sorted[-length(sorted)] - sorted[-1]
    # Group boundaries: a new rank group starts where the drop exceeds
    # numerical tolerance.
    new_grp    <- c(TRUE, gap > tol)
    grp_sorted <- cumsum(new_grp)
    grp_node   <- integer(p)
    grp_node[ord] <- grp_sorted

    for (j in seq_len(p)) {
      pos   <- which(grp_sorted == grp_node[j])
      first <- pos[1L]
      m     <- length(pos)
      # Tied nodes jointly occupy positions first..(first + m - 1).
      ranks <- first:(first + m - 1L)
      ranks <- ranks[ranks <= nRanks]
      if (length(ranks) > 0L)
        percent[j, ranks] <- percent[j, ranks] + 1 / m
    }
  }

  percent <- percent / B
  counts  <- percent * B
  colnames(counts) <- paste0("repeattop_", seq_len(nRanks))
  # Return a data.frame (not a matrix) so columns can be selected by name,
  # matching the layout of findMaxN().
  cbind(as.data.frame(percent), as.data.frame(counts))
}

# Assemble the bootstrap inference table, the ranking-stability table and the
# stability plot shared by the two stability functions.
analytical_stability_core <- function(resample_means, effects, ma_nodes,
                                      node_names, sd_original, nRanks,
                                      method) {
  B <- nrow(resample_means)
  resample_effects <- resample_means[, -1, drop = FALSE] -
                      resample_means[, "original"]
  colnames(resample_effects) <- node_names

  se_effect       <- apply(resample_effects, 2, stats::sd, na.rm = TRUE)
  ci_effect_lower <- apply(resample_effects, 2, stats::quantile,
                           probs = 0.025, na.rm = TRUE)
  ci_effect_upper <- apply(resample_effects, 2, stats::quantile,
                           probs = 0.975, na.rm = TRUE)
  ci_lower        <- apply(resample_means, 2, stats::quantile,
                           probs = 0.025, na.rm = TRUE)
  ci_upper        <- apply(resample_means, 2, stats::quantile,
                           probs = 0.975, na.rm = TRUE)

  # Approximate centered-bootstrap tail area for H0: effect = 0.  This is not
  # a pairwise target-comparison test, and its calibration depends on the
  # network estimator and bootstrap approximation.
  pvals <- sapply(seq_along(node_names), function(i) {
    d <- resample_effects[, i]
    (1 + sum(abs(d - effects[i]) >= abs(effects[i]))) / (B + 1)
  })

  maxn   <- rank_stability(resample_effects, nRanks)
  nRanks <- ncol(maxn) / 2L

  stat <- data.frame(
    node             = node_names,
    mean_activation  = ma_nodes,
    ci_lower         = unname(ci_lower[node_names]),
    ci_upper         = unname(ci_upper[node_names]),
    effect           = effects,
    glass_d          = effects / sd_original,
    se_effect        = se_effect,
    ci_effect_lower  = ci_effect_lower,
    ci_effect_upper  = ci_effect_upper,
    p                = pvals,
    p.adjust         = stats::p.adjust(pvals, method = method),
    row.names        = node_names,
    stringsAsFactors = FALSE
  )
  for (k in seq_len(nRanks))
    stat[[paste0("percenttop_", k)]] <- maxn[[paste0("percenttop_", k)]]
  stat$stars <- ifelse(is.na(stat$p.adjust), "",
                       ifelse(stat$p.adjust < 0.001, "***",
                              ifelse(stat$p.adjust < 0.01, "**",
                                     ifelse(stat$p.adjust < 0.05, "*", ""))))

  gp <- plotMaxN(
    maxn, low = 1, high = nRanks,
    xlab = "Cumulative occurrence percentage across bootstrap resamples")

  list(stat             = stat,
       resample_effects = resample_effects,
       maxn             = maxn,
       plot             = gp,
       B                = B)
}

# Exact mean-activation line plot without confidence intervals.  Shared by
# analyticalNIRAtest() and analyticalBridgeNIRAtest().
analytical_exact_plot <- function(ma_original, ma_nodes, node_names,
                                  perturbation_type, ylab) {
  plot_data <- data.frame(
    condition       = c("original", node_names),
    mean_activation = c(ma_original, ma_nodes),
    stringsAsFactors = FALSE
  )
  if (perturbation_type == "aggravating") {
    ord <- order(plot_data$condition != "original",
                 -plot_data$mean_activation)
  } else {
    ord <- order(plot_data$condition != "original",
                 plot_data$mean_activation)
  }
  plot_data   <- plot_data[ord, , drop = FALSE]
  plot_data$id <- seq_len(nrow(plot_data)) - 1

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$id, y = .data$mean_activation, group = 1)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(size = 3, color = "tomato") +
    ggplot2::scale_x_continuous(breaks = plot_data$id,
                                labels = plot_data$condition) +
    ggplot2::labs(x = "Condition", y = ylab, title = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x       = ggplot2::element_text(face = "bold", size = 12),
      axis.title.y       = ggplot2::element_text(face = "bold", size = 12),
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1,
                                                 face = "bold", size = 10),
      axis.text.y        = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.major.x = ggplot2::element_blank())
}
