#' Plot Node Means with Confidence Intervals and Significance Stars
#'
#' @description
#' Produces a publication-ready line plot showing the mean sum score of
#' each network condition (original + all intervention nodes) with 95 \%
#' confidence intervals, sorted by effect direction, and annotated with
#' significance stars from the multiple-comparison-corrected p-values
#' computed by \code{\link{permutationNIRAtest}}.
#'
#' Nodes are sorted so that the intervention with the strongest effect
#' appears at one end of the x-axis:
#' \itemize{
#'   \item \code{"alleviating"}: sorted in \strong{ascending} order of
#'     mean (lowest mean = most effective alleviating intervention).
#'   \item \code{"aggravating"}: sorted in \strong{descending} order of
#'     mean (highest mean = strongest aggravating effect).
#' }
#' The unperturbed \code{"original"} condition is always placed first
#' (leftmost).
#'
#' Significance stars reflect the adjusted p-values from
#' \code{permutationNIRAtest}:
#' \itemize{
#'   \item \code{***}: \eqn{p_\text{adj} < 0.001}
#'   \item \code{**}: \eqn{p_\text{adj} < 0.01}
#'   \item \code{*}: \eqn{p_\text{adj} < 0.05}
#'   \item (blank): \eqn{p_\text{adj} \geq 0.05}
#' }
#'
#' @param statresult A \code{list} produced by
#'   \code{\link{permutationNIRAtest}}.  Must contain:
#'   \describe{
#'     \item{\code{plot_data}}{A \code{data.frame} with columns
#'       \code{mean}, \code{ciLower}, \code{ciUpper}, \code{node}.}
#'     \item{\code{stat}}{A \code{data.frame} (rows = nodes) with a
#'       \code{p.adjust} column.}
#'   }
#' @param perturbation_type Character string: \code{"aggravating"} or
#'   \code{"alleviating"}.  Must match the type passed to
#'   \code{nodeIdentifyR::simulateResponses} and
#'   \code{nodeIdentifyR::plotSumScores}.
#'
#' @return
#' A \code{\link[ggplot2]{ggplot}} object that can be:
#' \itemize{
#'   \item Displayed immediately by printing.
#'   \item Saved with \code{ggplot2::ggsave()}.
#'   \item Further customised with additional \pkg{ggplot2} layers.
#' }
#'
#' The plot contains:
#' \describe{
#'   \item{Blue line}{Connects mean scores across conditions.}
#'   \item{Red points}{Mean score for each condition.}
#'   \item{Grey error bars}{95 \% confidence intervals.}
#'   \item{Black stars}{Significance annotation above each point.}
#' }
#'
#' @section Method Recommendations:
#' \tabular{ll}{
#'   \strong{Customisation} \tab \strong{Code} \cr
#'   Change title \tab
#'     \code{+ ggplot2::labs(title = "My Title")} \cr
#'   Save at high resolution \tab
#'     \code{ggplot2::ggsave("fig.pdf", width=8, height=5)} \cr
#'   Rotate x labels \tab
#'     \code{+ ggplot2::theme(axis.text.x = element_text(angle=90))} \cr
#'   Report in paper \tab
#'     Use the \code{p.adjust} column from \code{stat_result$stat} \cr
#' }
#'
#' @examples
#' ## Example 1: standard alleviating intervention plot
#' \dontrun{
#' data("single_gds")
#' gs_fit  <- bootnet::estimateNetwork(single_gds,
#'              default = "IsingFit", rule = "AND")
#' gs_samp <- nodeIdentifyR::simulateResponses(
#'              gs_fit$graph, gs_fit$intercepts, "alleviating", 2)
#' gs_long <- nodeIdentifyR::prepareDFforPlottingAndANOVA(
#'              nodeIdentifyR::calculateSumScores(gs_samp))
#'
#' stat_result <- permutationNIRAtest(gs_long, method = "bonferroni")
#' plotMeanNIRA(stat_result, "alleviating")
#' }
#'
#' ## Example 2: aggravating intervention
#' \dontrun{
#' gs_samp2 <- nodeIdentifyR::simulateResponses(
#'               gs_fit$graph, gs_fit$intercepts, "aggravating", 2)
#' gs_long2 <- nodeIdentifyR::prepareDFforPlottingAndANOVA(
#'               nodeIdentifyR::calculateSumScores(gs_samp2))
#' stat2 <- permutationNIRAtest(gs_long2, method = "BH")
#' plotMeanNIRA(stat2, "aggravating")
#' }
#'
#' ## Example 3: customise and save
#' \dontrun{
#' library(ggplot2)
#' p <- plotMeanNIRA(stat_result, "alleviating") +
#'   labs(title = "GAD-7: Alleviating Intervention Effects",
#'        subtitle = "Bonferroni-corrected permutation tests") +
#'   theme(plot.title = element_text(size = 14, face = "bold"))
#' ggsave("Figure_3.pdf", p, width = 8, height = 5, dpi = 300)
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{permutationNIRAtest}} — produces the required
#'     \code{statresult} input.
#'   \item \code{\link{permutationTest}} — the single-comparison engine
#'     underlying the significance stars.
#'   \item \code{\link[ggplot2]{ggsave}} — save the returned plot.
#'   \item \code{\link[nodeIdentifyR]{plotSumScores}} — the simpler
#'     nodeIdentifyR visualisation that this plot replaces with formal
#'     inference.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_errorbar
#'   geom_text scale_x_continuous labs theme_minimal theme
#'   element_text element_blank
#' @importFrom dplyr mutate arrange left_join case_when desc
#' @importFrom magrittr %>%
#' @export
plotMeanNIRA <- function(statresult, perturbation_type) {

  if (!perturbation_type %in% c("aggravating", "alleviating"))
    stop("perturbation_type must be 'aggravating' or 'alleviating'.")
  if (!all(c("plot_data", "stat") %in% names(statresult)))
    stop("statresult must contain 'plot_data' and 'stat' components.")

  # Sort by direction
  od <- if (perturbation_type == "aggravating")
    statresult$plot_data %>%
      dplyr::mutate(nt = ifelse(node == "original", "original", "other")) %>%
      dplyr::arrange(nt, dplyr::desc(mean)) %>%
      dplyr::mutate(node = factor(node, levels = unique(node)))
  else
    statresult$plot_data %>%
      dplyr::mutate(nt = ifelse(node == "original", "original", "other")) %>%
      dplyr::arrange(nt, mean) %>%
      dplyr::mutate(node = factor(node, levels = unique(node)))
  od$newid <- 0:(nrow(od) - 1)

  # Join significance stars
  pd <- od %>%
    dplyr::left_join(
      statresult$stat %>%
        dplyr::mutate(node = rownames(statresult$stat)),
      by = "node") %>%
    dplyr::mutate(stars = dplyr::case_when(
      p.adjust < 0.001 ~ "***",
      p.adjust < 0.01  ~ "**",
      p.adjust < 0.05  ~ "*",
      TRUE             ~ ""))

  ggplot2::ggplot(pd, ggplot2::aes(x = newid, y = mean, group = 1)) +
    ggplot2::geom_line(color = "steelblue", linewidth = 1) +
    ggplot2::geom_point(size = 3, color = "tomato") +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = ciLower, ymax = ciUpper),
      width = 0.2, color = "gray40") +
    ggplot2::geom_text(
      ggplot2::aes(label = stars, y = mean + 1.2 * (mean - ciLower)),
      color = "black", size = 5, vjust = 0.5) +
    ggplot2::scale_x_continuous(breaks = pd$newid, labels = pd$node) +
    ggplot2::labs(x = "Node", y = "Mean Sum Score", title = "") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x       = ggplot2::element_text(face = "bold", size = 12),
      axis.title.y       = ggplot2::element_text(face = "bold", size = 12),
      axis.text.x        = ggplot2::element_text(angle = 45, hjust = 1,
                                                  face = "bold", size = 10),
      axis.text.y        = ggplot2::element_text(face = "bold", size = 10),
      panel.grid.major.x = ggplot2::element_blank())
}