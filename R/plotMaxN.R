#' Stacked Bar Plot of Node Ranking Stability Across Simulation Repetitions
#'
#' @description
#' Visualises the output of \code{\link{findMaxN}} as a horizontal stacked
#' bar plot.  Each bar represents one intervention node; stacked segments
#' show what proportion of simulation repetitions that node occupied each
#' rank position (from rank \code{low} to rank \code{high}).
#'
#' The plot is designed to communicate \emph{stability} at a glance:
#' \itemize{
#'   \item A node with a long dark segment (rank 1) that dwarfs all others
#'     is a \strong{highly stable} intervention target — it consistently
#'     produced the largest difference from the original simulated network
#'     across independent simulation runs.
#'   \item A node whose colour segments are spread roughly equally across
#'     ranks is \strong{unstable} — its "best" position varies widely
#'     across runs.
#' }
#'
#' Bars are ordered from the node with the lowest total cumulative
#' percentage (top of y-axis) to the highest (bottom), so the most stable
#' node appears at the bottom.  Colours follow the Nature Publishing Group
#' (NPG) palette from \pkg{ggsci}.
#'
#' @param maxn A \code{data.frame} produced by \code{\link{findMaxN}}.
#'   Must contain columns named \code{percenttop_1}, \code{percenttop_2},
#'   …, up to the \code{n} used in \code{findMaxN}.
#' @param low Integer.  First (lowest) rank position to display.
#'   Default: \code{1}.
#' @param high Integer.  Last (highest) rank position to display.
#'   Must satisfy \code{high >= low} and
#'   \code{high <= n} where \code{n} was passed to \code{findMaxN}.
#'   Default: \code{3}.  The tutorial paper uses \code{low = 1, high = 7}.
#'
#' @return
#' A \code{\link[ggplot2]{ggplot}} object.  It can be:
#' \itemize{
#'   \item Displayed by printing.
#'   \item Saved with \code{ggplot2::ggsave()}.
#'   \item Extended with additional \pkg{ggplot2} layers.
#' }
#'
#' \strong{Plot elements:}
#' \describe{
#'   \item{Horizontal bars}{Each bar represents one node; total bar
#'     length = sum of \code{percenttop_low} through
#'     \code{percenttop_high}.}
#'   \item{Stacked segments}{Each segment's colour identifies the rank
#'     position; length = percentage of repetitions at that rank.}
#'   \item{X-axis}{Labelled in percentage (0 \%–100 \%).}
#'   \item{Y-axis}{Node names, ordered lexicographically: primarily by
#'     \code{percenttop_1} (descending), ties broken by
#'     \code{percenttop_2}, then \code{percenttop_3}, etc.}
#'   \item{Legend}{Maps colours to rank labels.}
#' }
#'
#' @section Method Recommendations:
#' \tabular{lll}{
#'   \strong{Network size} \tab \strong{Recommended low:high} \tab
#'     \strong{Rationale} \cr
#'   ≤ 5 nodes \tab 1 : p \tab Show all possible ranks \cr
#'   6–10 nodes \tab 1 : 5 \tab Readable legend; top-5 captures most signal \cr
#'   > 10 nodes \tab 1 : 3 \tab Top-3 is sufficient for reporting \cr
#' }
#'
#' If rank segments overlap badly, reduce \code{high}.  If the plot looks
#' empty, increase \code{nReps} in \code{\link{stabilityNIRAtest}} to produce
#' more reliable frequency estimates.
#'
#' @examples
#' ## Example 1: standard usage
#' \dontrun{
#' data("single_gds")
#' gs_fit <- bootnet::estimateNetwork(single_gds,
#'             default = "IsingFit", rule = "AND")
#' simResult <- stabilityNIRAtest(
#'   gs_fit$graph, gs_fit$intercepts,
#'   "alleviating", 2, nReps = 100)
#'
#' node_ranks <- findMaxN(simResult, n = 7)
#' plotMaxN(node_ranks, low = 1, high = 7)
#' }
#'
#' ## Example 2: show only top 3 ranks
#' \dontrun{
#' plotMaxN(node_ranks, low = 1, high = 3)
#' }
#'
#' ## Example 3: customise the plot
#' \dontrun{
#' library(ggplot2)
#' p <- plotMaxN(node_ranks, low = 1, high = 7) +
#'   labs(title = "GAD-7 Node Ranking Stability",
#'        subtitle = "1 000 independent simulation repetitions") +
#'   theme(plot.title = element_text(face = "bold", size = 14))
#' ggsave("Figure_4.pdf", p, width = 7, height = 5, dpi = 300)
#' }
#'
#' @author
#' Fei Wang \email{bjnwangfei0501@@outlook.com}
#'
#' @seealso
#' \itemize{
#'   \item \code{\link{findMaxN}} — produces the required input.
#'   \item \code{\link{stabilityNIRAtest}} — generates the simulation
#'     repetitions that \code{findMaxN} summarises.
#'   \item \code{\link[ggsci]{scale_fill_npg}} — the NPG colour palette
#'     applied to the stacked segments.
#'   \item \code{\link[ggplot2]{ggsave}} — save the returned plot.
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_bar scale_x_continuous labs
#'   theme_minimal theme element_text element_blank
#' @importFrom tidyr pivot_longer
#' @importFrom ggsci scale_fill_npg
#' @export
plotMaxN <- function(maxn, low = 1, high = 3) {

  if (!is.data.frame(maxn) && !is.matrix(maxn))
    stop("maxn must be the data.frame output of findMaxN().")

  pct_cols <- grep("percenttop", colnames(maxn), value = TRUE)[low:high]
  if (length(pct_cols) == 0)
    stop("No 'percenttop_*' columns found for the requested low:high range.")

  # Lexicographic ordering: sort by percenttop_low first, then percenttop_(low+1),
  # and so on.  This ensures the node most frequently ranked #1 sits at the
  # bottom (most prominent), with ties broken by rank 2, then rank 3, etc.
  sort_data <- maxn[, pct_cols, drop = FALSE]
  # lapply over rev(pct_cols): least-priority key first for do.call(order)
  sort_args <- lapply(rev(pct_cols), function(col) sort_data[[col]])
  ord       <- rownames(maxn)[do.call(order, c(sort_args, list(decreasing = FALSE)))]
  sub      <- maxn[ord, pct_cols, drop = FALSE]
  sub$node <- factor(rownames(sub), levels = ord)

  dl          <- tidyr::pivot_longer(sub, cols = -node,
                                     names_to = "Rank", values_to = "Pct")
  dl$Pct      <- dl$Pct * 100

  ggplot2::ggplot(dl, ggplot2::aes(x = Pct, y = node, fill = Rank)) +
    ggplot2::geom_bar(stat = "identity", position = "stack") +
    ggsci::scale_fill_npg() +
    ggplot2::scale_x_continuous(
      labels = function(x) paste0(x, "%"),
      expand = c(0, 0)) +
    ggplot2::labs(
      x    = "Cumulative occurrence percentage across simulation repetitions",
      y    = "Node",
      fill = "Rank") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      text               = ggplot2::element_text(size = 12),
      legend.position    = "right",
      panel.grid.major.y = ggplot2::element_blank())
}