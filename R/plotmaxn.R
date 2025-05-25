#' Visualize Node Occurrence Frequencies
#'
#' @description
#' Creates a line plot showing the frequency distribution of nodes appearing in top-n
#' most perturbed positions across bootstrap replications. Uses Nature Publishing Group
#' color palette for categorical coloring.
#'
#' @param maxn A matrix or data frame output from [`findmaxn()`] where:
#'   - Rows represent network nodes
#'   - Columns contain frequency counts (typically `percenttop_1` to `percenttop_n`)
#' @param low Integer specifying the first rank to include (default: 1)
#' @param high Integer specifying the last rank to include (default: 3)
#'   - Must satisfy: `1 ≤ low ≤ high ≤ ncol(maxn)`
#'
#' @return
#' A [ggplot2][ggplot2::ggplot2-package] object with:
#'
#' **Graphical Elements**:
#' - **Lines**: Connect frequency values across ranks for each node
#' - **Points**: Mark individual frequency values
#' - **Color Mapping**: Differentiates ranks using NPG palette
#' - **Axes**:
#'   - Y-axis: Node names (flipped coordinates)
#'   - X-axis: Occurrence percentage (0-100%)
#' - **Theme**: Minimal theme with right-aligned legend
#'
#' @examples
#' ```r
#' # Basic usage
#' library(ggplot2)
#' plot <- plotmaxn(top_nodes_results)
#'
#' # Customized plot
#' plot +
#'   labs(title = "Top 3 Most Perturbed Nodes",
#'        subtitle = "Bootstrap Frequency Distribution") +
#'   theme(legend.position = "bottom")
#' ```
#' @author
#' Wang Fei, State Key Laboratory of Cognitive Neuroscience and Learning,
#' Beijing Normal University (2025-03-22)
#'
#' @section Visual Encoding:
#' 1. **Position**:
#'    - Vertical: Nodes ordered by input data
#'    - Horizontal: Occurrence percentage
#' 2. **Color**: Represents different top-n ranks
#' 3. **Line Weight**: 1pt for clear visibility
#' 4. **Point Size**: 3pt for visual emphasis
#'
#' @seealso
#' * [`findmaxn()`] for required input data preparation
#' * [`ggsci::scale_color_npg()`] for color palette details
#' * [`ggplot2::coord_flip()`] for vertical label orientation
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_x_continuous
#' @importFrom ggplot2 labs theme_minimal theme coord_flip
#' @importFrom tidyr pivot_longer
#' @importFrom ggsci scale_color_npg
#' @export

plotmaxn <- function(maxn,low=1,high=3) {
  maxn_wide<-maxn[,low:high]
  #ID<-data.frame(id=rownames(maxn_wide))
  ID <- data.frame(id=1:nrow(maxn_wide))
  maxn_wide <- cbind(maxn_wide,ID)
  df_long <-  tidyr::pivot_longer(maxn_wide,
      cols = -id,
      names_to = "Node",
      values_to = "Value")

  ggplot2::ggplot(df_long, aes(x = id, y = Value, color = Node)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    ggsci::scale_color_npg() +
    scale_x_continuous(
      breaks = 1:length(rownames(maxn_wide)),
      labels = rownames(maxn_wide)
    )+
    labs(
      x = "Node",
      y = "Node occurrence percentage in total replication",
      color = "Category"
    ) +
    theme_minimal() +
    theme(
      text = element_text(size = 12),
      legend.position = "right"
    )+
    coord_flip()
}
