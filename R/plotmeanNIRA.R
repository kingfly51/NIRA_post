#' Visualize Mean Values with Confidence Intervals from NIRA Analysis
#'
#' Creates a line plot showing mean values with confidence intervals for original and
#' perturbed networks, annotated with statistical significance.
#'
#' @param statresult A list output from [`getstat()`] containing:
#'   * `plot_data`: Formatted data for plotting (must contain columns: `mean`, `ciLower`, `ciUpper`, `node`)
#'   * `stat`: Statistical results with `p.adjust` column
#' @param perturbation_type Intervention direction:
#'   * `"aggravating"` - Sorts nodes by descending mean values
#'   * `"alleviating"` - Sorts nodes by ascending mean values
#'
#' @return
#' A [ggplot2][ggplot2::ggplot2-package] object showing:
#'
#' **Plot Elements**:
#' * Line plot connecting node means
#' * Points indicating mean values
#' * Error bars for 95% confidence intervals
#' * Significance markers:
#'   - `***` p < 0.001
#'   - `**` p < 0.01
#'   - `*` p < 0.05
#'
#' @examples
#' ```r
#' # After running getstat()
#' library(ggplot2)
#'
#' # For aggravating interventions
#' p <- plotmeanNIRA(statresult, "aggravating")
#' p + labs(title = "Aggravating Intervention Effects")
#'
#' # For alleviating interventions
#' plotmeanNIRA(statresult, "alleviating")
#' ```
#'
#' @author
#' Wang Fei, State Key Laboratory of Cognitive Neuroscience and Learning,
#' Beijing Normal University (2025-03-22)
#'
#' @section Technical Details:
#' The plot automatically:
#' 1. Orders nodes by intervention effect size
#' 2. Adds significance stars based on adjusted p-values
#' 3. Formats axes for readability (45° x-axis labels)
#'
#' @seealso
#' * [`getstat()`] for generating the input data
#' * [`ggplot2::theme_minimal()`] for plot styling
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_errorbar
#' @importFrom ggplot2 geom_text scale_x_continuous labs theme_minimal
#' @importFrom ggplot2 element_text element_blank
#' @importFrom dplyr %>% mutate arrange left_join case_when if_else desc
#' @importFrom magrittr %>%
#' @export

plotmeanNIRA <- function(statresult, perturbation_type) {
  # Validate inputs
   if (!perturbation_type %in% c("aggravating", "alleviating")) {
     stop('perturbation_type must be either "aggravating" or "alleviating"')
   }

   if (!all(c("plot_data", "stat") %in% names(statresult))) {
     stop("statresult must contain 'plot_data' and 'stat' components")
   }

   if (perturbation_type=="aggravating") {
    ordered_data <- statresult$plot_data %>%
      mutate(node_type = ifelse(node == "original", "original", "other")) %>%
      arrange(node_type, desc(mean)) %>%
      mutate(node = factor(node, levels = unique(node)))
    ordered_data$newid<- 0:(nrow(ordered_data)-1)
  } else if (perturbation_type=="alleviating") {
    ordered_data <- statresult$plot_data %>%
      mutate(node_type = ifelse(node == "original", "original", "other")) %>%
      arrange(node_type, mean) %>%
      mutate(node = factor(node, levels = unique(node)))
    ordered_data$newid<- 0:(nrow(ordered_data)-1)
  }

  plot_data <- ordered_data %>%
    left_join(
      statresult$stat %>%
        mutate(node = rownames(statresult$stat)),
      by = "node"
    ) %>%
    mutate(
      stars = case_when(
        p.adjust < 0.001 ~ "***",
        p.adjust < 0.01 ~ "**",
        p.adjust < 0.05 ~ "*",
        TRUE ~ ""
      )
    )

  ggplot2::ggplot(plot_data, aes(x = newid, y = mean, group = 1)) +
    geom_line(color = "steelblue", linewidth = 1) +
    geom_point(size = 3, color = "tomato") +
    geom_errorbar(
      aes(ymin = ciLower, ymax = ciUpper),
      width = 0.2,
      color = "gray40"
    ) +
    geom_text(
      aes(label = stars, y = mean+1.2*(mean-ciLower)),
      color = "black",
      size = 5,
      vjust = 0.5
    ) +
    scale_x_continuous(
      breaks = plot_data$newid,
      labels = plot_data$node
    ) +
    labs(
      x = "Node",
      y = "Mean Value",
      title = ""
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.7, face = "bold", size = 14),
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 10),
      axis.text.y = element_text(face = "bold", size = 10),
      panel.grid.major.x = element_blank()
    )
}





