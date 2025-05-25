#' Identify Most Frequently Perturbed Nodes in Bootstrap Analysis
#'
#' Calculates the frequency and percentage of nodes showing the largest differences
#' compared to the original simulated sample across bootstrap networks.
#'
#' @param Bootresult List output from [`BootNIRAresult()`] containing:
#'   * `mean`: Matrix with columns:
#'     - `"original"`: Original network means
#'     - Other columns: Perturbed network means
#' @param n Integer specifying how many top differences to consider (default: 3).
#'   Returns frequencies for top 1 through top n differences.
#'
#' @return
#' A data frame with two sections per rank (1 to n):
#'
#' **Frequency Analysis**:
#' - `repeattop_*`: Absolute counts of nodes appearing in top n differences
#' - `percenttop_*`: Percentage occurrence (counts/total bootstrap samples)
#'
#' Rows represent network nodes, columns show statistics for each rank.
#'
#' @examples
#' ```r
#' # After running BootNIRAresult()
#' top_nodes <- findmaxn(Bootresult, n = 3)
#'
#' # Get most frequently perturbed nodes
#' head(top_nodes[order(top_nodes$percenttop_1, decreasing = TRUE), ])
#' ```
#' @author
#' Wang Fei, State Key Laboratory of Cognitive Neuroscience and Learning,
#' Beijing Normal University (2025-03-22)
#'
#' @section Algorithm Details:
#' 1. Computes absolute differences between original and perturbed means
#' 2. For each bootstrap sample, identifies top n nodes with largest differences
#' 3. Aggregates frequencies across all samples
#' 4. Calculates percentage occurrence
#'
#' @seealso
#' * [`BootNIRAresult()`] for required input data generation
#' * [`apply()`] for matrix operations core to the algorithm
#'
#' @importFrom stats setNames
#' @importFrom utils head
#' @export

findmaxn <- function(Bootresult,n=3) {
  original_col <- Bootresult$mean[,"original"]
  other_cols <-  Bootresult$mean[, !colnames(Bootresult$mean) %in% "original"]
  meandata <- apply(other_cols,2,function(x){abs(original_col-x)})

  get_top_n_colnames <- function(row) {
    sorted_indices <- order(row, decreasing = TRUE)
    colnames(meandata)[sorted_indices[1:n]]
  }

  top_n_cols <- apply(meandata, 1, get_top_n_colnames)
  top_n_cols_df <- as.data.frame(t(top_n_cols))
  colnames(top_n_cols_df) <- paste0("repeattop_",1:n,seq="")

  df <- data.frame(matrix(0,nrow = ncol(other_cols), ncol = n))
  rownames(df) = colnames(other_cols)
  colnames(df) = colnames(top_n_cols_df)
  for (i in 1:n){
    top_i <- top_n_cols_df[, i]
    counts <- table(top_i)
    for (var in names(counts)) {
      df[var, i] <- counts[var]
    }
  }
  df_2 <- df
  colnames(df_2)<-paste0("percenttop_",1:n,seq="")
  df_2<-df_2/nrow(Bootresult$mean)
  maxn<-cbind(df_2,df)
  return(maxn)
}
