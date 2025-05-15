#' @param Bootresult The output result of the BootNIRAresult function
#' @param n Extract the frequency and percentage of occurrences for each node corresponding to the 1st, 2nd, ..., nth largest differences compared to the original simulated sample across nboots networks.Default is 3
#' @examples findmaxn(Bootresult,n=3)
#' @author Written by Wang Fei 2025/03/23
#' NIRA
#' State Key Laboratory of Cognitive Neuroscience and Learning,Beijing Normal University
#' @email bjnwangfei0501@outlook.com

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
