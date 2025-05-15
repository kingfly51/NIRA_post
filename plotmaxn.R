#' @param maxn The output result of the findmaxn function(m*n matrix,m:subject,n:variables)
#' @param low Number, the range cannot exceed the number of columns in maxn, and the lowest column of the data to be plotted is set to 1 by default
#' @param high Number, the range cannot exceed the number of columns in maxn, and the highest column of the data to be plotted is set to 3 by default
#' @examples plotmaxn (maxn,low=1,high=3)
#' @author Written by Wang Fei 2025/03/23
#' NIRA
#' State Key Laboratory of Cognitive Neuroscience and Learning,Beijing Normal University
#' @email bjnwangfei0501@outlook.com


plotmaxn <- function(maxn,low=1,high=3) {
  if (!require("ggplot2")) install.packages("ggplot2")
  if (!require("ggsci")) install.packages("ggsci")
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    install.packages("tidyr")}
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

