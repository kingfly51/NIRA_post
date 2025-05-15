#' @param statresult The output of the getstat function
#' @param perturbation_type a string specifying a perturbation direction. Choose between "aggravating" (+) and "alleviating" (-). Should be equal to the argument passed to simulateReponses().
#' @examples plotmeanNIRA(statresult,"aggravating")
#' @author Written by Wang Fei 2025/03/22
#' NIRA
#' State Key Laboratory of Cognitive Neuroscience and Learning,Beijing Normal University
#' @email bjnwangfei0501@outlook.com

plotmeanNIRA <- function(statresult,perturbation_type) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
    library(dplyr)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
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
        mutate(node = rownames(stat_result$stat)),  
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






