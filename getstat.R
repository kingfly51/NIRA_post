#' @param gs_sumIsingSamplesLong The output result of the prepareDFforPlottingAndANOVA function
#' @param method method of p.adjust,c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")default:bonferroni
#' @examples getstat(gs_sumIsingSamplesLong)
#' @author Written by Wang Fei 2025/03/22
#' NIRA
#' State Key Laboratory of Cognitive Neuroscience and Learning,Beijing Normal University
#' @email bjnwangfei0501@outlook.com

getstat <- function(gs_sumIsingSamplesLong,method="bonferroni") {
  
    statresult <- list(
      stat = list(),
      plot_data = list()
    )
    
    wide_data <- gs_sumIsingSamplesLong %>%
        group_by(sample) %>%
        mutate(row = row_number()) %>%
        tidyr::pivot_wider(names_from = sample, values_from = sumscore) %>%
        select(-row)%>%
        mutate(across(everything(), as.numeric)) 
    
    categorie <- unique(gs_sumIsingSamplesLong[[2]])
    original_col <- as.matrix(wide_data[,"original"])
    other_cols <- as.matrix(wide_data[, !categorie %in% "original"])

    set.seed(2025)
    stat_list<-as.data.frame(NULL)
    for (i in(1:ncol(other_cols))){
      perlist<-permutation_test(original_col, other_cols[, i],5000)
      stat_list[i,1]<-perlist$mean_other
      stat_list[i,2]<-perlist$sd_other
      stat_list[i,3]<-perlist$se_other
      stat_list[i,4]<-perlist$ci_other
      stat_list[i,5]<-perlist$ciLower_other
      stat_list[i,6]<-perlist$ciUpper_other
      stat_list[i,7]<-perlist$cohens_d
      stat_list[i,8]<-perlist$p_value
    }
    rownames(stat_list) = colnames(other_cols)
    colnames(stat_list) =c("mean_other","sd_other","se_other","ci_other","ciLower_other","ciUpper_other","cohen_d","p")
    stat_list$p.adjust = p.adjust(stat_list$p,method = method)
    
    dat1<-stat_list[,c(1,5,6)]
    colnames(dat1) <- c("mean","ciLower","ciUpper")
    dat2<-cbind(perlist$mean_original,perlist$ciLower_original,perlist$ciUpper_original)
    colnames(dat2) <- c("mean","ciLower","ciUpper")
    rownames(dat2) <- "original"
    dat3<-rbind(dat2,dat1)
    dat3$id <-0:(nrow(dat3)-1)
    dat3$node <- rownames(dat3)

    statresult$stat <- stat_list
    statresult$plot_data <- dat3
    return(statresult)
}


