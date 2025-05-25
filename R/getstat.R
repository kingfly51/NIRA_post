#' Calculate Statistical Comparisons Between Original and Perturbed Networks
#'
#' @description
#' Performs permutation tests between original network scores and perturbed network scores,
#' with multiple comparison correction. Part of the NodeIdentifyR Algorithm (NIRA) framework.
#'
#' @param gs_sumIsingSamplesLong A long-format data frame output from 
#'   [`prepareDFforPlottingAndANOVA()`], containing sum scores for original and 
#'   perturbed networks. Must contain columns:
#'   * `sample`: Network condition identifier
#'   * `sumscore`: Computed sum scores
#' @param method Multiple testing correction method. One of:  
#'   `"holm"`, `"hochberg"`, `"hommel"`, `"bonferroni"`, `"BH"`, `"BY"`, `"fdr"`, `"none"`  
#'   (default: `"bonferroni"`)
#'
#' @return 
#' A list with two components:
#' 
#' **1. `stat`**  
#' Data frame containing:
#' * `mean_other` - Mean of perturbed network  
#' * `sd_other` - Standard deviation  
#' * `se_other` - Standard error  
#' * `ci_*` - Confidence interval columns  
#' * `cohen_d` - Cohen's d effect size  
#' * `p` - Raw p-value  
#' * `p.adjust` - Adjusted p-value  
#' 
#' **2. `plot_data`**  
#' Visualization-ready data with columns:
#' * `mean` - Mean score  
#' * `ciLower`/`ciUpper` - 95% CI bounds  
#' * `id` - Numeric ID  
#' * `node` - Node name  
#'
#' @examples
#' ```r
#' # After prepareDFforPlottingAndANOVA()
#' results <- getstat(gs_sumIsingSamplesLong, method = "BH")
#' 
#' # Access results
#' head(results$stat)
#' plot_data <- results$plot_data
#' ```
#' 
#' @author
#' Wang Fei, State Key Laboratory of Cognitive Neuroscience and Learning,  
#' Beijing Normal University (2025-03-22)
#'
#' @seealso
#' * [`prepareDFforPlottingAndANOVA()`] for data preparation
#' * [`p.adjust()`] for multiple testing correction
#' * [`permutation_test()`] for core statistical computation
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate select across row_number
#' @importFrom tidyr pivot_wider
#' @importFrom stats p.adjust sd
#' @importFrom rlang .data
#' @export

getstat <- function(gs_sumIsingSamplesLong,method="bonferroni") {
   # --- Input validation ---
    if (!is.data.frame(gs_sumIsingSamplesLong)) {
      stop("Input must be a data frame from prepareDFforPlottingAndANOVA()")
    }
   
    if (!all(c("sample", "sumscore") %in% colnames(gs_sumIsingSamplesLong))) {
      stop("Input must contain 'sample' and 'sumscore' columns")
    }
   
    valid_methods <- c("holm", "hochberg", "hommel", "bonferroni", 
                      "BH", "BY", "fdr", "none")
    if (!method %in% valid_methods) {
      stop("Invalid method. Choose from: ", paste(valid_methods, collapse = ", "))
    }

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
