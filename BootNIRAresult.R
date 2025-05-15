#' @param edge_weights Edge weight matrix in Ising network
#' @param thresholds Threshold parameter vector composed of intercept terms of all nodes in Ising network
#' @param perturbation_type a string specifying a perturbation direction. Choose between "aggravating" (+) and "alleviating" (-)
#' @param amount_of_SDs_perturbation an integer specifying with how many standard deviations of the threshold distribution to perturbate the threshold values
#' @param nboots Number of resampling,eg:1000
#' @param parallel Logical scalar,eg:TRUE, indicating whether to execute parallel running,default is False
#' @param ncores Number,eg:6, choose the corresponding number of cores based on the performance of your computer. If parallel running is set without adding ncores, the maximum number of cores allowed will be set by default, and ncores will be set to NULL by default
#' @examples BootNIRAresult(edgeWeightMatrix,thresholdVector,"aggravating",2,1000) for Serial
#' @examples Bootresult<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"aggravating",2,1000,parallel = TRUE, n_cores = 6) for parallel
#' @author Written by Wang Fei 2025/03/23
#' NIRA
#' State Key Laboratory of Cognitive Neuroscience and Learning,Beijing Normal University
#' @email bjnwangfei0501@outlook.com


BootNIRAresult <- function(edge_weights, thresholds, perturbation_type, 
                           amount_of_SDs_perturbation, nboots,
                           parallel = FALSE, ncores = NULL) {

  if (!requireNamespace("parallel", quietly = TRUE)) {
    install.packages("parallel")
  }
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  if (!requireNamespace("nodeIdentifyR", quietly = TRUE)) {
    devtools::install_github("JasperNaberman/nodeIdentifyR")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  if (parallel) {
    if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
    library(doParallel)
    } else {
      if (!requireNamespace("progress", quietly = TRUE)) {
        install.packages("progress")
      } 
    }

  
  Bootresult <- list(
    BootSamples = vector(mode = "list", length = nboots),
    mean = list(),
    sd = list()
  )
  

  if (parallel) {
    if (is.null(ncores)) {
      ncores <- parallel::detectCores() - 1  
    }
    
    cl <- if (.Platform$OS.type == "unix") {
      makeForkCluster(ncores)
    } else {
      makePSOCKcluster(ncores)
    }
    
    registerDoParallel(cl)
    on.exit(stopCluster(cl)) 
  }
  

  if (!parallel && interactive()) {
    pb <- progress::progress_bar$new(
      format = "[:bar] :percent completed: :eta remaining time",
      total = nboots,  
      clear = FALSE,  
      width = 60    
    )
  }
  

  combine_results <- function(list1, list2) {
    list(
      sample = c(list1$sample,list2$sample),
      means = rbind(list1$means, list2$means),
      sds = rbind(list1$sds, list2$sds)
    )
  }
  
  if (parallel){
    result <- foreach(i = 1:nboots, .combine = "combine_results") %dopar% {
      gs_IsingSamples <- nodeIdentifyR::simulateResponses(edge_weights, thresholds, 
                                                          perturbation_type, amount_of_SDs_perturbation)
      gs_sumIsingSamples <- nodeIdentifyR::calculateSumScores(gs_IsingSamples)
      library(dplyr)
      gs_sumIsingSamplesLong <- nodeIdentifyR::prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
      
      categories <- sort(unique(gs_sumIsingSamplesLong[[2]]))
      category_means <- sapply(categories, function(cat) {
        mean(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
      })
      category_sd <- sapply(categories, function(cat) {
        sd(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
      })
      
      
      if (!parallel && interactive()) {
        pb$tick()
      }
      
      list(
        sample = list(gs_sumIsingSamplesLong),
        means = category_means,
        sds = category_sd
      )
    }
    } else {
      result <- foreach(i = 1:nboots, .combine = "combine_results") %do% {
        gs_IsingSamples <- nodeIdentifyR::simulateResponses(edge_weights, thresholds, 
                                                            perturbation_type, amount_of_SDs_perturbation)
        gs_sumIsingSamples <- nodeIdentifyR::calculateSumScores(gs_IsingSamples)
        library(dplyr)
        gs_sumIsingSamplesLong <- nodeIdentifyR::prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
        
        categories <- sort(unique(gs_sumIsingSamplesLong[[2]]))
        category_means <- sapply(categories, function(cat) {
          mean(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
        })
        category_sd <- sapply(categories, function(cat) {
          sd(gs_sumIsingSamplesLong$sumscore[gs_sumIsingSamplesLong[[2]] == cat], na.rm = TRUE)
        })
        pb$tick()
        list(
          sample = list(gs_sumIsingSamplesLong),
          means = category_means,
          sds = category_sd
        )
      }
    }
  
  Bootresult$BootSamples <- result$sample
  Bootresult$mean <-result$means
  Bootresult$sd <-result$sd
  categories <- sort(unique(result[["sample"]][[1]]$sample))
  colnames(Bootresult$mean) <- categories
  colnames(Bootresult$sd) <- categories
  return(Bootresult)
}





