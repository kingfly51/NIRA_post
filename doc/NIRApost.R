## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----package_load-------------------------------------------------------------
library(NIRApost)#Load NIRApost Package
library(bootnet)#Load bootnet Package
library(nodeIdentifyR)#load nodeIdentifyR Package 
library(dplyr)#Load dplyr Package

## ----Ising_Network_Fit, echo=TRUE---------------------------------------------
#data(package="NIRApost")#View the dataset stored in the NIRAspost package
data("single_gds")#load data
gs_fit<-estimateNetwork(single_gds,default = "IsingFit",rule="AND")#Fit Ising model
edgeWeightMatrix<-gs_fit$graph#Extract the weight matrix of edges
thresholdVector<-gs_fit$intercepts#Extract the intercept term from logistic regression as the threshold parameter

## ----Simulated_alleviate_intervention-----------------------------------------
set.seed(2025)#set seed
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "alleviating", 2)#Generate N+1 simulated samples, where N represents simulated alleviate intervention samples for N nodes and 1 represents simulated origin samples
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)#Calculate the total score of each participant in each simulated sample
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)#Convert the total score calculated from the previously generated origin simulation samples and the simulation intervention samples of each node into a long format for plotting purposes
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="alleviating")#Draw an average total score figure


## ----permutation_test---------------------------------------------------------
stat_result<-getstat(gs_sumIsingSamplesLong,method="bonferroni")#Perform permutation test and extract statistical measures, and use Bonferroni for multiple comparison correction
head(stat_result$stat)#View statistics
plotmeanNIRA(stat_result,"alleviating")#Visualize Mean Values with Confidence Intervals from NIRA Analysis,Please note,The 'alleviating' needs to be consistent with the perturbation_type mentioned earlier

## ----stability_test-----------------------------------------------------------
#parallel code
Bootresult<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"alleviating",2,100,parallel = TRUE, ncores = 6)#Perform parallel and repeated simulation interventions with 6 cores for 100 times. For demonstration purposes, only 100 times will be repeated here

#serial code
#Bootresult<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"alleviating",2,100,parallel = FALSE, ncores = NULL)#It is worth noting that serial has a timeline, but parallel does not

repeat_node <- findmaxn(Bootresult,n=7)#Extract the percentage and frequency of each node appearing in the top n absolute mean difference rankings between post-intervention simulated samples and the original simulated sample.
plotmaxn(repeat_node,low=1,high=7)#Generate stability plots of node ranking orders.


## ----Moderation_Effect_Test---------------------------------------------------
sig_moder<-run_mgmm_analysis(data=as.matrix(single_gds), plot_results = FALSE,rule = "AND", lambdaGam = 0.25, nB = 10)#Iteratively calculate the moderation effect of each node in the network.The purpose of the demonstration is to set nB to 10
print(sig_moder$significant_moderators)#Print the results of significant moderation effects

