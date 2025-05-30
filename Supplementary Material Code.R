#####置换检验的必要性#####
library(readxl)
single_gds <- read_excel("single_gds.xlsx")#DATA1
PHQ9 <- read_excel("PHQ9.xlsx")#DATA2
GDS9 <- read_excel("GDS9.xlsx")#DATA3
ACE <- read_excel("ACE.xlsx")#DATA4
disease <- read_excel("disease.xlsx")#DATA5


library(bootnet)
sup_fit1<-estimateNetwork(single_gds,default = "IsingFit",rule="AND")
edgeWeightMatrix<-sup_fit1$graph
thresholdVector<-sup_fit1$intercepts
library(nodeIdentifyR)
set.seed(2025)
##alleviating+2SD
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "alleviating", 2)#alleviating
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)
library(dplyr)
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="alleviating")#
library(NIRApost)
stat_result<-getstat(gs_sumIsingSamplesLong,method="bonferroni")
p<-plotmeanNIRA(stat_result,"alleviating")
plot(p)
p$data
##aggravating +2SD
set.seed(2025)
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "aggravating", 2)#alleviating
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)
library(dplyr)
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="aggravating")#
stat_result<-getstat(gs_sumIsingSamplesLong,method="bonferroni")
p<-plotmeanNIRA(stat_result,"aggravating")
plot(p)
p$data


#####稳定性检验的必要性#######
library(bootnet)
sup_fit1<-estimateNetwork(disease,default = "IsingFit",rule="AND")
edgeWeightMatrix<-sup_fit1$graph
thresholdVector<-sup_fit1$intercepts
library(nodeIdentifyR)
set.seed(1111)
##alleviating+2SD
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "alleviating", 2)#alleviating
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)
library(dplyr)
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="alleviating")#
set.seed(2222)
##alleviating+2SD
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "alleviating", 2)#alleviating
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)
library(dplyr)
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="alleviating")#
set.seed(3333)
##alleviating+2SD
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "alleviating", 2)#alleviating
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)
library(dplyr)
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="alleviating")#

library(NIRApost)
#library(dplyr)
set.seed(2025)
#parallel
Bootresult<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"alleviating",2,100,parallel = TRUE, ncores = 6)
#serial
#Bootresult_serial<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"aggravating",2,1000,parallel = FALSE, ncores = NULL)
repeat_node <- findmaxn(Bootresult,n=7)#提取模拟干预后样本均值与模拟原始样本均值差排序前n组中各个节点的百分比和频次
plotmaxn(repeat_node,low=1,high=7)#绘制出各节点排序稳定性图


#####调节效应检验######
library(mgm)
mgm_mod <- mgm(
  data = as.matrix(single_gds),
  type = rep("c", 7),
  level = rep(2, 7),
  lambdaSel = "EBIC", 
  lambdaGam = 0.25,
  ruleReg = "AND",
  moderators = 1,
  scale = FALSE  
)
#mgm_mod
mgm_mod$interactions$indicator[[2]]
modnets::condPlot(mgm_mod,to="TR",from="UW")
modnets::condPlot(mgm_mod,to="UW",from="TR")



graph1<-plot(sup_fit)
library(qgraph)
centralityPlot(graph1,include = c("ExpectedInfluence","Strength","Closeness","Betweenness"),scale = c('z-score'))
edgeWeightMatrix<-sup_fit$graph#提取出边的权重矩阵
thresholdVector<-sup_fit$intercepts#提取出逻辑回归的截取项作为阈值参数
library(nodeIdentifyR)
set.seed(2025)
#加重，aggravating，缓解alleviating
#进行模拟干预，simulateResponses函数进行了"变量数+1"次模拟，
#基于原始ising模型的边权重和阈值参数拟合了一个新的样本，
#然后又分别拟合了“变量数个”满足原始ising模型的边权重
#(加重，aggravating，缓解alleviating)和某个节点2个标准差后的新样本
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "alleviating", 2)#alleviating
#计算每个模拟样本中每个被试的得分总和
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)
library(dplyr)
#将之前生成的原始模拟样本和每次模拟干预样本计算的总分转成长格式，以便绘图
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
#绘制平均总分图
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="alleviating")#

library(NIRApost)
#进行置换检验并提取统计量
stat_result<-getstat(gs_sumIsingSamplesLong,method="bonferroni")
#绘制图
p<-plotmeanNIRA(stat_result,"alleviating")
plot(p)
p$data




library(mgm)
mgm_mod <- mgm(
  data = as.matrix(disease),
  type = rep("c", 8),
  level = rep(2, 8),
  lambdaSel = "EBIC", 
  lambdaGam = 0.25,
  ruleReg = "AND",
  moderators = 2,
  scale = FALSE  
)
#mgm_mod
mgm_mod$interactions$indicator[[2]]
modnets::condPlot(fit1,to="IRR",from="UW")
modnets::condPlot(fit1,to="UW",from="IRR")


######调节效应检验########
source("run_mgmm_analysis.R")
library(mgm)
sig_moder<-run_mgmm_analysis(data=as.matrix(disease), plot_results = FALSE,rule="AND", lambdaGam = 0.25, nB = 10)
print(sig_moder$significant_moderators)

#install.packages(modnets)
fit1 <- modnets::fitNetwork(PHQ9,moderators = "PHQ1",rule="AND")
modnets::condPlot(fit1,to="PHQ3",from="PHQ4")

fit1 <- modnets::fitNetwork(disease,moderators = "Respiratory Disease",rule="AND")
modnets::condPlot(fit1,to=7,from=2)#diabetes

