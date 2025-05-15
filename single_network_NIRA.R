####### 1.导入数据集并进行清理#######
rm(list=ls())#清空工作空间
library(readxl)
single_gds <- read_excel("single_gds.xlsx")#导入数据集
GS_non_null<- na.omit(single_gds[,c(1:7)])#删除所有有缺失值的被试


####### 2.构建Ising模型#######
library(bootnet)
gs_fit<-estimateNetwork(GS_non_null,default = "IsingFit",rule="AND")#拟合Ising模型
type=rep(c("GAD"),c(7))
#plot.new()
#spring,circle
graph1<-plot(gs_fit, layout = 'spring',groups=type,color=c("#F9EB77"),legend=FALSE)#绘图
#legend("topright", legend = c("GAD"), fill = c("#F9EB77"))


###### 3.计算中心性指标 ######
library(qgraph)
#画中介性，紧密度，强度
centralityPlot(graph1,include = c("ExpectedInfluence","Strength","Closeness","Betweenness"),scale = c('z-score'))
#bootstrap检验节点稳定性
Cortral<-bootnet(gs_fit,statistics = c("ExpectedInfluence","strength","closeness","betweenness"),nBoots = 1000,nCores=8,type = "case",missing="pairwise")
corStability(Cortral)
plot(Cortral,statistics = c("expectedInfluence","strength","closeness","betweenness"))



###### 4.nodeIdentifyR ######
edgeWeightMatrix<-gs_fit$graph#提取出边的权重矩阵
thresholdVector<-gs_fit$intercepts#提取出逻辑回归的截取项作为阈值参数
#导入nodeIdentifyR包并载入到当前工作空间
library(devtools)
devtools::install_github("JasperNaberman/nodeIdentifyR")
library(nodeIdentifyR)
set.seed(2025)
#加重，aggravating，缓解alleviating
#进行模拟干预，simulateResponses函数进行了"变量数+1"次模拟，
#基于原始ising模型的边权重和阈值参数拟合了一个新的样本，
#然后又分别拟合了“变量数个”满足原始ising模型的边权重
#(加重，aggravating，缓解alleviating)和某个节点2个标准差后的新样本
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "aggravating", 2)#alleviating
#计算每个模拟样本中每个被试的得分总和
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)
library(dplyr)
#将之前生成的原始模拟样本和每次模拟干预样本计算的总分转成长格式，以便绘图
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
#绘制平均总分图
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="aggravating")#


####显著性检验：置换检验
source('getstat.R')
source('permutation_test.R')
source('plotmeanNIRA.R')
#进行置换检验并提取统计量
stat_result<-getstat(gs_sumIsingSamplesLong,method="bonferroni")
#绘制图
plotmeanNIRA(stat_result,"aggravating")


####稳定性检验
source('BootNIRAresult.R')
source('findmaxn.R')
source('plotmaxn.R')
#library(dplyr)
set.seed(2025)
#parallel
Bootresult<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"aggravating",2,100,parallel = TRUE, ncores = 6)
#serial
#Bootresult_serial<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"aggravating",2,1000,parallel = FALSE, ncores = NULL)
repeat_node <- findmaxn(Bootresult,n=7)#提取模拟干预后样本均值与模拟原始样本均值差排序前n组中各个节点的百分比和频次
plotmaxn(repeat_node,low=1,high=3)#绘制出各节点排序稳定性图

##alleviating
Bootresult_allev<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"alleviating",2,10,parallel = TRUE, ncores = 6)
#serial
#Bootresult_serial_allev<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"alleviating",2,1000,parallel = FALSE, ncores = NULL)
repeat_node_allev <- findmaxn(Bootresult_allev,n=7)#提取模拟干预后样本均值与模拟原始样本均值差排序前n组中各个节点的百分比和频次
plotmaxn(repeat_node_allev,low=1,high=7)#绘制出各节点排序稳定性图
stat_result_allev<-getstat(Bootresult_allev)
plotmeanNIRA(stat_result_allev,"alleviating")



#####有调节的网络分析
source("run_mgmm_analysis.R")
library(mgm)
sig_moder<-run_mgmm_analysis(data=as.matrix(GS_non_null),node_num=7, plot_results = FALSE, lambdaGam = 0.25, nB = 10)
print(sig_moder$significant_moderators)

#install.packages(modnets)
fit1 <- modnets::fitNetwork(GS_non_null,moderators = "ASH",rule="AND")
modnets::condPlot(fit1,to="IRR",from="UW")
modnets::condPlot(fit1,to="UW",from="IRR")




