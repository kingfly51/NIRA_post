####### 1.导入数据集并进行清理#######
rm(list=ls())#清空工作空间
library(readxl)
bridge_network <- read_excel("bridge_network.xlsx")#导入数据集
GS_non_null<- na.omit(bridge_network[,c(1:11)])#删除所有有缺失值的被试




####### 2.构建Ising模型#######
library(bootnet)
gs_fit<-estimateNetwork(GS_non_null,default = "IsingFit")#拟合Ising模型
type=rep(c("JKXW","MH"),c(8,3))
plot.new()
plot(gs_fit, layout = 'spring',groups=type,color=c("#04CED1","#F9EB77"),legend=FALSE)#绘图
legend("topright", legend = c("JKXW","MH"), fill = c("#04CED1", "#F9EB77"))


#计算桥接指标
library(networktools)
grapha1<-gs_fit$graph
b<-bridge(grapha1,communities= type)
c<-expectedInf(grapha1)
summary(b)
summary(c)
plot(b,order="value", zscore=T,include=c("Bridge Expected Influence (1-step)","Bridge Strength", "Bridge Betweenness","Bridge Closeness"),color = TRUE)
plot(c)
#建议桥接中心性指数的稳定性
Results1<-bootnet(gs_fit,statistics = c("bridgeExpectedInfluence","bridgeStrength","bridgeCloseness","bridgeBetweenness"),type="person",communities= type,nBoots = 100,nCores=7)#,caseN=80，nonparametric
corStability(Results1)
plot(Results1,statistics = c("bridgeExpectedInfluence","bridgeStrength","bridgeCloseness","bridgeBetweenness"))



edgeWeightMatrix<-gs_fit$graph#提取出边的权重矩阵
thresholdVector<-gs_fit$intercepts#提取出逻辑回归的截取项作为阈值参数




####### nodeIdentifyR#######
library(NIRApost)

# Calculate all stable moderation effects present in the network
set.seed(2025)
sig_moder <- runMgmmAnalysis(data = as.matrix(GS_non_null),
                             plotResults = FALSE,
                             rule = "AND",
                             lambdaGam = 0.25,
                             nB = 10)
# Print out the moderation effects
print(sig_moder$significant_moderators)

groups <- rep(c("JKXW","MH"),c(8,3))

res2 <- analyticalBridgeNIRAtest(data=GS_non_null,
                        groups=groups,
                        intervention_group="JKXW",
                        outcome_group="MH",
                        perturbation_type = "aggravating",
                        edge_weights = edgeWeightMatrix,
                        thresholds =thresholdVector,
                        amount_of_SDs_perturbation = 2,
                        nResample = 1000,
                        method = "bonferroni",
                        seed = 2025,
                        parallel = TRUE,
                        ncores = 6) 

res2$stat            # per-source-node inference on the outcome community

res$random          # target-vs-random-source p-values (uncorrected)

res2$exact           # exact outcome-community mean activation

res2$plot
