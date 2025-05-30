#####1.1 Install the required R packages#####
install.packages("readxl") # install the readxl package
install.packages("bootnet") # install the bootnet package
install.packages("qgraph")# install the qgraph package
install.packages("devtools")# install the devtools package
devtools::install_github("JasperNaberman/nodeIdentifyR")# install the nodeIdentifyR package
install.packages("dplyr")# install the dplyr package
install.packages("mgm")# install the mgm package
devtools::install_github("kingfly51/NIRA_post")# install the NIRApost package



#####1.2 Load the required R packages into the current workspace#####
library("readxl") # Load the readxl package
library("bootnet")# Load the bootnet package
library("qgraph")# Load the qgraph package
library("dplyr")# Load the dplyr package
library("nodeIdentifyR")# Load the nodeIdentifyR package
library("NIRApost")# Load the NIRApost package



#####1.3 Import dataset #####
#single_gds <- read_excel("single_gds.xlsx") #Import dataset,Method 1
#single_gds <- read_excel("D:/Rdaima/AMPPS_NIRA/NIRApost/data_raw/single_gds.xlsx")#,Method 2
data("single_gds",package="NIRApost")#Method 3
GS_non_null<- na.omit(single_gds[,c(1:7)])#Delete all subjects with missing values



#####1.4 Construct an Ising network and perform network analysis######
gs_fit<-estimateNetwork(GS_non_null,default = "IsingFit",rule="AND")#Fit Ising model
type=rep(c("GAD"),c(7))#Set variable category
#Draw the network diagram of the Ising model
graph1<-plot(gs_fit, layout ="spring",groups=type,color=c("#F9EB77"),legend=FALSE)
centralityPlot(graph1,include = c("ExpectedInfluence","Strength","Closeness","Betweenness"),scale = c('z-score'))#Plot the four centrality measures for all nodes in the network
set.seed(2025)
Cortral<-bootnet(gs_fit,statistics = c("ExpectedInfluence","strength","closeness","betweenness"),nBoots = 1000,nCores=8,type = "case",missing="pairwise")#Bootstrap checks node stability
corStability(Cortral)#View the CS coefficients of each indicator
#Plot the case-drop bootstrap plot for the stability test
plot(Cortral,statistics = c("expectedInfluence","strength","closeness","betweenness"))



#####1.5 Simulation intervention based on NodeIdentififyR algorithm ######
#Extract the edge weight matrix component from the previously fitted Ising model gs_fit, which is stored in the graph object
edgeWeightMatrix<-gs_fit$graph
#Extract the intercept of the logistic regression as the threshold parameter vector
thresholdVector<-gs_fit$intercepts
set.seed(2025)#Set seed
#Generate simulated origin samples and simulated post- aggravating intervention samples
gs_IsingSamples<-simulateResponses(edgeWeightMatrix, thresholdVector, "aggravating", 2)
#Calculate the total score of each subject in each simulated sample
gs_sumIsingSamples<-calculateSumScores(gs_IsingSamples)
#Wide format to long format
gs_sumIsingSamplesLong<-prepareDFforPlottingAndANOVA(gs_sumIsingSamples)
#Draw the average total score chart
plotSumScores(sum_scores_long=gs_sumIsingSamplesLong,perturbation_type="aggravating")



#####1.6 Permutation test for NIRA ######
#Perform Permutation Testing and Extract Test Statistics
stat_result<-getstat(gs_sumIsingSamplesLong,method="bonferroni")
#Print the statistical characteristics of the post-intervention sample means for each node.
print(stat_result$stat)
#Visualize Mean Values with Confidence Intervals from NIRA Analysis
plotmeanNIRA(stat_result,"aggravating")



#####1.7 Stability testing of NIRA ######
set.seed(2025)#set seed
#parallel
Bootresult<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"aggravating",2,100,parallel = TRUE, ncores = 6)
#serial
#Bootresult_serial<-BootNIRAresult(edgeWeightMatrix,thresholdVector,"aggravating",2,100,parallel = FALSE, ncores = NULL)
# Extracting Node Percentages and Frequencies in Top-7 Differential Groups
repeat_node <- findmaxn(Bootresult,n=7)
#View the stability indicators of each node in the stability test
View(repeat_node)
#Visualization of Node Ranking Stability
plotmaxn(repeat_node,low=1,high=7)



#####1.8 NIRA's test of moderation effects ######
#Calculate all stable moderation effects present in the network
sig_moder<-run_mgmm_analysis(data=as.matrix(GS_non_null), plot_results = FALSE, rule = "AND", lambdaGam = 0.25, nB = 100)
#Print out the moderation effects
print(sig_moder$significant_moderators)



#####1.9 Visualization of Moderation Effects (Optional) #####
#Install 'modnets' if not already installed (required only for first-time use of this R package).
install.packages(modnets)
#Constructing a moderated network using ASH as a moderating variable
fit1 <- modnets::fitNetwork(GS_non_null,moderators = "ASH",rule="AND")
#Plot the moderation effect diagram of UW on IRR
modnets::condPlot(fit1,to="IRR",from="UW")
#Plot the moderation effect diagram of IRR on UW
modnets::condPlot(fit1,to="UW",from="IRR")

