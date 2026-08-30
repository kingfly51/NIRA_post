#Note! The updated version NIRApost 1.2.0 introduces a workflow based on analytical solutions. In NIRApost 1.1.0, 
#Monte Carlo simulations were used for approximation, thus requiring stability checks. However, Cui et al. (2026)
#pointed out that there is no need to employ simulation methods that incorporate Monte Carlo error, because an analytical method 
#exists that provides an error-free solution. Moreover, this approach is computationally very fast. Therefore, the updated NIRApost 1.2.0 
#includes a new function, analyticalNIRAtest, for the analytical solution-based NIRA procedure. This function encompasses the computation
#of analytical solutions, uncertainty interval estimation (providing both bootstrap tests and permutation tests), and plotting capabilities.
#Consequently, NIRApost now offers two distinct NIRA workflows:
#The first is the four-step NIRA approach combined with moderation effect testing, permutation testing, and stability testing, as demonstrated in the single_gds_NIRA.R file.
#The second is the moderation effect testing combined with analytical solution-based NIRA, as demonstrated in the single_gad_analyticalNIRA.R file.


#####1.1 Install the required R packages#####
# Install and load required packages
if (!require("bootnet")) install.packages("bootnet")       # install the bootnet package
if (!require("qgraph")) install.packages("qgraph")         # install the qgraph package
if (!require("devtools")) install.packages("devtools")     # install the devtools package
if (!require("dplyr")) install.packages("dplyr")           # install the dplyr package
if (!require("mgm")) install.packages("mgm")               # install the mgm package

# Install GitHub packages if not already installed
if (!require("nodeIdentifyR")) {
  devtools::install_github("JasperNaberman/nodeIdentifyR") # install the nodeIdentifyR package
}
if (!require("NIRApost")) {
  devtools::install_github("kingfly51/NIRA_post")          # install the NIRApost package
}



#####1.2 Load the required R packages into the current workspace#####
library("bootnet")       # Load the bootnet package
library("qgraph")        # Load the qgraph package
library("dplyr")         # Load the dplyr package
library("nodeIdentifyR") # Load the nodeIdentifyR package
library("NIRApost")      # Load the NIRApost package



#####1.3 Import dataset#####
#single_gds <- read_excel("single_gds.xlsx")                                           # Import dataset, Method 1
#single_gds <- read_excel("D:/Rdaima/AMPPS_NIRA/NIRApost/data_raw/single_gds.xlsx")   # Import dataset, Method 2

# Load the built-in dataset from the NIRApost package
data("single_gds", package = "NIRApost")
# Explicitly retrieve the dataset to avoid unevaluated promises
data <- get("single_gds")      # Method 3
GS_non_na <- na.omit(data)     # Delete all subjects with missing values



#####1.4 Construct an Ising network and perform network analysis######
fit <- estimateNetwork(GS_non_na, default = "IsingFit", rule = "AND") # Fit Ising model
type <- rep(c("GAD"), c(7))    # Set variable category
# Draw the network diagram of the Ising model
graph1 <- plot(fit, layout = "spring", groups = type, color = c("#F9EB77"), legend = FALSE)
centralityPlot(Objects = graph1,
               include = c("ExpectedInfluence", "Strength"),
               scale = c("z-score"))   # Plot the two centrality measures for all nodes in the network
set.seed(2025)
Cortral <- bootnet(data = fit,
                   statistics = c("ExpectedInfluence", "strength"),
                   nBoots = 1000,
                   nCores = 8,
                   type = "case",
                   missing = "pairwise") # Bootstrap checks node stability
corStability(Cortral)                    # View the CS coefficients of each indicator
# Plot the case-drop bootstrap plot for the stability test
plot(Cortral, statistics = c("expectedInfluence", "strength"))



#####1.5 NIRA's test of moderation effects######
# Calculate all stable moderation effects present in the network
set.seed(2025)
sig_moder <- runMgmmAnalysis(data = as.matrix(GS_non_na),
                             plotResults = FALSE,
                             rule = "AND",
                             lambdaGam = 0.25,
                             nB = 10)
# Print out the moderation effects
print(sig_moder$significant_moderators)



#####1.6 Visualization of Moderation Effects (Optional)#####
# Install 'modnets' if not already installed (required only for first-time use of this R package).
if (!require("modnets")) install.packages("modnets")
library(modnets)
# Constructing a moderated network using ASH as a moderating variable
fit1 <- fitNetwork(GS_non_na, moderators = "ASH", rule = "AND")
# Plot the moderation effect diagram of UW on IRR
condPlot(fit1, to = "IRR", from = "UW")
# Plot the moderation effect diagram of IRR on UW
condPlot(fit1, to = "UW", from = "IRR")



#####1.7 Simulation intervention based on analytical NodeIdentifyR algorithm######
source('D:/Rdaima/AMPPS_NIRA/fourth_revise/Cui_comments/NIRA_post-main/R/analyticalNIRAtest.r')
# Extract the edge weight matrix component from the previously fitted Ising model
edgeWeightMatrix <- fit$graph
# Extract the intercept of the logistic regression as the threshold parameter vector
thresholdVector <- fit$intercepts

res <- analyticalNIRAtest(
   data              = single_gds,
   perturbation_type = "aggravating",
   edge_weights      = edgeWeightMatrix,
   thresholds        = thresholdVector,
   amount_of_SDs_perturbation = 2,
   nResample = 5000,
   method = "bonferroni",
   seed = 2025,
   parallel = TRUE,
   ncores = 6
)

res$stat
res$random
res$plot

res$exact










