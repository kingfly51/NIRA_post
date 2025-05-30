# NIRA_post
The full set of code from the Simulating Interventions for Cross-Sectional Network Models: An R Tutorial article is stored here and is mainly used to perform single network NIRA and dual network NIRA. Among them, the single_gds_CIRA.R file is the complete set of code appearing in the article, and the Supplementary Material Code. R is the code used in the supplementary materials of the article.


##1.What are the reasons for using NIRA? ##

(1)Traditional centrality indices (e.g., strength centrality) only reflect structural information, failing to capture causal directional effects between nodes (e.g., in-influence vs. out-influence). NIRA resolves this by simulating interventions to directly evaluate causal impacts.

(2)By simulating node removal/adjustment (e.g., setting node activation probability to 0), NIRA quantifies each node’s global effect on the network, avoiding reliance on centrality indices alone (e.g., high-centrality nodes may not be optimal targets).

(3)Cross-sectional models lack explicit causal pathways. NIRA infers latent causal relationships through simulated interventions, aiding in identifying core/bridging symptoms.




##2. How to use NIRApost ##


There is a NIRApost.html file in the doc folder, which contains a simple tutorial on how to use NIRApost.




##3.How to install NIRApost ##

Users can install the NIRPostpost package using the following command. Of course, the prerequisite is that the user needs to install the devtools package first.In addition, users can also download and install locally.

install.packages("devtools")

devtools::install_github("kingfly51/NIRA_post")




##4.Introduction to each file in NIRApost##

The data/ directory contains five built-in datasets distributed with the NIRApost package, stored in .rda format for immediate accessibility upon package loading. Users can directly load these datasets into their R environment—for example, the single_gds dataset can be accessed via:

data("single_gds")


The data_raw/ directory stores the original .xlsx files for all five datasets.


The doc/ directory contains three key files: (1)NIRApost.html: A self-contained webpage providing a comprehensive tutorial on using the NIRApost package.(2) NIRApost.Rmd: The source R Markdown document used to generate the tutorial, containing all code, text, and formatting instructions.(3) NIRApost.R: An R script version of the tutorial code (without Markdown text), suitable for direct execution or adaptation.


The man/ directory contains the help documentation files for NIRApost's core functions. Taking the getstat function as an example, users can access its documentation in R using the following command:

help(getstat,package = "NIRApost")


The Meta/ directory contains vignette metadata files, while the vignettes/ folder stores the source R Markdown documents for vignettes. These directories house intermediate files generated during the construction of NIRApost's tutorial materials.


The R/ directory contains all core functions of the NIRApost package.

