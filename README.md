# NIRA_post
The full set of code from the Simulating Interventions for Cross-Sectional Network Models: An R Tutorial article is stored here and is mainly used to perform single network NIRA and dual network NIRA.

##1.How to install NIRApost ##

Users can install the NIRPostpost package using the following command. Of course, the prerequisite is that the user needs to install the devtools package first.In addition, users can also download and install locally.

install.packages("devtools")

devtools::install_github("kingfly51/NIRA_post")


##2.Introduction to each file in NIRApost##

The data/ directory contains five built-in datasets distributed with the NIRApost package, stored in .rda format for immediate accessibility upon package loading. Users can directly load these datasets into their R environment—for example, the single_gds dataset can be accessed via:

data("single_gds")

The data_raw/ directory stores the original .xlsx files for all five datasets.



