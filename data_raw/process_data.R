# library package
library(readxl)
library(usethis)

# Read and process data
single_gds <- as.data.frame(readxl::read_excel("data_raw/single_gds.xlsx"))
ACE <- as.data.frame(readxl::read_excel("data_raw/ACE.xlsx"))
disease <- as.data.frame(readxl::read_excel("data_raw/disease.xlsx"))
GDS9 <- as.data.frame(readxl::read_excel("data_raw/GDS9.xlsx"))
PHQ9 <- as.data.frame(readxl::read_excel("data_raw/PHQ9.xlsx"))

# Save the processed data in R format
usethis::use_data(single_gds, ACE, disease, GDS9, PHQ9, overwrite = TRUE)
