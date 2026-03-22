## data-raw/prepare_single_gds.R
## 用于从原始数据生成 data/single_gds.rda
## 运行方法: usethis::use_data_raw("single_gds") 后替换本文件内容，
##           然后执行 source("data-raw/prepare_single_gds.R")

# ---------- 方法 1: 从 Excel 文件读取 ----------
# library(readxl)
# single_gds <- readxl::read_excel("data-raw/single_gds.xlsx")
# single_gds <- as.data.frame(single_gds)

# ---------- 方法 2: 从 CSV 文件读取 ----------
# single_gds <- read.csv("data-raw/single_gds.csv")

# ---------- 保存为包内数据 ----------
# usethis::use_data(single_gds, overwrite = TRUE, compress = "xz")

