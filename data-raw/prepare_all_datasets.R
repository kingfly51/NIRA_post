## data-raw/prepare_all_datasets.R
## 从 xlsx 读取原始数据，二值化，保存为 data/*.rda
## 执行: source("data-raw/prepare_all_datasets.R")
library(readxl)
DATA_PATH <- "D:/Rdaima/AMPPS_NIRA/fourth_revise/NIRApost/data-raw"
PKG_PATH  <- "D:/Rdaima/AMPPS_NIRA/fourth_revise/NIRApost"
save_rda <- function(obj, name) {
  assign(name, obj)
  save(list = name,
       file = file.path(PKG_PATH, "data", paste0(name, ".rda")),
       compress = "xz")
  message("Saved: data/", name, ".rda  [", nrow(obj), " x ", ncol(obj), "]")
}
# 1. single_gds — GAD-7 (N=2404): 0="从不", 其余=1
raw <- as.data.frame(readxl::read_excel(file.path(DATA_PATH,"single_gds.xlsx")))
single_gds <- as.data.frame(lapply(raw, function(x) ifelse(is.na(x),NA,ifelse(x==0,0L,1L))))
save_rda(single_gds, "single_gds")
# 2. PHQ9 — PHQ-9 (N=45829): 0="完全不会", 其余=1
raw <- as.data.frame(readxl::read_excel(file.path(DATA_PATH,"PHQ9.xlsx")))
PHQ9 <- as.data.frame(lapply(raw, function(x) ifelse(is.na(x),NA,ifelse(x==0,0L,1L))))
save_rda(PHQ9, "PHQ9")
# 3. GDS9 — 已二值 (N=3097)
raw <- as.data.frame(readxl::read_excel(file.path(DATA_PATH,"GDS9.xlsx")))
GDS9 <- as.data.frame(lapply(raw, function(x) ifelse(is.na(x),NA,as.integer(x))))
save_rda(GDS9, "GDS9")
# 4. ACE — 已二值 (N=525)
raw <- as.data.frame(readxl::read_excel(file.path(DATA_PATH,"ACE.xlsx")))
ACE <- as.data.frame(lapply(raw, function(x) ifelse(is.na(x),NA,as.integer(x))))
save_rda(ACE, "ACE")
# 5. disease — 已二值 (N=27353)
raw <- as.data.frame(readxl::read_excel(file.path(DATA_PATH,"disease.xlsx")))
disease <- as.data.frame(lapply(raw, function(x) ifelse(is.na(x),NA,as.integer(x))))
save_rda(disease, "disease")
message("\n全部数据集保存至: ", file.path(PKG_PATH, "data"))

