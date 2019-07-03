BiocManager::install("tximport")
library(tximport)
setwd("./airway2/inst/extdata/quants")
srrs<-list.files("./airway2/inst/extdata/quants")
files<- file.path(srrs, "quant.sf.gz")
file.exists(files)

txi<-tximport(files, type = "salmon", txOut = TRUE)
