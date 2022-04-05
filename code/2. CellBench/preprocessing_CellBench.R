remove(list = ls())

library(Seurat)

#===============load data
Data <- read.csv("Combined_10x_CelSeq2_5cl_data.csv",row.names = 1)
Labels <- as.matrix(read.csv("Labels.csv"))
load("CV_folds.RData")


#===============normalization:
obj=CreateSeuratObject(counts = t(Data))
obj=NormalizeData(obj, normalization.method = "LogNormalize")

Data_norm=as.matrix(obj@assays$RNA@data)
Data_norm=t(Data_norm)

#===============split:
cb_10x=data.frame(yy=Labels,
                  Data_norm)
cb_10x=cb_10x[Train_Idx[[1]],]
saveRDS(cb_10x,file = "cellbench_10x.rds")

cb_celseq2=data.frame(yy=Labels,
                      Data_norm)
cb_celseq2=cb_celseq2[Train_Idx[[2]],]
saveRDS(cb_celseq2,file = "cellbench_celseq2.rds")
