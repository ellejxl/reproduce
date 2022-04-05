remove(list = ls())

library(Seurat)

#===============load data
Data <- read.csv("MouseALM_HumanMTG.csv",row.names = 1)
Labels3 <- as.matrix(read.csv("MouseALM_HumanMTG_Labels3.csv"))
load("MouseALM_HumanMTG_folds.RData")


#===============normalization:
obj=CreateSeuratObject(counts = t(Data))
obj=NormalizeData(obj, normalization.method = "LogNormalize")

Data_norm=as.matrix(obj@assays$RNA@data)
Data_norm=t(Data_norm)

#===============split:
#label 3
ALM3=data.frame(yy=Labels3,
                  Data_norm)
ALM3=ALM3[Train_Idx[[1]],]
saveRDS(ALM3,file = "ALM3.rds")

MTG3=data.frame(yy=Labels3,
                   Data_norm)
MTG3=MTG3[Train_Idx[[2]],]
saveRDS(MTG3,file = "MTG3.rds")
