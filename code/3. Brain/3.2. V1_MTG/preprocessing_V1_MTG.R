remove(list = ls())

library(Seurat)

#===============load data
Data <- read.csv("MouseV1_HumanMTG.csv",row.names = 1)
Labels3 <- as.matrix(read.csv("MouseV1_HumanMTG_Labels3.csv"))
load("MouseV1_HumanMTG_folds.RData")


#===============normalization:
obj=CreateSeuratObject(counts = t(Data))
obj=NormalizeData(obj, normalization.method = "LogNormalize")

Data_norm=as.matrix(obj@assays$RNA@data)
Data_norm=t(Data_norm)

#===============split:
#label 3
V13=data.frame(yy=Labels3,
                  Data_norm)
V13=V13[Train_Idx[[1]],]
saveRDS(V13,file = "V13.rds")

MTG3=data.frame(yy=Labels3,
                   Data_norm)
MTG3=MTG3[Train_Idx[[2]],]
saveRDS(MTG3,file = "MTG3.rds")

