remove(list = ls())
library(SeuratData)
library(Seurat)
#InstallData("pbmcsca")
#loading data
data("pbmcsca")

#===============split
#remove the Unassigned cells
pbmcsca=subset(pbmcsca, subset = CellType!="Unassigned")

#combine 10x Chromium (v2) A and 10x Chromium (v2) B as 10x Chromium (v2)
table(pbmcsca$Method)
pbmcsca$Method[which(pbmcsca$Method=="10x Chromium (v2) A")]="10x-v2"
pbmcsca$Method[which(pbmcsca$Method=="10x Chromium (v2) B")]="10x-v2"
pbmcsca$Method[which(pbmcsca$Method=="10x Chromium (v2)")]="10x-v2"
pbmcsca$Method[which(pbmcsca$Method=="10x Chromium (v3)")]="10x-v3"

table(pbmcsca$Method)


#normalization:
pbmcsca=NormalizeData(pbmcsca, normalization.method = "LogNormalize")

pbmc_list=SplitObject(pbmcsca, split.by = "Method")

pbmc_methods=names(pbmc_list)

remove(pbmcsca)


#split by methods
for(ii in 1:length(pbmc_methods)){
  #===============splite
  data1=as.matrix(pbmc_list[[ii]]@assays$RNA@data)
  celltype=pbmc_list[[ii]]@meta.data$CellType
  
  ref.df=data.frame(yy=celltype,t(data1))
  saveRDS(ref.df,paste0("pbmc_",pbmc_methods[ii],".rds"))
}



