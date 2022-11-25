remove(list = ls())

#install packages
#install.packages("devtools")
#devtools::install_github("immunogenomics/harmony")
#devtools::install_github("powellgenomicslab/scPred")

library("Seurat")
library("scPred")
library("magrittr")
set.seed(0)

args=commandArgs(trailingOnly = TRUE)
train_inp=as.character(args[1])
test_inp=as.character(args[2])

print(paste0("============ scPred: train=",train_inp," test=",test_inp,"======"))

#=====load data
train_dat=readRDS(paste0(train_inp,".rds"))
test_dat=readRDS(paste0(test_inp,".rds"))

#!!!preprocessing data
#=================================================
#remove the celltype which is not present in the training
train_cellnames=names(table(train_dat[,1]))
test_dat=test_dat[which(test_dat[,1]%in%train_cellnames),]

train_xx=train_dat[,-1]
test_xx=test_dat[,-1]

#===============transfer matrix to seurat
#train
#col data include cell type 
cData <- data.frame(cell_type= train_dat[,1])
rownames(cData) <- colnames(train_dat[,1])
#transfer matrix to SingleCellAssay
reference <- CreateSeuratObject(counts=as.matrix(t(train_xx)))
reference[["cell_type"]]=cData$cell_type

#test
#col data include cell type 
cData <- data.frame(cell_type= test_dat[,1])
rownames(cData) <- colnames(test_dat[,1])
#transfer matrix to SingleCellAssay
query <- CreateSeuratObject(counts=as.matrix(t(test_xx)))
query[["cell_type"]]=cData$cell_type


reference <- reference %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)


reference <- getFeatureSpace(reference, "cell_type")

reference <- trainModel(reference)


query <- scPredict(query, reference)


labels=data.frame(cell_label=query@meta.data$cell_type,
                  prediction=query@meta.data$scpred_prediction)

saveRDS(labels,paste0("labels_scPred_",train_inp,"_",test_inp,".rds"))

