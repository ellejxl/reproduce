remove(list = ls())

#install packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("CHETAH")


library("SingleCellExperiment")
library("CHETAH")
set.seed(0)
args=commandArgs(trailingOnly = TRUE)
train_inp=as.character(args[1])
test_inp=as.character(args[2])

print(paste0("============ CHETAH: train=",train_inp," test=",test_inp,"======"))


# !!! 1. Load Data
#=================================================
#read datasets
train_dat=readRDS(paste0(train_inp,".rds"))
test_dat=readRDS(paste0(test_inp,".rds"))

print(paste0("initial genes number: ",ncol(train_dat)))

#!!! 2. pre-processing data
#=================================================
#!!!a) remove the celltype which is not present in the training
train_cellnames=names(table(train_dat[,1]))
test_dat=test_dat[which(test_dat[,1]%in%train_cellnames),]


#training data
train_label=train_dat[,1]
sce=SingleCellExperiment(assay=list(counts=as.matrix(t(train_dat[,-1]))),
                         colData = DataFrame(train_label))
# use gene names as feature symbols
colData(sce)$celltypes <- as.character(colData(sce)$train_label)

#testing data
test_label=test_dat[,1]
sce_test=SingleCellExperiment(assay=list(counts=as.matrix(t(test_dat[,-1]))),
                              colData = DataFrame(test_label))
colData(sce_test)$celltypes <- as.character(colData(sce_test)$test_label)

sce_test=CHETAHclassifier(input = sce_test,ref_cells = sce)


labels=data.frame(cell_label=colData(sce_test)$celltypes,
                  prediction=sce_test$celltype_CHETAH)


saveRDS(labels,paste0("labels_CHETAH_",train_inp,"_",test_inp,".rds"))


