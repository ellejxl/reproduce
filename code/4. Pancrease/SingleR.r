remove(list = ls())

#install packages

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("SingleR")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("scran")

library("SingleR")
library("scran")
set.seed(0)

args=commandArgs(trailingOnly = TRUE)
train_inp=as.character(args[1])
test_inp=as.character(args[2])

print(paste0("============ SingleR: train=",train_inp," test=",test_inp,"======"))

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

train_yy=train_dat[,1]
exprsMat_train=t(train_dat[,-1])

test_yy=test_dat[,1]
exprsMat_test=t(test_dat[,-1])


pred.grun=SingleR(test = exprsMat_test,
                  ref = exprsMat_train,
                  labels = train_yy,
                  de.method = "wilcox")


labels=data.frame(cell_label=test_yy,
                  prediction=pred.grun$labels)


saveRDS(labels,paste0("labels_SingleR_",train_inp,"_",test_inp,".rds"))


