remove(list=ls())

#+++++++++++++++++++++++++++++++++++++++++
args=commandArgs(trailingOnly = TRUE)
train_inp=as.character(args[1])
test_inp=as.character(args[2])

print(paste0("============mixture model: train=",train_inp," test=",test_inp,"======"))

library(scAnnotate)

set.seed(0)

# !!! 1. Load Data
#=================================================
#read datasets
train_dat=readRDS(paste0(train_inp,".rds"))
test_dat=readRDS(paste0(test_inp,".rds"))

print(paste0("initial genes number: ",ncol(train_dat)))

#arrange the dataset with cell type
train_dat=train_dat[order(train_dat[,1]),]
train_cellnames=names(table(train_dat[,1]))

#remove the celltype which is not present in the training
test_dat=test_dat[which(test_dat[,1]%in%train_cellnames),]

test_yy=test_dat[,1]
test_cellnames=names(table(test_yy))

test_xx=test_dat[,-1]
remove(test_dat)

predict=scAnnotate(train=train_dat,
                   test=test_xx,
                   distribution="normal",
                   correction="auto",
                   threshold=0)

labels=data.frame(cell_label=test_yy,
                  prediction=predict)

saveRDS(labels,paste0("labels_scAnnotate_",train_inp,"_",test_inp,".rds"))

