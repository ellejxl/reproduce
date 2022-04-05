remove(list = ls())

#install packages
#install.packages("devtools")
#devtools::install_github("BatadaLab/scID")
library("scID")
set.seed(0)

args=commandArgs(trailingOnly = TRUE)
train_inp=as.character(args[1])
test_inp=as.character(args[2])

print(paste0("============ scID: train=",train_inp," test=",test_inp,"======"))



# !!! 1. Load Data
#=================================================
#read datasets
train_dat=readRDS(paste0(train_inp,".rds"))
test_dat=readRDS(paste0(test_inp,".rds"))

#transfer back 
train_dat[,-1]=exp(train_dat[,-1])-1
test_dat[,-1]=exp(test_dat[,-1])-1


print(paste0("initial genes number: ",ncol(train_dat)))

#!!! 2. pre-processing data
#=================================================
#!!!a) remove the celltype which is not present in the training
train_cellnames=names(table(train_dat[,1]))
test_dat=test_dat[which(test_dat[,1]%in%train_cellnames),]


#target_gem:Data frame of gene expression (rows) per cell (columns) in target data
target_gem=as.data.frame(t(test_dat[,-1]))
#reference_gem: Data frame of gene expression (rows) per cell (columns) in reference data
reference_gem=as.data.frame(t(train_dat[,-1]))
#reference_clusters: Named list of cluster IDs of reference cells
reference_clusters=as.factor(train_dat[,1])
names(reference_clusters)=colnames(reference_gem)

scid_pred=scid_multiclass(target_gem = target_gem,
                          reference_gem = reference_gem,
                          reference_clusters = reference_clusters,
                          normalize_reference = FALSE)

labels=data.frame(cell_label=test_dat[,1],
                  prediction=scid_pred$labels)


saveRDS(labels,paste0("labels_scID_",train_inp,"_",test_inp,".rds"))


