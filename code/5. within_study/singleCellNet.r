remove(list = ls())

#install packages
#install.packages("devtools")
#devtools::install_github("pcahan1/singleCellNet")
library("singleCellNet")
set.seed(0)

args=commandArgs(trailingOnly = TRUE)
train_inp=as.character(args[1])
test_inp=as.character(args[2])

print(paste0("============ singleCellNet: train=",train_inp," test=",test_inp,"======"))

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

exprsMat_train=t(train_dat[,-1])

exprsMat_test=t(test_dat[,-1])

train_yy=data.frame(label=train_dat[,1],sample=colnames(exprsMat_train))
rownames(train_yy)=colnames(exprsMat_train)
test_yy=data.frame(label=test_dat[,1],sample=colnames(exprsMat_test))
rownames(test_yy)=colnames(exprsMat_test)

#train the classifier with training data
class_info<-scn_train(stTrain = train_yy,
                      expTrain = exprsMat_train, 
                      nTopGenes = 10,
                      nRand = 70, 
                      nTrees = 1000, 
                      nTopGenePairs = 25,
                      dLevel = "label", 
                      colName_samp = "sample")

#make predition with testing data
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']],
                               expDat=exprsMat_test,
                               nrand = 50)
test_yy=get_cate(classRes = classRes_val_all, 
                 sampTab = test_yy, 
                 dLevel = "label",
                 sid = "sample",
                 nrand = 50)



labels=data.frame(cell_label=test_yy$label,
                  prediction=test_yy$category)

saveRDS(labels,paste0("labels_singleCellNet_",train_inp,"_",test_inp,".rds"))

