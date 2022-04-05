remove(list = ls())

#install packages
#install.packages("BiocManager")
#BiocManager::install(c("S4Vectors", "hopach", "limma"))
#install.packages("devtools")
#devtools::install_github("SydneyBioX/scClassify")

library("scClassify")
set.seed(0)

args=commandArgs(trailingOnly = TRUE)
train_inp=as.character(args[1])
test_inp=as.character(args[2])


print(paste0("============ scClassify: train=",train_inp," test=",test_inp,"======"))


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
train_yy=train_dat[,1]

exprsMat_test=t(test_dat[,-1])
test_yy=test_dat[,1]


scClassify_res=scClassify(exprsMat_train = exprsMat_train,
                          cellTypes_train = train_yy,
                          exprsMat_test = list(test=exprsMat_test),
                          cellTypes_test=list(test=test_yy),
                          tree = "HOPACH",
                          algorithm = "WKNN",
                          selectFeatures = c("limma"),
                          similarity = c("pearson"),
                          returnList = FALSE,
                          verbose = FALSE)


pred.label<-scClassify_res$testRes$test$pearson_WKNN_limma$predRes

labels=data.frame(cell_label=test_yy,
                  prediction=pred.label)


saveRDS(labels,paste0("labels_scClassify_",train_inp,"_",test_inp,".rds"))



