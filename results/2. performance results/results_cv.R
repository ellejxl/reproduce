remove(list = ls())

library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)

#======================================================================
#Function to calculate the f1 score and accuracy
#inputs:
#1. prediction:classification result from model 
#2. cell_label:label cell-type for each cell 
#3. cellnames:names of each celltype 
#
#outputs:
#f1 score matrix: f1 score for each cell-type, mean of f1 score,and accuracy
#========================================================================
eva.cal=function(prediction,cell_label,cellnames){
  F1=matrix(NA,nrow = 1,ncol = length(cellnames)+2)
  colnames(F1)=c(cellnames,"mean","accuracy")
  for(jj in 1:length(cellnames)){
    tp=length(which(cell_label==cellnames[jj]&prediction==cellnames[jj]))#true positive
    fp=length(which(cell_label!=cellnames[jj]&prediction==cellnames[jj]))#false positive
    fn=length(which(cell_label==cellnames[jj]&prediction!=cellnames[jj]))#false negative
    precision=tp/(tp+fp)#precision
    recall=tp/(tp+fn)#recall
    F1[1,jj]=2*(precision*recall)/(precision+recall)
  }
  #average F1 error rate
  F1[1,][is.nan(F1[1,])]=0
  F1[1,length(cellnames)+1]=mean(F1[1,1:length(cellnames)])
  
  #accuracy:
  lab_c=table(prediction==cell_label)["TRUE"]
  accuracy=lab_c/length(cell_label)
  F1[1,length(cellnames)+2]=accuracy
  return(F1)
}


#1.reading label prediction for all methods in BarH dataset
methods.names=c("scAnnotate","singleCellNet","scPred",
                "scClassify","CaSTLe","SingleR",
                "scmapCluster","scmapCell","CHETAH",
                "scID")

label_list=list()

for(jj in 1:10){
  temp=readRDS(paste0("labels_",methods.names[1],"_train_",jj,"_test_",jj,".rds"))
  temp=temp[,c(1,2)]
  label_list[[jj]]=matrix(data=NA,nrow = nrow(temp),ncol = length(methods.names)+1)
  colnames(label_list[[jj]])=c("cell_label",methods.names)
  label_list[[jj]][,1]=temp[,1]
  label_list[[jj]][,2]=temp[,2]
  for(kk in 2:length(methods.names)){
    temp=readRDS(paste0("labels_",methods.names[kk],"_train_",jj,"_test_",jj,".rds"))
    label_list[[jj]][,kk+1]=temp[,2]
    print(identical(temp[,1],temp[,1]))
  }
}

#saveRDS(label_list,file = "CV_label_predicition.rds")
#2. performance for each methods:

perform.list=list()

for(ii in 1:10){
  test_cellnames=names(table(label_list[[ii]][,1]))
  temp=eva.cal(prediction = label_list[[ii]][,2],cell_label = label_list[[ii]][,1],cellnames = test_cellnames)
  perform.list[[ii]]=matrix(data=NA,nrow = length(methods.names),ncol = ncol(temp))
  rownames(perform.list[[ii]])=methods.names
  colnames(perform.list[[ii]])=colnames(temp)
  perform.list[[ii]][1,]=temp
  for(kk in 2:length(methods.names)){
    temp=eva.cal(prediction = label_list[[ii]][,kk+1],cell_label = label_list[[ii]][,1],cellnames = test_cellnames)
    perform.list[[ii]][kk,]=temp
  }
}


saveRDS(perform.list,file="CV_performance.rds")

