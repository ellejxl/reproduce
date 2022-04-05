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
  colnames(F1)=c(cellnames,"mean F1","accuracy")
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


#1.reading label prediction for all methods in crossspecies dataset
methods.names=c("scAnnotate","singleCellNet","scPred",
                "scClassify","CaSTLe","SingleR",
                "scmapCluster","scmapCell","CHETAH",
                "scID")

#1.1) pancrease
train_dat_list=c("panc_human","panc_mouse")
test_dat_list=c("panc_human","panc_mouse")
label_list=list()
for(ii in 1:length(train_dat_list)){
  for(jj in 1:length(test_dat_list)){
    if(train_dat_list[ii]!=test_dat_list[jj]){
      dat_combine=paste0(train_dat_list[ii],"_",test_dat_list[jj])
      temp=readRDS(paste0("labels_",methods.names[1],"_",dat_combine,".rds"))
      label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]]=matrix(data=NA,nrow = nrow(temp),ncol = length(methods.names)+1)
      colnames(label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]])=c("cell_label",methods.names)
      label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]][,1]=temp[,1]
      label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]][,2]=temp[,2]
      for(kk in 2:length(methods.names)){
        temp=readRDS(paste0("labels_",methods.names[kk],"_",dat_combine,".rds"))
        label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]][,kk+1]=temp[,2]
        print(identical(label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]][,1],temp[,1]))
      }
    }
  }
}


#1.2) cellbench dataset (lung cancer cell)
train_dat_list=c("V13","ALM3","MTG3")
test_dat_list=c("V13","ALM3","MTG3")

for(ii in 1:length(train_dat_list)){
  for(jj in 1:length(test_dat_list)){
    if(train_dat_list[ii]!=test_dat_list[jj]){
      dat_combine=paste0(train_dat_list[ii],"_",test_dat_list[jj])
      if(dat_combine!="V13_ALM3"&dat_combine!="ALM3_V13"){
        temp=readRDS(paste0("labels_",methods.names[1],"_",dat_combine,".rds"))
        label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]]=matrix(data=NA,nrow = nrow(temp),ncol = length(methods.names)+1)
        colnames(label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]])=c("cell_label",methods.names)
        label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]][,1]=temp[,1]
        label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]][,2]=temp[,2]
        for(kk in 2:length(methods.names)){
          temp=readRDS(paste0("labels_",methods.names[kk],"_",dat_combine,".rds"))
          label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]][,kk+1]=temp[,2]
          print(identical(label_list[[paste0(train_dat_list[ii],"_",test_dat_list[jj])]][,1],temp[,1]))
        }
      }
    }
  }
}

#saveRDS(label_list,file = "crossspecies_label_predicition.rds")
#2. performance for each methods:

perform.list=list()
dat_combine_list=names(label_list)
for(ii in 1:length(dat_combine_list)){
  test_cellnames=names(table(label_list[[ii]][,1]))
  temp=eva.cal(prediction = label_list[[ii]][,2],cell_label = label_list[[ii]][,1],cellnames = test_cellnames)
  perform.list[[dat_combine_list[ii]]]=matrix(data=NA,nrow = length(methods.names),ncol = ncol(temp))
  rownames(perform.list[[dat_combine_list[ii]]])=methods.names
  colnames(perform.list[[dat_combine_list[ii]]])=colnames(temp)
  perform.list[[dat_combine_list[ii]]][1,]=temp
  for(kk in 2:length(methods.names)){
    temp=eva.cal(prediction = label_list[[ii]][,kk+1],cell_label = label_list[[ii]][,1],cellnames = test_cellnames)
    perform.list[[dat_combine_list[ii]]][kk,]=temp
  }
}

saveRDS(perform.list,file="crossspecies_performance.rds")
