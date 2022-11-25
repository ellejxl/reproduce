remove(list = ls())


library("MTPS")
library("Seurat")
set.seed(0)
#load dataset
input=readRDS("cellbench_10x.rds")

#arrange the dataset with cell type
input=input[order(input[,1]),]

celltype=input[,1]


#straitified sampling
#80% of data are training(idx==0), 20% of data are test(idx==1).
sample.eachtype=table(input[,1])
cellnames=names(table(input[,1]))
idx.matrix=matrix(NA,nrow = nrow(input),ncol = 10)
for(ii in 1:ncol(idx.matrix)){
  idx=vector()
  for(jj in 1:length(cellnames)){
    temp=vector()
    temp[1:sample.eachtype[jj]]=createFolds(input[,1][1:sample.eachtype[jj]],k=5,list = F)
    idx=c(idx,temp)
  }# out for jj
  idx.matrix[,ii]=idx
}

idx.matrix[which(idx.matrix!=1)]=0


for(ii in 1:ncol(idx.matrix)){
  idx=idx.matrix[,ii]
  
  train_dat=input[which(idx==0),]
  test_dat=input[which(idx==1),]
  saveRDS(train_dat,paste0("train_",ii,".rds"))
  saveRDS(test_dat,paste0("test_",ii,".rds"))
  
}
