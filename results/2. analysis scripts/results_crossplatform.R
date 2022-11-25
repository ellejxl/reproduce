remove(list = ls())

library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)
library(ggdendro)

set.seed(0)
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


#1.reading label prediction for all methods in crossplatform dataset
methods.names=c("scAnnotate","singleCellNet","scPred",
                "scClassify","CaSTLe","SingleR",
                "scmapCluster","scmapCell","CHETAH",
                "scID")

#1.1) PBMC dataset
train_dat_list=c("pbmc_10x-v3", "pbmc_10x-v2", "pbmc_Smart-seq2", "pbmc_Seq-Well", "pbmc_inDrops", "pbmc_Drop-seq", "pbmc_CEL-Seq2")
test_dat_list=c("pbmc_10x-v3", "pbmc_10x-v2", "pbmc_Smart-seq2", "pbmc_Seq-Well", "pbmc_inDrops", "pbmc_Drop-seq", "pbmc_CEL-Seq2")

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
train_dat_list=c("cellbench_celseq2","cellbench_10x")
test_dat_list=c("cellbench_celseq2","cellbench_10x")

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

saveRDS(label_list,file = "crossplatform_label_predicition.rds")
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

saveRDS(perform.list,file="crossplatform_performance.rds")


#3.acc summary
acc.oup.dat=matrix(NA,ncol =length(methods.names),nrow = length(dat_combine_list) )
colnames(acc.oup.dat)=methods.names
#tissue=c(rep(" (PBMC)",length(dat_combine_list)-length(train_dat_list)),
#         rep("(Lung Cancer)",length(train_dat_list)))
#rownames(acc.oup.dat)=paste0(str_split_fixed(dat_combine_list,"_",4)[,2],", ",str_split_fixed(dat_combine_list,"_",4)[,4],tissue)
rownames(acc.oup.dat)=paste0(str_split_fixed(dat_combine_list,"_",4)[,2],", ",str_split_fixed(dat_combine_list,"_",4)[,4])
rownames(acc.oup.dat)[43]="CEL-Seq2, 10x (*)"
rownames(acc.oup.dat)[44]="10x, CEL-Seq2 (*)"

for(ii in 1:length(dat_combine_list)){
  for(kk in 1:length(methods.names)){
    acc.oup.dat[ii,kk]=perform.list[[ii]][kk,ncol(perform.list[[ii]])]
  }
}



acc.oup.dat[is.na(acc.oup.dat)]=0

#4. test significant difference by Wilcoxon-Sum test
wilcox.matrix=matrix(NA,nrow = length(methods.names)-1,ncol =1)
colnames(wilcox.matrix)=methods.names[1]
rownames(wilcox.matrix)=methods.names[-1]
for(ii in 1:nrow(wilcox.matrix)){
  wilcox.matrix[ii,1]=wilcox.test(acc.oup.dat[,1],acc.oup.dat[,ii+1],paired = TRUE)$p.value
}


#5.make plots:
inv.oup=-(acc.oup.dat)
Ranking=t(apply(inv.oup,1,min_rank))
Ranking=round(Ranking,0)

Ranking[which(is.na(acc.oup.dat))]=NA

df=as.data.frame(Ranking)
df$row.names=rownames(df)
long.df=melt(df,id=c("row.names"))

df2=as.data.frame(acc.oup.dat)
df2$row.names=rownames(df2)
long.df.2=melt(df2,id=c("row.names"))
colnames(long.df.2)=c("row.names","variable","Accuracy")

df.oup=full_join(long.df,long.df.2,by=c("row.names","variable"))
colnames(df.oup)=c("dataset","methods","Ranking","Accuracy")
df.oup$Ranking=factor(df.oup$Ranking)

#specific the order of y axis
df.oup$dataset=as.character(df.oup$dataset) #turn the variable into a character vector
df.oup$dataset=factor(df.oup$dataset,levels = unique(df.oup$dataset) ) #then turn it back into a factor with the levels in the correct order

#ranking plot
ggplot(df.oup, aes(x=dataset,y = methods))+
  theme_classic()+
  geom_point(aes(col=Ranking,size=Accuracy))+
  labs(x="Datasets (training, test)")+
  scale_colour_manual(values = c("#D30000","#FF2800","#FA8072",
                                 "#FF9529","#FFDF00",
                                 "#0171B9","#AFDCEC",
                                 "#107C10","#3CB043","#d0f0c0"))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(text = element_text(size=19),legend.title = element_text(size=19),legend.position = "left")+
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave("crossplatform_acc.png",width = 50, height = 17, units = "cm")  

#geom_violin 1: compared by each methods

ggplot(long.df.2,aes(x=variable,y=Accuracy))+
  geom_violin()+
  labs(y="Classification accuracy",
       x="methods")+
#  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.05)+
#  geom_jitter(shape=16, position=position_jitter(0.2))+
  theme_bw()+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(axis.text.y=element_blank())+
  theme(text = element_text(size=20),legend.text = element_text(size=11),legend.title = element_text(size=11))+
  coord_flip()


ggsave("violin_crossplatform_acc.png",width = 12, height = 12, units = "cm")

#geom_violin 2: compared by each dataset combination
ggplot(df.oup, aes(x=dataset,y = Accuracy))+
  geom_violin()+
  labs(y="Classification accuracy",
       x="dataset")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.35))+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(text = element_text(size=22),legend.title = element_text(size=22))

ggsave("violin_crossplatform_dataset.png",width = 50, height = 15, units = "cm")   

#6. ensemble plot

#heatmap functions:

#only well performed methods
for(pp in 1:length(label_list)){
  
  selected.results=label_list[[pp]][,c(1:7)]
  
  tf.matrix=matrix(NA,nrow = nrow(selected.results),ncol = ncol(selected.results)-1)
  colnames(tf.matrix)=methods.names[1:6]
  
  for(ii in 1:ncol(tf.matrix)){
    tf.matrix[which(selected.results[,1]==selected.results[,ii+1]),ii]="1"
    tf.matrix[which(selected.results[,1]!=selected.results[,ii+1]),ii]="0"
  }
  
  tf.matrix=apply(tf.matrix,2,as.numeric)
  
  tf.matrix=tf.matrix[-which(rowSums(tf.matrix)==0),]
  tf.matrix=tf.matrix[-which(rowSums(tf.matrix)==ncol(tf.matrix)),]
  
  
  heatmap(t(tf.matrix),
          scale = "none",
          col=c("black","white"))
}


#ggplot


#only well performed methods

#add dataset name in each plot
dataset.names=names(label_list)
dataset.names
dataset.names=str_replace_all(dataset.names,"pbmc_10x-v3","PBMC.10Xv3")
dataset.names=str_replace_all(dataset.names,"pbmc_10x-v2","PBMC.10Xv2")
dataset.names=str_replace_all(dataset.names,"pbmc_Smart-seq2","PBMC.SS")
dataset.names=str_replace_all(dataset.names,"pbmc_Seq-Well","PBMC.SW")
dataset.names=str_replace_all(dataset.names,"pbmc_inDrops","PBMC.ID")
dataset.names=str_replace_all(dataset.names,"pbmc_Drop-seq","PBMC.DS")
dataset.names=str_replace_all(dataset.names,"pbmc_CEL-Seq2","PBMC.CS")
dataset.names=str_replace_all(dataset.names,"_",", ")

for(pp in 1:length(label_list)){
  
  #clustering all
  selected.results=label_list[[pp]][,c(1:7)]
  
  tf.matrix=matrix(NA,nrow = nrow(selected.results),ncol = ncol(selected.results)-1)
  colnames(tf.matrix)=methods.names[1:6]
  
  for(ii in 1:ncol(tf.matrix)){
    tf.matrix[which(selected.results[,1]==selected.results[,ii+1]),ii]="1"
    tf.matrix[which(selected.results[,1]!=selected.results[,ii+1]),ii]="0"
  }
  
  tf.matrix=apply(tf.matrix,2,as.numeric)
  
  tf.matrix=tf.matrix[-which(rowSums(tf.matrix)==0),]
  tf.matrix=tf.matrix[-which(rowSums(tf.matrix)==ncol(tf.matrix)),]
  
  # k-means clustering of samples
  set.seed(0)
  km=kmeans(tf.matrix,centers = 20)
  # sort the order of clustering by proportion of false prediction by scAnnotate
  sort.cluster=matrix(NA,nrow=length(km$cluster),ncol = 1)
  colnames(sort.cluster)="sort_cluster"
  cc.correct=matrix(NA,nrow = length(table(km$cluster)),ncol = 2) #proprtion of false prediction by scAnnotate for each clustering
  colnames(cc.correct)=c("proprotion_false","sort_cluster")
  for(cc in 1:nrow(cc.correct)){
    cc.correct[cc,1]=sum(tf.matrix[which(km$cluster==cc),1]==0)/sum(km$cluster==cc)
  }
  
  cc.correct[,2]=order(cc.correct[,1],decreasing = TRUE)
  
  
  for(cc in 1:nrow(cc.correct)){
    sort.cluster[which(km$cluster==cc.correct[cc,2]),1]=cc
  }
  

  sample.idx=order(sort.cluster[,1])
  
  
  tf.matrix.sort=tf.matrix[sample.idx,]
  
  df3=as.data.frame(tf.matrix.sort)
  
  for(ii in 1:ncol(df3)){
    df3[,ii]=factor(df3[,ii])
  }
  
  df3$row.names=rownames(df3)
  #reorder factors
  df3$row.names=factor(df3$row.names,levels = df3$row.names)
  
  long.df.3=melt(df3,id=c("row.names"))
  colnames(long.df.3)=c("sample","methods","annotation")
  
  long.df.3$annotation=str_replace_all(long.df.3$annotation,"0","False")
  long.df.3$annotation=str_replace_all(long.df.3$annotation,"1","True")
  
  ggplot(long.df.3,aes(x=sample,y=methods,fill=annotation))+
    geom_tile()+
    labs(title=paste0("dataset (training, test): ",dataset.names[pp]),y="methods",x="Cell")+
    theme(axis.text.x = element_blank())+
    theme(axis.title.y = element_blank())+
    scale_fill_manual(values = c("black","#e5e1e1"))+
    theme(text = element_text(size=22),legend.title = element_text(size=22))
  
  ggsave(paste0("ensemble_top_",pp,".png"),width = 30, height = 8, units = "cm") 
  
}

for(pp in 19){
  #clustering all
  selected.results=label_list[[pp]][,c(1:7)]
  
  tf.matrix=matrix(NA,nrow = nrow(selected.results),ncol = ncol(selected.results)-1)
  colnames(tf.matrix)=methods.names[1:6]
  
  for(ii in 1:ncol(tf.matrix)){
    tf.matrix[which(selected.results[,1]==selected.results[,ii+1]),ii]="1"
    tf.matrix[which(selected.results[,1]!=selected.results[,ii+1]),ii]="0"
  }
  
  tf.matrix=apply(tf.matrix,2,as.numeric)
  
  all.t=length(which(rowSums(tf.matrix)==0))
  all.f=length(which(rowSums(tf.matrix)==ncol(tf.matrix)))
  sample.select=nrow(tf.matrix)
  
  tf.matrix=tf.matrix[-which(rowSums(tf.matrix)==0),]
  tf.matrix=tf.matrix[-which(rowSums(tf.matrix)==ncol(tf.matrix)),]
  
  some.t=nrow(tf.matrix)
  
  hc=hclust(dist(t(tf.matrix)),"ave")
  ggdendrogram(hc)+
    theme(text = element_text(size=22),legend.title = element_text(size=22))+
    coord_flip()+
    scale_y_reverse()
  
  ggsave(paste0("dendrogram_",pp,".png"),width = 5, height = 8, units = "cm")
  
  # k-means clustering of samples
  set.seed(0)
  km=kmeans(tf.matrix,centers = 20)
  # sort the order of clustering by proportion of false prediction by scAnnotate
  sort.cluster=matrix(NA,nrow=length(km$cluster),ncol = 1)
  colnames(sort.cluster)="sort_cluster"
  cc.correct=matrix(NA,nrow = length(table(km$cluster)),ncol = 2) #proprtion of false prediction by scAnnotate for each clustering
  colnames(cc.correct)=c("proprotion_false","sort_cluster")
  for(cc in 1:nrow(cc.correct)){
    cc.correct[cc,1]=sum(tf.matrix[which(km$cluster==cc),1]==0)/sum(km$cluster==cc)
  }
  
  cc.correct[,2]=order(cc.correct[,1],decreasing = TRUE)
  
  
  for(cc in 1:nrow(cc.correct)){
    sort.cluster[which(km$cluster==cc.correct[cc,2]),1]=cc
  }
  
  
  sample.idx=order(sort.cluster[,1])
  
  
  tf.matrix.sort=tf.matrix[sample.idx,]
  
  df3=as.data.frame(tf.matrix.sort)
  
  for(ii in 1:ncol(df3)){
    df3[,ii]=factor(df3[,ii])
  }
  
  df3$row.names=rownames(df3)
  #reorder factors
  df3$row.names=factor(df3$row.names,levels = df3$row.names)
  
  long.df.3=melt(df3,id=c("row.names"))
  colnames(long.df.3)=c("sample","methods","annotation")
  
  long.df.3$annotation=str_replace_all(long.df.3$annotation,"0","False")
  long.df.3$annotation=str_replace_all(long.df.3$annotation,"1","True")
  
  #reorder the factors
  if(pp==19){
    long.df.3$methods=factor(long.df.3$methods,levels = c("scAnnotate","SingleR",
                                                          "CaSTLe","scClassify",
                                                          "singleCellNet","scPred"))
  }
  
  ggplot(long.df.3,aes(x=sample,y=methods,fill=annotation))+
    geom_tile()+
    theme(axis.text.x = element_blank())+
    scale_fill_manual(values = c("black","#e5e1e1"))+
    theme(text = element_text(size=22),legend.title = element_text(size=22))+
    theme(axis.title.y = element_blank())+
    theme(axis.title.x = element_blank())
  
  ggsave(paste0("ensemble_top_",pp,".png"),width = 30, height = 8, units = "cm") 
}




