remove(list = ls())

library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

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

saveRDS(label_list,file = "crossspecies_label_predicition.rds")
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


#3.acc summary
acc.oup.dat=matrix(NA,ncol =length(methods.names),nrow = length(dat_combine_list) )
colnames(acc.oup.dat)=methods.names
dat_combine_list
rownames(acc.oup.dat)=c("human, mouse (Pancreas)",
                        "mouse, human (Pancreas)",
                        "mouse[VISp], human (Brain)",
                        "mouse[ALM], human (Brain)",
                        "human, mouse[VISp] (Brain)",
                        "human, mouse[ALM] (Brain)")

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
  wilcox.matrix[ii,1]=wilcox.test(acc.oup.dat[,1],acc.oup.dat[,ii+1])$p.value
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
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_x_discrete(labels=function(x) str_wrap(x,width = 19))

ggsave("crossspecies_acc.png",width = 25, height = 16, units = "cm")  


#violin plot
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
  theme(text = element_text(size=18),legend.text = element_text(size=10),legend.title = element_text(size=11))+
  coord_flip()
  


ggsave("violin_crossspecies_acc.png",width = 8, height = 12, units = "cm")


