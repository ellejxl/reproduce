remove(list = ls())

library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)

#======================================================================
#Function to calculate the accuracy of cell population and overall auucracy
#inputs:
#1. prediction:classification result from model 
#2. cell_label:label cell-type for each cell 
#3. cellnames:names of each celltype 
#
#outputs:
#accuracy score matrix: accuracy score for each cell-type, mean of accuracy score,and accuracy
#========================================================================
eva.cal=function(prediction,cell_label,cellnames){
  acc=matrix(NA,nrow = 1,ncol = length(cellnames)+2)
  colnames(acc)=c(cellnames,"mean of accuracy","overall accuracy")
  for(jj in 1:length(cellnames)){
    tp=length(which(cell_label==cellnames[jj]&prediction==cellnames[jj]))#true positive
    acc[1,jj]=tp/length(which(cell_label==cellnames[jj]))
  }
  #average acc error rate
  acc[1,][is.nan(acc[1,])]=0
  acc[1,length(cellnames)+1]=mean(acc[1,1:length(cellnames)])
  
  #accuracy:
  lab_c=table(prediction==cell_label)["TRUE"]
  accuracy=lab_c/length(cell_label)
  acc[1,length(cellnames)+2]=accuracy
  return(acc)
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

saveRDS(label_list,file = "CV_label_predicition.rds")
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




#3.acc summary
acc.oup.dat=matrix(NA,ncol =length(methods.names),nrow = length(label_list) )
colnames(acc.oup.dat)=methods.names
rownames(acc.oup.dat)=c(1:10)

for(ii in 1:10){
  for(kk in 1:length(methods.names)){
    acc.oup.dat[ii,kk]=perform.list[[ii]][kk,ncol(perform.list[[ii]])]
  }
}



#acc.oup.dat[is.na(acc.oup.dat)]=0


df=as.data.frame(acc.oup.dat)
df$row.names=rownames(df)
long.df=melt(df,id=c("row.names"))
colnames(long.df)=c("replications","Methods","Accuracy")

#boxplot
ggplot(long.df,aes(x=Methods,y=Accuracy))+
  geom_boxplot(colour="grey50")+
  labs(y="Classification accuracy",
       x="methods")+
  theme_bw()+
  theme(axis.title.y = element_blank())+
  theme(axis.title.x = element_blank())+
  theme(text = element_text(size=22),legend.text = element_text(size=12),legend.title = element_text(size=12))+
  coord_flip()

ggsave("boxplot_cv_acc.png",width = 25, height = 25, units = "cm")

#4. some analysis
mean.matrix=matrix(NA,nrow = length(perform.list),ncol = ncol(perform.list[[1]]))
colnames(mean.matrix)=colnames(perform.list[[1]])
for(ii in 1:nrow(mean.matrix)){
  mean.matrix[ii,]=perform.list[[ii]][1,]
}
round(colMeans(mean.matrix),4)

