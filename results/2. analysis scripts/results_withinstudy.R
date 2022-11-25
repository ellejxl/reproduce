remove(list = ls())

library(stringr)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

file.name=list.files(pattern = "median_")
acc.oup.dat=readRDS(file = file.name[1])

for(ii in 2:length(file.name)){
  temp=readRDS(file = file.name[ii])
  acc.oup.dat=rbind(acc.oup.dat,temp)
}

rownames(acc.oup.dat)=c("ALM",
                        "Baron (Human)",
                        "Baron (Mouse)",
                        "CellBench 10X",
                        "CellBench Cel-seq2",
                        "MTG",
                        "PBMC.10Xv2",
                        "PBMC.10Xv3",
                        "PBMC.CS",
                        "PBMC.DS",
                        "PBMC.ID",
                        "PBMC.SS",
                        "PBMC.SW",
                        "VISp")

acc.oup.dat=acc.oup.dat[c(7,8,10,13,11,12,9,4,5,14,1,6,3,2),]


acc.oup.dat[is.na(acc.oup.dat)]=0

methods.names=colnames(acc.oup.dat)

#3.median, mean, and range of performance on scAnnotate
round(colMeans(acc.oup.dat),4)
summary(acc.oup.dat)$median

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
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_x_discrete(labels=function(x) str_wrap(x,width = 20))

ggsave("withinstudy_acc.png",width = 25, height = 15, units = "cm")   


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
  theme(text = element_text(size=20),legend.text = element_text(size=11),legend.title = element_text(size=11))+
  coord_flip()

ggsave("violin_withinstudy_acc.png",width = 8, height = 12, units = "cm")


