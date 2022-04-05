remove(list = ls())

#install.packages("GEOquery")

###1.download and subset the dataset
#=====================================
library(GEOquery)
geo.acc <- "GSE84133"
untar(paste0(geo.acc,"_RAW.tar"))
files<-dir(pattern="gz$")
sapply(files,gunzip)

#1) human 
filelist <- dir(pattern="human.*\\.csv$")
df=read.csv(filelist[1],row.names = NULL)
for(ii in 2:length(filelist)){
  temp=read.csv(filelist[ii],row.names = NULL)
  df=rbind(df,temp)
}

rownames(df)=make.names(df[,2],unique = TRUE)
df=df[,-c(1,2)]

#8569 samples, 20125 feature genes + 1 assigned_cluster

saveRDS(df,"Baron_human.rds")


remove(list = ls())
#2) mouse 
filelist <- dir(pattern="mouse.*\\.csv$")
df=read.csv(filelist[1])
for(ii in 2:length(filelist)){
  temp=read.csv(filelist[ii])
  df=rbind(df,temp)
}

rownames(df)=make.names(df[,2],unique = TRUE)
df=df[,-c(1,2)]

#1886 samples, 14878 feature genes + 1 assigned_cluster

saveRDS(df,"Baron_mouse.rds")



remove(list = ls())
###2. preprocessing
#1)normalization
#2)Load the ortholog table: convert test data (human) gene names to mouse ortholog names
#  and limit analysis to genes in common between the training and test data.
#3)top 2000 highly variable genes
#==================================
library("singleCellNet")
library("Seurat")

# !!! 1. Load Data
#=================================================
#read datasets
train_dat=readRDS("Baron_mouse.rds")
test_dat=readRDS("Baron_human.rds")

#!!! 2. pre-processing data
#=================================================
#!!!a)normalization
#transfer matrix to SingleCellAssay
reference <- CreateSeuratObject(counts=as.matrix(t(train_dat[,-1])))
reference <- NormalizeData(reference, normalization.method = "LogNormalize")

expTrain=GetAssayData(reference)
expTrain=as.matrix(expTrain)
train.yy=train_dat[,1]

#!!!b)testing dataset normalization
qury=CreateSeuratObject(counts=as.matrix(t(test_dat[,-1])))
qury=NormalizeData(qury, normalization.method = "LogNormalize")

expQuery=GetAssayData(qury)
expQuery=as.matrix(expQuery)
test.yy=test_dat[,1]

#!!!c)Load the ortholog table: convert gene
oTab = utils_loadObject("human_mouse_genes_Jul_24_2018.rda")
aa = csRenameOrth(expQuery, expTrain, oTab)
expTrainOrth = aa[['expTrain']]
expQueryOrth = aa[['expQuery']]


B_mouse=data.frame(yy=train.yy,t(expTrainOrth))

B_human=data.frame(yy=test.yy,t(expQueryOrth))


#make the T_cell are same label on both humna and mouse dataset.
B_human[which(B_human[,1]=="t_cell"),1]="T_cell"

saveRDS(B_mouse,"panc_mouse.rds")
saveRDS(B_human,"panc_human.rds")
