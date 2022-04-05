remove(list = ls())

#install packages
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("scater")
#install.packages("xgboost")
#install.packages("igraph")

library(scater)  
library(xgboost) 
library(igraph) 
set.seed(0)

args=commandArgs(trailingOnly = TRUE)
train_inp=as.character(args[1])
test_inp=as.character(args[2])

print(paste0("============ CaSTLe: train=",train_inp," test=",test_inp,"======"))


runCastlePerCellType = function(train_dat, test_dat){
  # 1. Load datasets in scater format: loaded files expected to contain "Large SingleCellExperiment" object
  source = t(train_dat[,-1])
  target = t(test_dat[,-1])
  ds1 = train_dat[,-1] 
  ds2 = test_dat[,-1]
  sourceCellTypes = as.factor(train_dat[,1])
  targetCellTypes = as.factor(test_dat[,1])
  
  # remove unlabeled from source
  unknownLabels = levels(sourceCellTypes)[grep("not applicable|unclassified|contaminated|unknown", levels(sourceCellTypes))]
  if (length(unknownLabels)>0) {
    hasKnownLabel = is.na(match(sourceCellTypes, unknownLabels))
    sourceCellTypes = sourceCellTypes[hasKnownLabel]
    sourceCellTypes = as.factor(as.character(sourceCellTypes))
    ds1 = ds1[hasKnownLabel,]
  }
  
  # 2. Unify sets, excluding low expressed genes
  source_n_cells_counts = apply(source, 1, function(x) { sum(x > 0) } )
  target_n_cells_counts = apply(target, 1, function(x) { sum(x > 0) } )
  common_genes = intersect( colnames(ds1)[source_n_cells_counts>10], 
                            colnames(ds2)[target_n_cells_counts>10]
  )
  remove(source_n_cells_counts, target_n_cells_counts)
  ds1 = ds1[, colnames(ds1) %in% common_genes]
  ds2 = ds2[, colnames(ds2) %in% common_genes]
  ds = rbind(ds1[,common_genes], ds2[,common_genes])
  isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
  remove(ds1, ds2)
  
  # 3. Highest mean in both source and target
  topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
  
  L = length(levels(sourceCellTypes))
  summary=data.frame(InstancesInSource=character(L),InstancesInTarget=character(L),Accuracy=character(L),Sensitivity=character(L),Specificity=character(L),AUC=character(L), stringsAsFactors = FALSE)
  rownames(summary) = levels(sourceCellTypes)
  
  # for each cell - what is the most probable classification?
  targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(sourceCellTypes))
  
  
  # iterate over all source cell types
  for (cellType in levels(sourceCellTypes)) {
    
    inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))
    
    # 4. Highest mutual information in source
    topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))
    
    # 5. Top n genes that appear in both mi and avg
    selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
    
    # 6. remove correlated features
    tmp = cor(ds[,selectedFeatures], method = "pearson")
    tmp[!lower.tri(tmp)] = 0
    selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
    remove(tmp)
    
    # 7,8. Convert data from continous to binned dummy vars
    # break datasets to bins
    dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
    # use only bins with more than one value
    nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
    # convert to dummy vars
    ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
    remove(dsBins, nUniq)
    
    cat(paste0("<h2>Classifier for ",cellType,"</h2>"))
    
    inTypeSource = sourceCellTypes == cellType
    # 9. Classify
    xg=xgboost(data=ds0[isSource,] , 
               label=inTypeSource,
               objective="binary:logistic", 
               eta=0.7 , nthread=1, nround=20, verbose=0,
               gamma=0.001, max_depth=5, min_child_weight=10)
    
    # 10. Predict
    inTypeProb = predict(xg, ds0[!isSource, ])
    inTypePred = factor(ifelse(inTypeProb > 0.5,1,0), 0:1)
    
    targetClassification[cellType,] = inTypeProb
    
  }# out  for (cellType in levels(sourceCellTypes))
  
  True_Labels_Castle <- list(test_dat[,1])
  Pred_Labels_Castle <- list(rownames(targetClassification)[apply(targetClassification,2,which.max)])
  
  True_Labels_Castle <- as.vector(unlist(True_Labels_Castle))
  Pred_Labels_Castle <- as.vector(unlist(Pred_Labels_Castle))
  
  labels=data.frame(cell_label=True_Labels_Castle,
                    prediction=Pred_Labels_Castle)
  
  return(labels)
}# out for function runCastlePerCellType


BREAKS=c(-1, 0, 1, 6, Inf)
nFeatures = 100

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


labels=runCastlePerCellType(train_dat = train_dat,
                            test_dat = test_dat)

saveRDS(labels,paste0("labels_CaSTLe_",train_inp,"_",test_inp,".rds"))

