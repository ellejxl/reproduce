# Reproduce

This repository contains scripts and information to reproduce all results in our manuscript.

# About this project 

We build an automated cell annotation tool based on supervised machine learning
algorithms. We use a marginal mixture model to describe both the dropout proportion and the non-dropout expression level distribution of a gene. We develop a marginal model based ensemble learning approach to avoid having to specify and estimate a high-dimensional joint distribution for all genes. First, we build a ‘weak’ classifier using the mixture model for each gene. Then, we combine ‘weak’ classifiers of all genes into a single ‘strong’ classifier to annotate cells.

---

# Directory Layout

<details><summary>code</summary>

    ├── code                                               # All of scripts to reproduce the results in scAnnotate manuscript.
    │     ├── 1. PBMC
    │            ├── preprocessing_PBMC_crossplatform.sh   # Shell script to run the preprocessing_PBMC_crossplatform.R
    │            ├── preprocessing_PBMC_crossplatform.R    # R script to preprocessing the PBMC dataset
    │            ├── run_pbmc.sh                           # Shell script with loops to submit all shell scripts for each methods on PBMC dataset at once
    │            ├── run.sh                                # shell script to run each corresponding R script on selected dataset
    │            ├── scAnnotate.r                          # R scripts to run our method scAnnotate with selected dataset
    │            ├── CaSTLe.r                              # R scripts to run the competing method CaSTLe with selected dataset
    │            ├── CHETAH.r                              # R scripts to run the competing method CHETAH with selected dataset
    │            ├── scClassify.r                          # R scripts to run the competing method scClassify with selected dataset
    │            ├── scID.r                                # R scripts to run the competing method scID with selected dataset
    │            ├── scmapCell.r                           # R scripts to run the competing method scmapCell with selected dataset
    │            ├── scmapCluster.r                        # R scripts to run the competing method scmapCluster with selected dataset
    │            ├── scPred.r                              # R scripts to run the competing method scPred with selected dataset
    │            ├── singleCellNet.r                       # R scripts to run the competing method singleCellNet with selected dataset
    │            └── SingleR.r                             # R scripts to run the competing method SingleR with selected dataset
    |
    │     ├── 2. CellBench
    │            ├── preprocessing_CellBench.sh            # Shell script to run the preprocessing_CellBench.R
    │            ├── preprocessing_CellBench.R             # R script to preprocessing the CellBench dataset
    │            ├── run_cellbench.sh                      # Shell script with loops to submit all shell scripts for each methods on cellbench dataset at once
    │            ├── run.sh                                # Shell script to run each corresponding R script on selected dataset
    │            ├── scAnnotate.r                          # R scripts to run our method scAnnotate with selected dataset
    │            ├── CaSTLe.r                              # R scripts to run the competing method CaSTLe with selected dataset
    │            ├── CHETAH.r                              # R scripts to run the competing method CHETAH with selected dataset
    │            ├── scClassify.r                          # R scripts to run the competing method scClassify with selected dataset
    │            ├── scID.r                                # R scripts to run the competing method scID with selected dataset
    │            ├── scmapCell.r                           # R scripts to run the competing method scmapCell with selected dataset
    │            ├── scmapCluster.r                        # R scripts to run the competing method scmapCluster with selected dataset
    │            ├── scPred.r                              # R scripts to run the competing method scPred with selected dataset
    │            ├── singleCellNet.r                       # R scripts to run the competing method singleCellNet with selected dataset
    │            └── SingleR.r                             # R scripts to run the competing method SingleR with selected dataset  
    |
    │     ├── 3. Brain
    │            ├── 3.1. ALM_MTG
    |                    ├── preprocessing_ALM_MTG.sh      # Shell script to run the preprocessing_ALM_MTG.R
    │                    ├── preprocessing_ALM_MTG.R       # R script to preprocessing the ALM and MTG dataset
    │                    ├── run_ALM_MTG.sh                # Shell script with loops to submit all shell scripts for each methods on ALM and MTG dataset at once
    │                    ├── run.sh                        # Shell script to run each corresponding R script on selected dataset
    │                    ├── scAnnotate.r                  # R scripts to run our method scAnnotate with selected dataset
    │                    ├── CaSTLe.r                      # R scripts to run the competing method CaSTLe with selected dataset
    │                    ├── CHETAH.r                      # R scripts to run the competing method CHETAH with selected dataset
    │                    ├── scClassify.r                  # R scripts to run the competing method scClassify with selected dataset
    │                    ├── scID.r                        # R scripts to run the competing method scID with selected dataset
    │                    ├── scmapCell.r                   # R scripts to run the competing method scmapCell with selected dataset
    │                    ├── scmapCluster.r                # R scripts to run the competing method scmapCluster with selected dataset
    │                    ├── scPred.r                      # R scripts to run the competing method scPred with selected dataset
    │                    ├── singleCellNet.r               # R scripts to run the competing method singleCellNet with selected dataset
    │                    └── SingleR.r                     # R scripts to run the competing method SingleR with selected dataset   
    |
    │            └── 3.2. V1_MTG   
    |                    ├── preprocessing_V1_MTG.sh      # Shell script to run the preprocessing_V1_MTG.R
    │                    ├── preprocessing_V1_MTG.R       # R script to preprocessing the V1 and MTG dataset
    │                    ├── run_V1_MTG.sh                # Shell script with loops to submit all shell scripts for each methods on V1 and MTG dataset at once
    │                    ├── run.sh                        # Shell script to run each corresponding R script on selected dataset
    │                    ├── scAnnotate.r                  # R scripts to run our method scAnnotate with selected dataset
    │                    ├── CaSTLe.r                      # R scripts to run the competing method CaSTLe with selected dataset
    │                    ├── CHETAH.r                      # R scripts to run the competing method CHETAH with selected dataset
    │                    ├── scClassify.r                  # R scripts to run the competing method scClassify with selected dataset
    │                    ├── scID.r                        # R scripts to run the competing method scID with selected dataset
    │                    ├── scmapCell.r                   # R scripts to run the competing method scmapCell with selected dataset
    │                    ├── scmapCluster.r                # R scripts to run the competing method scmapCluster with selected dataset
    │                    ├── scPred.r                      # R scripts to run the competing method scPred with selected dataset
    │                    ├── singleCellNet.r               # R scripts to run the competing method singleCellNet with selected dataset
    │                    └── SingleR.r                     # R scripts to run the competing method SingleR with selected dataset   
    |
    │     ├── 4. Pancrease
    │            ├── processing_Baron.sh                   # Shell script to run the processing_Baron.R
    │            ├── processing_Baron.R                    # R script to preprocessing the Baron pancrease dataset
    │            ├── run_panc.sh                           # Shell script with loops to submit all shell scripts for each methods on Baron pancrease dataset at once
    │            ├── run.sh                                # Shell script to run each corresponding R script on selected dataset
    │            ├── scAnnotate.r                          # R scripts to run our method scAnnotate with selected dataset
    │            ├── CaSTLe.r                              # R scripts to run the competing method CaSTLe with selected dataset
    │            ├── CHETAH.r                              # R scripts to run the competing method CHETAH with selected dataset
    │            ├── scClassify.r                          # R scripts to run the competing method scClassify with selected dataset
    │            ├── scID.r                                # R scripts to run the competing method scID with selected dataset
    │            ├── scmapCell.r                           # R scripts to run the competing method scmapCell with selected dataset
    │            ├── scmapCluster.r                        # R scripts to run the competing method scmapCluster with selected dataset
    │            ├── scPred.r                              # R scripts to run the competing method scPred with selected dataset
    │            ├── singleCellNet.r                       # R scripts to run the competing method singleCellNet with selected dataset
    │            └── SingleR.r                             # R scripts to run the competing method SingleR with selected dataset 
    │
    │     └── 5. CV_Baron
    │            ├── cv_Baronhuman.sh                      # Shell script to run the cv_Baronhuman.R
    │            ├── cv_Baronhuman.R                       # R script to processing the Baron human pancrease dataset
    │            ├── run_cellbench.sh                      # Shell script with loops to submit all shell scripts for each methods on Baron human pancrease dataset at once
    │            ├── run.sh                                # Shell script to run each corresponding R script on selected dataset
    │            ├── scAnnotate.r                          # R scripts to run our method scAnnotate with selected dataset
    │            ├── CaSTLe.r                              # R scripts to run the competing method CaSTLe with selected dataset
    │            ├── CHETAH.r                              # R scripts to run the competing method CHETAH with selected dataset
    │            ├── scClassify.r                          # R scripts to run the competing method scClassify with selected dataset
    │            ├── scID.r                                # R scripts to run the competing method scID with selected dataset
    │            ├── scmapCell.r                           # R scripts to run the competing method scmapCell with selected dataset
    │            ├── scmapCluster.r                        # R scripts to run the competing method scmapCluster with selected dataset
    │            ├── scPred.r                              # R scripts to run the competing method scPred with selected dataset
    │            ├── singleCellNet.r                       # R scripts to run the competing method singleCellNet with selected dataset
    │            └── SingleR.r                             # R scripts to run the competing method SingleR with selected dataset   
</details>

<details><summary>results</summary>
   
    ├── results                                            # All results used in scAnnotate manuscript.
    │     ├── 1. annotation results
    │            ├── 1.1 PBMC                              # All annotation results for each methods on PBMC dataset
    │            ├── 1.2 CellBench                         # All annotation results for each methods on CellBench dataset 
    │            ├── 1.3 Brain                             # All annotation results for each methods on Brain dataset
    │            ├── 1.4 Pancrease                         # All annotation results for each methods on pancrease dataset
    │            └── 1.5 CV_Baron                          # All annotation results for each methods on Baron human dataset
    |
    │     └── 2. analysis scripts        
    │            ├── results_crossplatform.R               # R script to evaluate the performance of each methods on cross-platform dataset by accuracy of each cell population, mean of accuracy, and overall accuracy
    │            ├── results_crossspecies.R                # R script to evaluate the performance of each methods on cross-species dataset by accuracy of each cell population, mean of accuracy, and overall accuracy
    │            └──results_cv.R                          # R script to evaluate the performance of each methods on intra-dataset (Baron_human) by accuracy of each cell population, mean of accuracy, and overall accuracy 
</details>
  

---

# Package dependencies
- To reproduce the results, you should download and install the required software and package with the corresponding version in the table below:

#### Tools tables:
|software|version|
|--------|-------|
|R|4.1.2|
|scAnnotate|0.0.1|
|scran|1.22.1|
|SingleR|1.8.0|
|singleCellNet|0.1.0|
|scPred|1.9.2|
|Seurat|4.0.5|
|scmap|1.16.0|
|scID|2.2|
|scClassify|1.5.1|
|CHETAH|1.9.0|
|xgboost|1.5.0.2|
|scater|1.22.0|
|SingleCellExperiment|1.16.0|
|igraph|1.2.10|
|magrittr|2.0.1|
|MTPS|1.0.1|
|SeuratData|0.2.1|
|stringr|1.7.6|
|ggplot2|3.3.5|
|reshape2|1.4.4|
|dplyr|1.0.7|
|ggpubr|0.4.0|
|ggdendro|0.1.23|


- scAnnotate will be available to download from CRAN later. Before it is available in CRAN, you can directly download the package file scAnnotate_0.0.1.tar from here.

- To install the package, type the following at the command line:

```
R CMD INSTALL scAnnotate_0.0.1.tar
```

---

# Running files
- We have arrange the scripts by each dataset (1. PBMC 2. CellBench 3. Brain 4. Pancrease 5. CV_Baron), so you could follow the order of each dataset to reproduce all results. 
- All of our scripts were running on the ComputeCanada. Before you submit the shell script, please change the account ID in each shell script (#SBATCH --account=) to your ComputeCanada accout. 
- All of the raw datasets used in our study were publicly available. Due to the size limitation of Github, we didn't upload the raw data files here. You should download the raw dataset directly from its original sources before running the scripts.

<details><summary>1. PBMC </summary>
  
  1). Download the raw data: PBMC dataset was downloaded from the SeuratData package with dataset name "pbmcsca". 
~~~
  R
  library(SeuratData)
  InstallData("pbmcsca")
~~~  
  2). Process the dataset: submit the shell script to processing the PBMC dataset
~~~
 sbatch preprocessing_PBMC_crossplatform.sh
~~~
  3). Run the methods: submit the shell script to run each methods R scripts on selected dataset.
~~~
 sbatch run_pbmc.sh
~~~
</details>

<details><summary>2. CellBench </summary> 
  
  1). Download the raw data: CellBench was download from Zenodo https://zenodo.org/record/3357167#.YjpANOfMKw4
  
  2). Process the dataset: submit the shell script to processing the CellBench dataset
~~~
 sbatch preprocessing_CellBench.sh
~~~
  
  3). Run the methods: submit the shell script to run each methods R scripts on selected dataset.
~~~
 sbatch run_cellbench.sh
~~~
</details>

<details><summary>3. Brain </summary>
  
  3.1 ALM_MTG
  
  1). Download the raw data: CellBench was download from Zenodo https://zenodo.org/record/3357167#.YjpANOfMKw4
  
  2). Process the dataset: submit the shell script to processing the ALM and MTG dataset
~~~
 sbatch preprocessing_ALM_MTG.sh
~~~
  
  3). Run the methods: submit the shell script to run each methods R scripts on selected dataset.
~~~
 sbatch run_ALM_MTG.sh
~~~
  
  3.2 V1_MTG
 
  1). Download the raw data: CellBench was download from Zenodo https://zenodo.org/record/3357167#.YjpANOfMKw4
  
  2). Process the dataset: submit the shell script to processing the V1 and MTG dataset
~~~
 sbatch preprocessing_V1_MTG.sh
~~~
  
  3). Run the methods: submit the shell script to run each methods R scripts on selected dataset.
~~~
 sbatch run_V1_MTG.sh
~~~

</details>


<details><summary>4. Pancrease </summary>
  
  1). Download the raw data: Pancrease data was downloaded from National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) for GSE84133. https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file

  2). Process the dataset: submit the shell script to processing the Baron pancrease dataset
~~~
 sbatch processing_Baron.sh
~~~
  3). Run the methods: submit the shell script to run each methods R scripts on selected dataset.
~~~
 sbatch run_panc.sh
~~~
</details>


<details><summary>5. CV_Baron </summary>
  
  1). move the select the dataset into the folder: selected dataset "Baron_human.rds" was output from processing_Baron.R

  2). Process the dataset: submit the shell script to processing the Baron human pancrease dataset
~~~
 sbatch cv_Baronhuman.sh
~~~
  3). Run the methods: submit the shell script to run each methods R scripts on selected dataset.
~~~
 sbatch run_baron10.sh
~~~
</details>



