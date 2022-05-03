# Reproduce

This repository contains scripts and information to reproduce all results in the scAnnotate manuscript.

# About this project 

We build an automated cell annotation tool based on supervised machine learning
algorithms. We use a marginal mixture model to describe both the dropout proportion and the non-dropout expression level distribution of a gene. We develop a marginal model based ensemble learning approach to avoid having to specify and estimate a high-dimensional joint distribution for all genes. First, we build a ‘weak’ classifier using the mixture model for each gene. Then, we combine ‘weak’ classifiers of all genes into a single ‘strong’ classifier to annotate cells.

---

# Directory Layout

<details><summary>code</summary>

    ├── code                                               # All of the scripts to reproduce the results in the scAnnotate manuscript.
    │     ├── 1. PBMC
    │            ├── preprocessing_PBMC_crossplatform.sh   # Shell script to run preprocessing_PBMC_crossplatform.R
    │            ├── preprocessing_PBMC_crossplatform.R    # R script to preprocess the PBMC dataset
    │            ├── run_pbmc.sh                           # Shell script with loops to submit all shell scripts for each method on the PBMC dataset at once
    │            ├── run.sh                                # Shell script to run each corresponding R script on the selected dataset
    │            ├── scAnnotate.r                          # R script to run our method, scAnnotate, on the selected dataset
    │            ├── CaSTLe.r                              # R script to run competing method CaSTLe on the selected dataset
    │            ├── CHETAH.r                              # R script to run competing method CHETAH on the selected dataset
    │            ├── scClassify.r                          # R script to run competing method scClassify on the selected dataset
    │            ├── scID.r                                # R script to run competing method scID on the selected dataset
    │            ├── scmapCell.r                           # R script to run competing method scmapCell on the selected dataset
    │            ├── scmapCluster.r                        # R script to run competing method scmapCluster on the selected dataset
    │            ├── scPred.r                              # R script to run competing method scPred on the selected dataset
    │            ├── singleCellNet.r                       # R script to run competing method singleCellNet on the selected dataset
    │            └── SingleR.r                             # R script to run competing method SingleR on the selected dataset
    |
    │     ├── 2. CellBench
    │            ├── preprocessing_CellBench.sh            # Shell script to run preprocessing_CellBench.R
    │            ├── preprocessing_CellBench.R             # R script to preprocess the CellBench dataset
    │            ├── run_cellbench.sh                      # Shell script with loops to submit all shell scripts for each method on the CellBench dataset at once
    │            ├── run.sh                                # Shell script to run each corresponding R script on the selected dataset
    │            ├── scAnnotate.r                          # R script to run our method, scAnnotate, on the selected dataset
    │            ├── CaSTLe.r                              # R script to run competing method CaSTLe on the selected dataset
    │            ├── CHETAH.r                              # R script to run competing method CHETAH on the selected dataset
    │            ├── scClassify.r                          # R script to run competing method scClassify on the selected dataset
    │            ├── scID.r                                # R script to run competing method scID on the selected dataset
    │            ├── scmapCell.r                           # R script to run competing method scmapCell on the selected dataset
    │            ├── scmapCluster.r                        # R script to run competing method scmapCluster on the selected dataset
    │            ├── scPred.r                              # R script to run competing method scPred on the selected dataset
    │            ├── singleCellNet.r                       # R script to run competing method singleCellNet on the selected dataset
    │            └── SingleR.r                             # R script to run competing method SingleR on the selected dataset  
    |
    │     ├── 3. Brain
    │            ├── 3.1. ALM_MTG
    |                    ├── preprocessing_ALM_MTG.sh      # Shell script to run preprocessing_ALM_MTG.R
    │                    ├── preprocessing_ALM_MTG.R       # R script to preprocess the ALM and MTG datasets
    │                    ├── run_ALM_MTG.sh                # Shell script with loops to submit all shell scripts for each method on the ALM and MTG datasets at once
    │                    ├── run.sh                        # Shell script to run each corresponding R script on the selected dataset
    │                    ├── scAnnotate.r                  # R script to run our method, scAnnotate, on the selected dataset
    │                    ├── CaSTLe.r                      # R script to run competing method CaSTLe on the selected dataset
    │                    ├── CHETAH.r                      # R script to run competing method CHETAH on the selected dataset
    │                    ├── scClassify.r                  # R script to run competing method scClassify on the selected dataset
    │                    ├── scID.r                        # R script to run competing method scID on the selected dataset
    │                    ├── scmapCell.r                   # R script to run competing method scmapCell on the selected dataset
    │                    ├── scmapCluster.r                # R script to run competing method scmapCluster on the selected dataset
    │                    ├── scPred.r                      # R script to run competing method scPred on the selected dataset
    │                    ├── singleCellNet.r               # R script to run competing method singleCellNet on the selected dataset
    │                    └── SingleR.r                     # R script to run competing method SingleR on the selected dataset   
    |
    │            └── 3.2. V1_MTG   
    |                    ├── preprocessing_V1_MTG.sh      # Shell script to run preprocessing_V1_MTG.R
    │                    ├── preprocessing_V1_MTG.R       # R script to preprocess the V1 and MTG datasets
    │                    ├── run_V1_MTG.sh                # Shell script with loops to submit all shell scripts for each method on the V1 and MTG datasets at once
    │                    ├── run.sh                        # Shell script to run each corresponding R script on the selected dataset
    │                    ├── scAnnotate.r                  # R script to run our method, scAnnotate, on the selected dataset
    │                    ├── CaSTLe.r                      # R script to run competing method CaSTLe on the selected dataset
    │                    ├── CHETAH.r                      # R script to run competing method CHETAH on the selected dataset
    │                    ├── scClassify.r                  # R script to run competing method scClassify on the selected dataset
    │                    ├── scID.r                        # R script to run competing method scID on the selected dataset
    │                    ├── scmapCell.r                   # R script to run competing method scmapCell on the selected dataset
    │                    ├── scmapCluster.r                # R script to run competing method scmapCluster on the selected dataset
    │                    ├── scPred.r                      # R script to run competing method scPred on the selected dataset
    │                    ├── singleCellNet.r               # R script to run competing method singleCellNet on the selected dataset
    │                    └── SingleR.r                     # R script to run competing method SingleR on the selected dataset   
    |
    │     ├── 4. Pancreas
    │            ├── processing_Baron.sh                   # Shell script to run processing_Baron.R
    │            ├── processing_Baron.R                    # R script to preprocess the Baron pancreas dataset
    │            ├── run_panc.sh                           # Shell script with loops to submit all shell scripts for each method on the Baron pancreas dataset at once
    │            ├── run.sh                                # Shell script to run each corresponding R script on the selected dataset
    │            ├── scAnnotate.r                          # R scripts to run our method, scAnnotate, on the selected dataset
    │            ├── CaSTLe.r                              # R scripts to run competing method CaSTLe on the selected dataset
    │            ├── CHETAH.r                              # R scripts to run competing method CHETAH on the selected dataset
    │            ├── scClassify.r                          # R scripts to run competing method scClassify on the selected dataset
    │            ├── scID.r                                # R scripts to run competing method scID on the selected dataset
    │            ├── scmapCell.r                           # R scripts to run competing method scmapCell on the selected dataset
    │            ├── scmapCluster.r                        # R scripts to run competing method scmapCluster on the selected dataset
    │            ├── scPred.r                              # R scripts to run competing method scPred on the selected dataset
    │            ├── singleCellNet.r                       # R scripts to run competing method singleCellNet on the selected dataset
    │            └── SingleR.r                             # R scripts to run competing method SingleR on the selected dataset 
    │
    │     └── 5. CV_Baron
    │            ├── cv_Baronhuman.sh                      # Shell script to run cv_Baronhuman.R
    │            ├── cv_Baronhuman.R                       # R script to process the Baron human pancreas dataset
    │            ├── run_cellbench.sh                      # Shell script with loops to submit all shell scripts for each method on the Baron human pancreas dataset at once
    │            ├── run.sh                                # Shell script to run each corresponding R script on the selected dataset
    │            ├── scAnnotate.r                          # R script to run our method, scAnnotate, on the selected dataset
    │            ├── CaSTLe.r                              # R script to run competing method CaSTLe on the selected dataset
    │            ├── CHETAH.r                              # R script to run competing method CHETAH on the selected dataset
    │            ├── scClassify.r                          # R script to run competing method scClassify on the selected dataset
    │            ├── scID.r                                # R script to run competing method scID on the selected dataset
    │            ├── scmapCell.r                           # R script to run competing method scmapCell on the selected dataset
    │            ├── scmapCluster.r                        # R script to run competing method scmapCluster on the selected dataset
    │            ├── scPred.r                              # R script to run competing method scPred on the selected dataset
    │            ├── singleCellNet.r                       # R script to run competing method singleCellNet on the selected dataset
    │            └── SingleR.r                             # R script to run competing method SingleR on the selected dataset   
</details>

<details><summary>results</summary>
   
    ├── results                                            # All results used in the scAnnotate manuscript.
    │     ├── 1. annotation results
    │            ├── 1.1 PBMC                              # All annotation results for each method on the PBMC dataset
    │            ├── 1.2 CellBench                         # All annotation results for each method on the CellBench dataset 
    │            ├── 1.3 Brain                             # All annotation results for each method on the Brain dataset
    │            ├── 1.4 Pancreas                          # All annotation results for each method on the pancreas dataset
    │            └── 1.5 CV_Baron                          # All annotation results for each method on the Baron human dataset
    |
    │     └── 2. analysis scripts        
    │            ├── results_crossplatform.R               # R script to evaluate the performance of each method on cross-platform datasets by accuracy for each cell population, mean accuracy, and overall accuracy
    │            ├── results_crossspecies.R                # R script to evaluate the performance of each method on cross-species datasets by accuracy for each cell population, mean accuracy, and overall accuracy
    │            └──results_cv.R                          # R script to evaluate the performance of each method on intra-datasets (Baron_human) by accuracy for each cell population, mean accuracy, and overall accuracy 
</details>
  

---

# Package dependencies
- To reproduce the results, please download and install the required software and packages with the corresponding versions as listed in the table below.

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
- We have arranged the scripts by dataset type (1. PBMC 2. CellBench 3. Brain 4. Pancreas 5. CV_Baron); you can follow this order to reproduce all of our results. 
- All of our scripts were run on ComputeCanada. Before you submit our provided shell scripts, please change the account ID in each shell script (#SBATCH --account=) to that of your own ComputeCanada account. 
- All of the raw datasets used in our study are publicly available. Due to GitHub's file size limitation, we are unable to upload the raw data files here. Please download the raw datasets used in this study directly from their original sources before running the scripts.

<details><summary>1. PBMC </summary>
  
  1). Download the raw data: the PBMC dataset was downloaded from the SeuratData package with dataset name "pbmcsca". 
~~~
  R
  library(SeuratData)
  InstallData("pbmcsca")
~~~  
  2). Process the dataset: submit the shell script to process the PBMC dataset
~~~
 sbatch preprocessing_PBMC_crossplatform.sh
~~~
  3). Run the methods: submit the shell script to run each method's R script on the selected dataset.
~~~
 sbatch run_pbmc.sh
~~~
</details>

<details><summary>2. CellBench </summary> 
  
  1). Download the raw data: CellBench was download from Zenodo (https://zenodo.org/record/3357167#.YjpANOfMKw4)
  
  2). Process the dataset: submit the shell script to process the CellBench dataset
~~~
 sbatch preprocessing_CellBench.sh
~~~
  
  3). Run the methods: submit the shell script to run each method's R script on the selected dataset.
~~~
 sbatch run_cellbench.sh
~~~
</details>

<details><summary>3. Brain </summary>
  
  3.1 ALM_MTG
  
  1). Download the raw data: MouseALM_HumanMTG(**) was download from Zenodo (https://zenodo.org/record/3357167#.YjpANOfMKw4)
  
  2). Process the dataset: submit the shell script to process the ALM and MTG datasets
~~~
 sbatch preprocessing_ALM_MTG.sh
~~~
  
  3). Run the methods: submit the shell script to run each method's R script on the selected dataset.
~~~
 sbatch run_ALM_MTG.sh
~~~
  
  3.2 V1_MTG
 
  1). Download the raw data: MouseV1_HumanMTG(**) was download from Zenodo (https://zenodo.org/record/3357167#.YjpANOfMKw4)
  
  2). Process the dataset: submit the shell script to process the V1 and MTG datasets
~~~
 sbatch preprocessing_V1_MTG.sh
~~~
  
  3). Run the methods: submit the shell script to run each method's R script on the selected dataset.
~~~
 sbatch run_V1_MTG.sh
~~~

</details>


<details><summary>4. Pancreas </summary>
  
  1). Download the raw data: the pancreas data was downloaded from National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) for GSE84133. https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84133&format=file

  2). Process the dataset: submit the shell script to process the Baron pancreas dataset
~~~
 sbatch processing_Baron.sh
~~~
  3). Run the methods: submit the shell script to run each method's R script on the selected dataset.
~~~
 sbatch run_panc.sh
~~~
</details>


<details><summary>5. CV_Baron </summary>
  
  1). Move the selected dataset into the (**home directory?) folder: the selected dataset "Baron_human.rds" is outputted from processing_Baron.R (which is already run in Step 4 above)

  2). Process the dataset: submit the shell script to process the Baron human pancreas dataset
~~~
 sbatch cv_Baronhuman.sh
~~~
  3). Run the methods: submit the shell script to run each method's R script on the selected dataset.
~~~
 sbatch run_baron10.sh
~~~
</details>



