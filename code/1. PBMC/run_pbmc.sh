#!/bin/bash
#SBATCH --account=rrg-ubcxzh

for methods in "CaSTLe" "CHETAH" "scClassify" "scmapCell" "scmapCluster" "scPred" "singleCellNet" "SingleR" "scID" "scAnnotate"
do
	for train in "pbmc_10x-v3" "pbmc_10x-v2" "pbmc_Smart-seq2" "pbmc_Seq-Well" "pbmc_inDrops" "pbmc_Drop-seq" "pbmc_CEL-Seq2"
	do
		for test in "pbmc_10x-v3" "pbmc_10x-v2" "pbmc_Smart-seq2" "pbmc_Seq-Well" "pbmc_inDrops" "pbmc_Drop-seq" "pbmc_CEL-Seq2"	
		do
			if [ $train != $test ]
			then
				sbatch run.sh $methods $train $test
			fi
		done
	done
done



