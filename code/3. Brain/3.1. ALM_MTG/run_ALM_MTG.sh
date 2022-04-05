#!/bin/bash
#SBATCH --account=rrg-ubcxzh

for methods in "CaSTLe" "CHETAH" "scClassify" "scmapCell" "scmapCluster" "scPred" "singleCellNet" "SingleR" "scID" "scAnnotate"
do
	for train in "ALM3" "MTG3"
	do
		for test in "ALM3" "MTG3"
		do
			if [ $train != $test ]
			then
				sbatch run.sh $methods $train $test
			fi
		done
	done
done
	

