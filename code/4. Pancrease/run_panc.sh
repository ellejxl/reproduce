#!/bin/bash
#SBATCH --account=rrg-ubcxzh

for methods in "CaSTLe" "CHETAH" "scClassify" "scmapCell" "scmapCluster" "scPred" "singleCellNet" "SingleR" "scID" "scAnnotate"
do
	for train in "panc_mouse" "panc_human"
	do
		for test in "panc_mouse" "panc_human"

		do
			if [ $train != $test ]
			then
				sbatch run.sh $methods $train $test
			fi
		done
	done
		
done




