#!/bin/bash
#SBATCH --account=def-ubcxzh

for methods in "CaSTLe" "CHETAH" "scClassify" "scmapCell" "scmapCluster" "scPred" "singleCellNet" "SingleR" "scID" "scAnnotate"
do
	for index in {1..10}
	do
		train="train_"$index
		test="test_"$index
		
		sbatch run.sh $methods $train $test
		
	done
	
done



