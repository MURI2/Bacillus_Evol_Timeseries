#!/bin/bash

declare -a strains=("B" "S")
declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")
declare -a times=("100" "200" "300")

annotated_path='/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq_annotated'
evidence_path='/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq_evidence'
output_path='/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq_output'
mkdir -p $annotated_path
mkdir -p $evidence_path
mkdir -p $output_path

for treat in "${treats[@]}"
do
  for strain in "${strains[@]}"
  do
    for rep in "${reps[@]}"
    do
      for time in "${times[@]}"
      do
        sample='Sample_L${treat}${strain}${rep}-${time}'
        annotated=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq/Sample_L${treat}${strain}${rep}-${time}/output/evidence/annotated.gd
        evidence=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq/Sample_L${treat}${strain}${rep}-${time}/output/evidence/evidence.gd
        output=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq/Sample_L${treat}${strain}${rep}-${time}/output/output.gd
        #echo ${annotated_path}/
        cp $annotated ${annotated_path}/Sample_L${treat}${strain}${rep}-${time}.gd
        cp $evidence ${evidence_path}/Sample_L${treat}${strain}${rep}-${time}.gd
        cp $output ${output_path}/Sample_L${treat}${strain}${rep}-${time}.gd
      done
    done
  done
done
