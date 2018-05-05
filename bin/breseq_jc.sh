#!/bin/bash

module load breseq

#declare -a strains=("B" "C" "D" "F" "J" "P" "S")
#declare -a treats=("0" "1" "2")
#declare -a reps=("1" "2" "3" "4" "5")
declare -a strains=("B" "S")
declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")
declare -a times=("100" "200" "300")

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_jc
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_jc/all
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_jc/merged



for treat in "${treats[@]}"
do
  for strain in "${strains[@]}"
  do
    for rep in "${reps[@]}"
    do
      pop="Sample_L${treat}${strain}${rep}"
      for time in "${times[@]}"
      do
        evidence_time="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq/${pop}-${time}/output/evidence/evidence.gd"
        if [ -f $evidence_time ] ; then
          junction_output="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_jc/all/${pop}-${time}.gd"
          cat $evidence_time | grep 'JC\|#' > $junction_output
        fi
      done
      junction_merged_output="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_jc/merged/${pop}.gd"
      if [ -f $junction_merged_output ] ; then
        rm $junction_merged_output
      fi
      evidence="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_jc/all/${pop}"*".gd"
      gdtools UNION -o $junction_merged_output -e $evidence
    done
  done
done
