#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=24:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -m n
#PBS -j oe

# convert BAM files from the repbreseq step to an ordered SAM file
# this is for iRep/bPTR
module load samtools

declare -a strains=("B" "S")
declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")
declare -a times=("100" "200" "300")

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq_sam

for treat in "${treats[@]}"
do
  for strain in "${strains[@]}"
  do
    for rep in "${reps[@]}"
    do
      for time in "${times[@]}"
      do
        sample='Sample_L${treat}${strain}${rep}-${time}'
        bam=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq/Sample_L${treat}${strain}${rep}-${time}/data/reference.bam
        sam=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq_sam/Sample_L${treat}${strain}${rep}-${time}.sam
        if [ -f $bam ]; then
          samtools sort $bam | samtools view -h -o $sam
        fi
      done
    done
  done
done
