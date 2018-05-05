#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=15:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -m n
#PBS -j oe

module load python

#declare -a taxa=("B")
declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")
#declare -a treats=("0")
#declare -a reps=("1")

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/timecourse_merged
ref=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq/Sample_L0B1-100/data/reference.fasta

for treat in "${treats[@]}"
do
  for rep in "${reps[@]}"
  do
    pop="Sample_L${treat}B${rep}"
    out_timecourse="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/timecourse/${pop}_timecourse.txt"
    out_timecourse_merged="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/timecourse_merged/${pop}_merged_timecourse.bz2"
    gd_files="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq/${pop}-"*"/output/evidence/evidence.gd"
    cat ${out_timecourse} | python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/create_breseq_timecourse_no_repeat.py ${ref} ${pop} ${gd_files} | bzip2 -c > ${out_timecourse_merged}
  done
done
