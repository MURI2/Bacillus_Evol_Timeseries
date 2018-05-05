#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=15:00:00
#PBS -M wrshoema@iu.edu
#PBS -m abe
#PBS -m n
#PBS -j oe

module load samtools
module load python

#declare -a strains=("B" "C" "D" "F" "J" "P" "S")
#declare -a treats=("0" "1" "2")
#declare -a reps=("1" "2" "3" "4" "5")

declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/pileup
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/timecourse

#ref=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff
#ref=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.fna
ref=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq/Sample_L0B1-100/data/reference.fasta
samtools faidx $ref

for treat in "${treats[@]}"
do
  for rep in "${reps[@]}"
  do
    pop="Sample_L${treat}B${rep}"
    bams="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq/${pop}"*"/data/reference.bam"
    out="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/pileup/${pop}.pileup"
    samtools mpileup -q10 -f ${ref} ${bams} > ${out}
    out_timecourse="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/timecourse/${pop}_timecourse.txt"
    cat ${out} | python /N/dc2/projects/muri2/Task2/PoolPopSeq/bin/create_timecourse.py ${pop} > ${out_timecourse}
    rm ${out}
  done
done
