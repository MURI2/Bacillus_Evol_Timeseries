#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=14:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

# iRep and bPTR uses paired-end information, which breseq doesn't use
# so we'll need to map the reads using bwa and then clean them a little
# with samtools

module load bwa
module load samtools

ref=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.fna

bwa index $ref
samtools faidx $ref

mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_bam
#mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_sam_merged
mkdir -p /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_bam_merged

#declare -a strains=("B")
#declare -a treats=("0")
#declare -a reps=("2")
#declare -a times=("100")

declare -a strains=("B" "S")
declare -a treats=("0" "1" "2")
declare -a reps=("1" "2" "3" "4" "5")
declare -a times=("100" "200" "300")

for treat in "${treats[@]}"
do
  for strain in "${strains[@]}"
  do
    for rep in "${reps[@]}"
    do
      for time in "${times[@]}"
      do
        sample=Sample_L${treat}${strain}${rep}-${time}
        if [ ! -d /N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${sample} ]; then
          continue
        fi
        declare -a seq_reps=()
        for file in /N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${sample}/*R1_*_paired.fastq.gz
        do
          seq_rep="$(  echo "$file" | cut -d'.' -f 1 | rev | cut -d"_" -f3-4 | rev)"
          seq_reps=("${seq_reps[@]}" "$seq_rep")
        done
        for seq in "${seq_reps[@]}"
        do
          bwa_bam_out=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_bam/${sample}_${seq}.bam
          R1=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${sample}/*_R1_${seq}*
          R2=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/${sample}/*_R2_${seq}*
          bwa mem -t 4 $ref $R1 $R2 | samtools view -F 4 -bT $ref - \
            | samtools sort -o $bwa_bam_out
        done
        bwa_bam_merged_out=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_bam_merged/${sample}.bam
        samtools merge - /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bwa_bam/${sample}_*.bam \
        | samtools sort -o $bwa_bam_merged_out
      done
    done
  done
done
