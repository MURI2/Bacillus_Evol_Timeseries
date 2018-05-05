#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=4:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe

module load python
module load gcc
module load bwa/0.7.2
module load samtools/0.1.19

REF=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Plasmid_Genome.fa
R1=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/D100/Sample_L0B1-100/L0B1-100_GTAGAGGA-CTAAGCCT_L002_R1_001_clean_paired.fastq.gz
R2=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/D100/Sample_L0B1-100/L0B1-100_GTAGAGGA-CTAAGCCT_L002_R2_001_clean_paired.fastq.gz
Out=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test

#bwa index $REF
#samtools faidx $REF

#bwa mem -t 4 $REF $R1 $R2 > "${Out}_mapped.sam"
# mapped reads
#samtools view -F 4 -bT $REF "${Out}_mapped.sam" > "${Out}_mapped.bam"
# unmapped reads
#samtools view -f 4 -bT $REF "${Out}_mapped.sam" > "${Out}_unmapped.bam"

#samtools sort "${Out}_mapped.bam" "${Out}_mapped_sort"
#samtools index "${Out}_mapped_sort.bam"
#samtools rmdup "${Out}_mapped_sort.bam" "${Out}_mapped_sort_NOdup.bam"
#samtools index "${Out}_mapped_sort_NOdup.bam"
#samtools sort "${Out}_mapped_sort_NOdup.bam" "${Out}_mapped_sort_NOdup_sort"
#samtools index "${Out}_mapped_sort_NOdup_sort.bam"
#samtools view -h -o "${Out}_mapped_sort_NOdup_sort.sam" "${Out}_mapped_sort_NOdup_sort.bam"

# same thing for unmapped reads
#samtools sort "${Out}_unmapped.bam" "${Out}_unmapped_sort"
#samtools index "${Out}_unmapped_sort.bam"
#samtools rmdup "${Out}_unmapped_sort.bam" "${Out}_unmapped_sort_NOdup.bam"
#samtools index "${Out}_unmapped_sort_NOdup.bam"
#samtools sort "${Out}_unmapped_sort_NOdup.bam" "${Out}_unmapped_sort_NOdup_sort"
#samtools index "${Out}_unmapped_sort_NOdup_sort.bam"
#samtools view -h -o "${Out}_unmapped_sort_NOdup_sort.sam" "${Out}_unmapped_sort_NOdup_sort.bam"
# bam to fastq
#samtools bam2fq "${Out}_mapped_sort_NOdup_sort.bam" > "${Out}_mapped_sort_NOdup_sort.fastq"
#samtools bam2fq "${Out}_unmapped_sort_NOdup_sort.bam" > "${Out}_unmapped_sort_NOdup_sort.fastq"
# random sample of fastq
#fastq_map=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_mapped_sort_NOdup_sort.fastq
#fastq_map_sample=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_mapped_sort_NOdup_sort_sample.fastq
#fastq_unmap=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_unmapped_sort_NOdup_sort.fastq
#fastq_unmap_sample=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_unmapped_sort_NOdup_sort_sample.fastq

#cat /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_mapped_sort_NOdup_sort.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("\t");} }' \
#  | awk -v k=100 ' BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x); if(s<k)R[s]=$0}END{for(i in R)print R[i]}' \
#  | awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4 > "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_mapped_sort_NOdup_sort_sample.fastq"}'

#cat /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_unmapped_sort_NOdup_sort.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("\t");} }' \
#  | awk -v k=100 ' BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x); if(s<k)R[s]=$0}END{for(i in R)print R[i]}' \
#  | awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4 > "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_unmapped_sort_NOdup_sort_sample.fastq"}'
#convert fastq to fasta
#sed -n '1~4s/^@/>/p;2~4p' "${Out}_mapped_sort_NOdup_sort_sample.fastq" > "${Out}_mapped_sort_NOdup_sort_sample.fasta"
#sed -n '1~4s/^@/>/p;2~4p' "${Out}_unmapped_sort_NOdup_sort_sample.fastq" > "${Out}_unmapped_sort_NOdup_sort_sample.fasta"


# plasmid size = 84215
# genome size = 4215607
# ratio = 0.0200

# bam plasmid size = 7376238
# unmapped bam size = 362466002
# ratio = 0.0204

# average coverage = 248.909
#samtools depth "${Out}_mapped_sort_NOdup_sort.bam" |  awk '{sum+=$3} END { print "Average = ",sum/NR}'

REF1=/N/dc2/projects/muri2/Task2/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_16S_rRNA.fna

bwa index $REF1
samtools faidx $REF1

bwa mem -t 4 $REF1 $R1 $R2 > "${Out}_16S_mapped.sam"
# mapped reads
samtools view -F 4 -bT $REF1 "${Out}_16S_mapped.sam" > "${Out}_16S_mapped.bam"
# unmapped reads
samtools view -f 4 -bT $REF1 "${Out}_16S_mapped.sam" > "${Out}_16S_unmapped.bam"

samtools sort "${Out}_16S_mapped.bam" "${Out}_16S_mapped_sort"
samtools index "${Out}_16S_mapped_sort.bam"
samtools rmdup "${Out}_16S_mapped_sort.bam" "${Out}_16S_mapped_sort_NOdup.bam"
samtools index "${Out}_16S_mapped_sort_NOdup.bam"
samtools sort "${Out}_16S_mapped_sort_NOdup.bam" "${Out}_16S_mapped_sort_NOdup_sort"
samtools index "${Out}_16S_mapped_sort_NOdup_sort.bam"
samtools view -h -o "${Out}_16S_mapped_sort_NOdup_sort.sam" "${Out}_16S_mapped_sort_NOdup_sort.bam"

# same thing for unmapped reads
samtools sort "${Out}_16S_unmapped.bam" "${Out}_16S_unmapped_sort"
samtools index "${Out}_16S_unmapped_sort.bam"
samtools rmdup "${Out}_16S_unmapped_sort.bam" "${Out}_16S_unmapped_sort_NOdup.bam"
samtools index "${Out}_16S_unmapped_sort_NOdup.bam"
samtools sort "${Out}_16S_unmapped_sort_NOdup.bam" "${Out}_16S_unmapped_sort_NOdup_sort"
samtools index "${Out}_16S_unmapped_sort_NOdup_sort.bam"
samtools view -h -o "${Out}_16S_unmapped_sort_NOdup_sort.sam" "${Out}_16S_unmapped_sort_NOdup_sort.bam"
# bam to fastq
samtools bam2fq "${Out}_16S_mapped_sort_NOdup_sort.bam" > "${Out}_16S_mapped_sort_NOdup_sort.fastq"
samtools bam2fq "${Out}_16S_unmapped_sort_NOdup_sort.bam" > "${Out}_16S_unmapped_sort_NOdup_sort.fastq"
# random sample of fastq
fastq_map=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_16S_mapped_sort_NOdup_sort.fastq
fastq_map_sample=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_16S_mapped_sort_NOdup_sort_sample.fastq
fastq_unmap=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_16S_unmapped_sort_NOdup_sort.fastq
fastq_unmap_sample=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_16S_unmapped_sort_NOdup_sort_sample.fastq

cat /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_16S_mapped_sort_NOdup_sort.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("\t");} }' \
  | awk -v k=100 ' BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x); if(s<k)R[s]=$0}END{for(i in R)print R[i]}' \
  | awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4 > "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_16S_mapped_sort_NOdup_sort_sample.fastq"}'

cat /N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_16S_unmapped_sort_NOdup_sort.fastq | awk '{ printf("%s",$0); n++; if(n%4==0) {printf("\n");} else { printf("\t");} }' \
  | awk -v k=100 ' BEGIN{srand(systime() + PROCINFO["pid"]);}{s=x++<k?x-1:int(rand()*x); if(s<k)R[s]=$0}END{for(i in R)print R[i]}' \
  | awk -F"\t" '{print $1"\n"$2"\n"$3"\n"$4 > "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/bacillus_test/bacillus_test_16S_unmapped_sort_NOdup_sort_sample.fastq"}'
#convert fastq to fasta
sed -n '1~4s/^@/>/p;2~4p' "${Out}_16S_mapped_sort_NOdup_sort_sample.fastq" > "${Out}_16S_mapped_sort_NOdup_sort_sample.fasta"
sed -n '1~4s/^@/>/p;2~4p' "${Out}_16S_unmapped_sort_NOdup_sort_sample.fastq" > "${Out}_16S_unmapped_sort_NOdup_sort_sample.fasta"
