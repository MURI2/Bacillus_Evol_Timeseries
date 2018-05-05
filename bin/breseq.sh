#!/bin/bash

bin_breseq_scripts="/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/breseq_scripts"
data_breseq="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq"
data_breseq_out="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_out"
data_breseq_err="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_err"

mkdir -p $bin_breseq_scripts
mkdir -p $data_breseq
mkdir -p $data_breseq_out
mkdir -p $data_breseq_err

ref=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff

for line_path in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/"*;
do
  line="$(echo "$line_path" | cut -d "/" -f10-10)"
  bash_out="${bin_breseq_scripts}/${line}_breseq.sh"
  if [ -f $bash_out ]; then
    rm $bash_out
  fi
  taxon="$(echo "$line_path" | grep -Po ".(?=.{5}$)")"
  if [[ $taxon == "B" ]] || [[ $taxon == "S" ]]; then

      reads="${line_path}/"*_clean_paired.fastq.gz
      OUT_breseq="${data_breseq}/${line}"
      mkdir -p $OUT_breseq
      OUT_breseq_out="${data_breseq_out}/${line}.out"
      if [ -f $OUT_breseq_out ]; then
        rm $OUT_breseq_out
      fi
      OUT_breseq_err="${data_breseq_err}/${line}.err"
      if [ -f $OUT_breseq_err ]; then
        rm $OUT_breseq_err
      fi

      echo '#!/bin/bash' >> $bash_out
      echo '#PBS -k o' >> $bash_out
      echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=30:00:00' >> $bash_out
      echo '#PBS -M wrshoema@iu.edu' >> $bash_out
      echo '#PBS -m abe' >> $bash_out
      echo '#PBS -j oe' >> $bash_out
      echo '' >> $bash_out
      echo 'module load breseq' >> $bash_out
      echo "breseq -j 8 -p --brief-html-output --polymorphism-reject-indel-homopolymer-length 0 --polymorphism-reject-surrounding-homopolymer-length 0 --polymorphism-score-cutoff 2 -o ${OUT_breseq} -r ${ref} ${reads} > ${OUT_breseq_out} 2> ${OUT_breseq_err}" >> $bash_out
      echo '' >> $bash_out

      #qsub $bash_out
  else
    continue
  fi
done
