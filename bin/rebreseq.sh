#!/bin/bash

bin_rebreseq_scripts="/N/dc2/projects/muri2/Task2/PoolPopSeq/bin/rebreseq_scripts"
data_rebreseq="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq"
data_rebreseq_out="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq_out"
data_rebreseq_err="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/rebreseq_err"

mkdir -p $bin_rebreseq_scripts
mkdir -p $data_rebreseq
mkdir -p $data_rebreseq_out
mkdir -p $data_rebreseq_err

ref=/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reference_assemblies_task2/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff

for line_path in "/N/dc2/projects/muri2/Task2/PoolPopSeq/data/reads_clean_trimmomatic/"*;
do
  line="$(echo "$line_path" | cut -d "/" -f10-10)"
  pop="$(echo "$line_path" | cut -d "/" -f10-10 | cut -d "-" -f1-1)"
  pop_gd="/N/dc2/projects/muri2/Task2/PoolPopSeq/data/breseq_jc/merged/${pop}.gd"
  bash_out="${bin_rebreseq_scripts}/${line}_rebreseq.sh"
  taxon="$(echo "$line_path" | grep -Po ".(?=.{5}$)")"
  if [ -f $bash_out ]; then
    rm $bash_out
  fi
  taxon="$(echo "$line_path" | grep -Po ".(?=.{5}$)")"
  if [[ $taxon == "B" ]] || [[ $taxon == "S" ]]; then
    reads="${line_path}/"*_clean_paired.fastq.gz
    OUT_rebreseq="${data_rebreseq}/${line}"
    mkdir -p $OUT_rebreseq
    OUT_rebreseq_out="${data_rebreseq_out}/${line}.out"
    if [ -f $OUT_rebreseq_out ]; then
      rm $OUT_rebreseq_out
    fi
    OUT_rebreseq_err="${data_rebreseq_err}/${line}.err"
    if [ -f $OUT_rebreseq_err ]; then
      rm $OUT_rebreseq_err
    fi

    echo '#!/bin/bash' >> $bash_out
    echo '#PBS -k o' >> $bash_out
    echo '#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=10:00:00' >> $bash_out
    echo '#PBS -M wrshoema@iu.edu' >> $bash_out
    echo '#PBS -m abe' >> $bash_out
    echo '#PBS -j oe' >> $bash_out
    echo '' >> $bash_out
    echo 'module load breseq' >> $bash_out
    # > data/breseq_output/${sample_name}.out 2> data/breseq_output/${sample_name}.err
    echo "breseq -j 8 -p --user-evidence-gd ${pop_gd} -o ${OUT_rebreseq} -r ${ref} ${reads} > ${OUT_rebreseq_out} 2> ${OUT_rebreseq_err}" >> $bash_out
    echo '' >> $bash_out

    qsub $bash_out
  else
    continue
  fi
done
