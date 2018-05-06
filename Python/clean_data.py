from __future__ import division
import os
import numpy as np
from Bio import SeqIO
import bacillus_tools as bt

mydir = os.path.expanduser("~/GitHub/Bacillus_Evol_Timeseries")


def clean_GBK():
    IN_path = mydir + '/data/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff'
    genome = SeqIO.parse(IN_path, "genbank")
    # protein_id
    df_out = open(mydir + '/data/gene_table.txt', 'w')
    header = ['LocusTag', 'protein_id' , 'Gene', 'Type', 'Size', 'GC', 'Sequence', 'Fold_1', \
            'Fold_2', 'Fold_2_S', 'Fold_2_V', 'Fold_3', 'Fold_4', 'N', 'S']
    df_out.write('\t'.join(header) + '\n')
    types_keep = ['CDS', 'rRNA', 'tRNA', 'tmRNA']
    total = []
    total1 = []
    for record in genome:
        if 'chromosome' in record.description:
            descript = record.description
            descript_split = descript.split(' ')
            descript_split_index = descript_split.index('chromosome')
            chrom = descript_split[descript_split_index] + '_' + descript_split[descript_split_index + 1].strip(',')
        elif 'plasmid' in record.description:
            descript = record.description
            descript_split = descript.split(' ')
            descript_split_index = descript_split.index('plasmid')
            chrom = descript_split[descript_split_index] + '_' + descript_split[descript_split_index + 1].strip(',')
        else:
            chrom = 'Genome'
        for f in record.features:
            total.append(f)
            if f.type not in types_keep:
                continue
            total1.append(f)
            if 'gene' in f.qualifiers:
                gene = f.qualifiers["gene"][0]
                gene = gene.replace(" ", "_")
            else:
                gene = 'nan'
            locus_tag = f.qualifiers["locus_tag"][0]
            if 'protein_id' in f.qualifiers:
                protein_id = f.qualifiers["protein_id"][0]
            else:
                protein_id = 'nan'
            size = f.location.end - f.location.start
            seq = str(f.extract(record.seq))
            if f.strand == -1:
                seq = seq[::-1]
            GC = round((seq.count('G') + seq.count('C')) / len(seq), 4)
            if f.type == "CDS":
                start_rf = int(f.qualifiers['codon_start'][0]) - 1
                codons = [seq[i + start_rf: i + start_rf + 3 ] for i in range(0, len(seq), 3)]
                nuc_list = ['A', 'C', 'G', 'T']
                codons = [x for x in codons if (len(x) == 3) and (len(np.setdiff1d(list(x),nuc_list)) == 0)]
                fold_1 = 0
                fold_2 = 0
                fold_3 = 0
                fold_4 = 0
                fold_2_V =0
                fold_2_S =0
                N = 0
                for codon in codons:
                    codon_list = list(codon)
                    N_codon = 0
                    for g in range(3):
                        fold_count = 0
                        fold_2_S_i = 0
                        fold_2_V_i = 0
                        for nuc in nuc_list:
                            codon_mut_list = list(codon_list)
                            if codon_mut_list[g] == nuc:
                                continue
                            codon_mut_list[g] = nuc
                            codon_mut = "".join(codon_mut_list)
                            S_V = bt.get_ts_tv_dict()[(codon_mut[g], codon[g])]
                            if bt.get_codon_dict()[codon_mut] == bt.get_codon_dict()[codon]:
                                fold_count += 1
                                if S_V == 'S':
                                    fold_2_S_i += 1
                                elif S_V == 'V':
                                    fold_2_V_i += 1
                        if fold_count == 0:
                            fold_1 += 1
                        elif fold_count == 1:
                            fold_2 += 1
                            if fold_2_S_i == 1 and fold_2_V_i == 0:
                                fold_2_S += 1
                            elif fold_2_S_i == 0 and fold_2_V_i == 1:
                                fold_2_V += 1
                            else:
                                #print(fold_S_count, fold_V_count)
                                continue
                        elif fold_count == 2:
                            fold_3 += 1
                        elif fold_count == 3:
                            fold_4 += 1
                        N_codon += (3 - fold_count) / 3
                    N += N_codon
                # synonymous sites.
                # calculated using http://bioinformatics.cvr.ac.uk/blog/calculating-dnds-for-ngs-datasets/
                S = (3*len(codons)) - N
                N = round(N, 2)
                S = round(S, 2)
                # fold_2_S and fold_2_V calculated using Comeron, 1995
                out_line = [locus_tag, protein_id, gene, f.type, size, GC, chrom, fold_1, fold_2, \
                        fold_2_S, fold_2_V, fold_3, fold_4, N, S]
                #out_line = [str(x) for x in out_line]

            else:
                out_line = [locus_tag, protein_id, gene, f.type, size, GC, chrom, 'nan', 'nan', \
                        'nan', 'nan', 'nan', 'nan', 'nan', 'nan']
                #out_line = [str(x) for x in out_line]
            print(locus_tag)
            df_out.write('\t'.join([str(x) for x in out_line]) + '\n')


    df_out.close()


def get_pop_by_gene_matrix():
    # just bother with day 100 for now
    gene_pop_matrix = {}
    to_keep = ['SNP', 'INS', 'DEL']
    directory = os.fsencode(mydir + '/data/pool_pop_seq/rebreseq_annotated')
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('-100.gd'):
            in_df = open(os.path.join(str(directory, 'utf-8'), filename), 'r')
            for line in in_df:
                line_split = line.strip().split()
                if line_split[0] not in to_keep:
                    continue
                #print(line_split)
                if line_split[8].split('=')[1] == 'intergenic':
                    continue
                # choose by locus_tag
                #print(line_split)
                #print(line_split[10])
                locus_tag_list = [s for s in line_split if 'locus_tag=' in s]

                #gene_name =
                #print(line_split)
                #gene_name = line_split[]
            #print(type(directory))
            #print(str(directory, 'utf-8'))
            #print(type())
            #print(bfilename.decode("utf-8") )
            #in_df =

#clean_GBK()
get_pop_by_gene_matrix()
