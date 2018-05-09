from __future__ import division
import os, subprocess, re
import numpy as np
import pandas as pd
from Bio import SeqIO
import bacillus_tools as bt


####
# add line to check whethe a gene is annotated as spore-related or not
####

def clean_GBK():
    IN_path = bt.get_path() + '/data/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff'
    genome = SeqIO.parse(IN_path, "genbank")
    # protein_id
    df_out = open(bt.get_path() + '/data/gene_table.txt', 'w')
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
            #print(f)
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
            else:
                #m = G + C -> A + T / A + T -> G + C
                out_line = [locus_tag, protein_id, gene, f.type, size, GC, chrom, 'nan', 'nan', \
                        'nan', 'nan', 'nan', 'nan', 'nan', 'nan']
            #print(locus_tag)
            df_out.write('\t'.join([str(x) for x in out_line]) + '\n')
    df_out.close()


def get_iRep():
    directory = os.fsencode(bt.get_path() + '/data/bwa_sam_merged')
    fasta = bt.get_path() + '/data/Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.fna'
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('-100.sam'):
            print(filename)
            out_file = bt.get_path() + '/data/bPTR/' + filename.split('.')[0]
            sam = os.path.join(str(directory, 'utf-8'), filename)
            #plot = bt.get_path() + '/data/bPTR/' + filename.split('.')[0]
            #subprocess.call(['irep', '-f', fasta, '-s', sam, '-o', out_file, '--sort'])
            subprocess.call(['bPTR', '-m', 'gc_skew', '-f', fasta, '-s', sam, '-o', str(out_file) + '.txt', '--sort', '-plot', str(out_file)])



def get_pop_by_gene_matrix():
    # just bother with day 100 for now
    to_exclude = bt.mutations_to_exclude()
    gene_pop_matrix = {}
    directory = os.fsencode(bt.get_path() + '/data/pool_pop_seq/rebreseq_annotated')
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        if filename.endswith('-100.gd'):
            in_df = open(os.path.join(str(directory, 'utf-8'), filename), 'r')
            pop = filename.split('.')[0]
            for line in in_df:
                line_split = line.strip().split()
                if (line_split[0] not in bt.get_to_keep()) or \
                    (line_split[8].split('=')[1] == 'intergenic') or \
                    (line_split[3] + '_' + line_split[4] in to_exclude):
                    continue
                # clean locus tags
                locus_tag = [s for s in line_split if 'locus_tag=' in s][0].split('=')[1]
                locus_tag_clean = re.sub('[][]', '', locus_tag)#.split(';')
                locus_tag_clean_split = re.findall(r"[\w']+", locus_tag_clean)
                for locus in locus_tag_clean_split:
                    if locus in gene_pop_matrix:
                        if pop in gene_pop_matrix[locus]:
                            gene_pop_matrix[locus][pop] += 1
                        else:
                            gene_pop_matrix[locus][pop] = 1
                    else:
                        gene_pop_matrix[locus] = {}
                        gene_pop_matrix[locus][pop] = 1

    df = pd.DataFrame.from_dict(gene_pop_matrix)
    df = df.fillna(0)
    df_out = bt.get_path() + '/data/pool_pop_seq/gene_by_pop.txt'
    df.to_csv(df_out, sep = '\t', index = True)





#clean_GBK()
#get_pop_by_gene_matrix()
get_iRep()
