from __future__ import division
import os, math, numbers, itertools, re
import pandas as pd
import numpy as np
import warnings
from Bio import SeqIO
warnings.filterwarnings('ignore')

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/data/")

# translation table 11
codon_dict = {
    "TTT":"F", "TCT":"S", "TAT":"Y", "TGT":"C",
    "TTC":"F", "TCC":"S", "TAC":"Y", "TGC":"C",
    "TTA":"L", "TCA":"S", "TAA":"*", "TGA":"*",
    "TTG":"L", "TCG":"S", "TAG":"*", "TGG":"W",

    "CTT":"L", "CCT":"P", "CAT":"H", "CGT":"R",
    "CTC":"L", "CCC":"P", "CAC":"H", "CGC":"R",
    "CTA":"L", "CCA":"P", "CAA":"Q", "CGA":"R",
    "CTG":"L", "CCG":"P", "CAG":"Q", "CGG":"R",

    "ATT":"I", "ACT":"T", "AAT":"N", "AGT":"S",
    "ATC":"I", "ACC":"T", "AAC":"N", "AGC":"S",
    "ATA":"I", "ACA":"T", "AAA":"K", "AGA":"R",
    "ATG":"M", "ACG":"T", "AAG":"K", "AGG":"R",

    "GTT":"V", "GCT":"A", "GAT":"D", "GGT":"G",
    "GTC":"V", "GCC":"A", "GAC":"D", "GGC":"G",
    "GTA":"V", "GCA":"A", "GAA":"E", "GGA":"G",
    "GTG":"V", "GCG":"A", "GAG":"E", "GGG":"G"
    }

ts_tv_dict = {
    ("A", "C"):"V", ("A", "G"):"S", ("A", "T"):"V",
    ("C", "A"):"V", ("C", "G"):"V", ("C", "T"):"S",
    ("G", "A"):"S", ("G", "C"):"V", ("G", "T"):"V",
    ("T", "A"):"V", ("T", "C"):"S", ("T", "G"):"V",
    }

merged_on = ['seq_id', 'position', 'gene_list', \
    'gene_name', 'gene_position', 'gene_product', \
    'locus_tag', 'gene_strand', 'transl_table', 'reference']

to_rename = ['mutation', 'codon_ref_seq', 'codon_new_seq', 'codon_number', \
    'codon_position', 'aa_position', 'aa_ref_seq', 'aa_new_seq', \
    'frequency', 'total_cov', 'number', 'file_number', \
    'prediction', 'consensus_score', 'polymorphism_score', \
    'fisher_strand_p_value', 'ks_quality_p_value', 'bias_e_value', \
    'bias_p_value', 'reject', 'snp_type', 'type', \
    'major_base', 'minor_base', 'sample']

def cleanGBK(strain):
    if strain == 'B':
        #file_name = 'Bacillus_subtilis_168/GCA_000009045.1_ASM904v1_genomic.gbff'
        file_name = 'Bacillus_subtilis_NCIB_3610/GCA_002055965.1_ASM205596v1_genomic.gbff'
    elif strain == 'C':
        file_name = 'Caulobacter_crescentus_NA1000/GCA_000022005.1_ASM2200v1_genomic.gbff'
    elif strain == 'D':
        file_name = 'Deinococcus_radiodurans_BAA816/GCA_000008565.1_ASM856v1_genomic.gbff'
    elif strain == 'F':
        file_name = 'Pedobacter_sp_KBS0701/G-Chr1.gbk'
    elif strain == 'J':
        file_name = 'Janthinobacterium_sp_KBS0711/KBS0711_2015_SoilGenomes_Annotate/G-Chr1.gbk'
    elif strain == 'P':
        file_name = 'Pseudomonas_sp_KBS0710/G-Chr1.gbk'
    else:
        print("Strain not recognized")
    IN_path = mydir + 'reference_assemblies_task2/' + file_name
    genome = SeqIO.parse(IN_path, "genbank")
    # protein_id
    OUT = open(mydir + 'reference_assemblies_task2/reference_assemblies_task2_table/' + strain + '.txt', 'w')
    print>> OUT, 'LocusTag', 'protein_id' , 'Gene', 'Type', 'Size', 'GC', 'Sequence', 'Fold_1', \
            'Fold_2', 'Fold_2_S', 'Fold_2_V', 'Fold_3', 'Fold_4', 'N', 'S'
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
        #print len(record.features)
        for f in record.features:
            total.append(f)
            if f.type not in types_keep:
                continue
            total1.append(f)
            #if f.type == "CDS" :
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
                            S_V = ts_tv_dict[(codon_mut[g], codon[g])]
                            if codon_dict[codon_mut] == codon_dict[codon]:
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
                                print fold_S_count, fold_V_count
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
                print>> OUT, locus_tag, protein_id, gene, f.type, size, GC, chrom, fold_1, fold_2, \
                        fold_2_S, fold_2_V, fold_3, fold_4, N, S
            else:
                print>> OUT, locus_tag, protein_id, gene, f.type, size, GC, chrom, 'nan', 'nan', \
                        'nan', 'nan', 'nan', 'nan', 'nan', 'nan'

    #print set(types)
    print len(total)
    print len(total1)
    OUT.close()



class cleanBreseq_evidence:

    def __init__(self, path):
        self.path = path

    def RA_line(self, row, columns):
        columns = columns[8:]
        row_0_8 = row[:8]
        row_8_inf = row[8:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row

    def MC_line(self, row, columns):
        columns = columns[8:]
        row_0_8 = row[:8]
        row_8_inf = row[8:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row

    def JC_line(self, row, columns):
        columns = columns[10:]
        row_0_8 = row[:10]
        row_8_inf = row[10:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row


    def split_evidence(self):
        with open(self.path) as f:
            for x ,line in enumerate(f):
                line = line.split()
        path_split = self.path.split('/')
        path = '/'.join(path_split[:-4]) + '/breseq_output_gbk_essentials_split/' + '/'.join(path_split[8:10])
        try:
            os.stat(path)
        except:
            os.mkdir(path)
        OUT_RA = open(path + '/evidence_RA.txt', 'w')
        OUT_MC = open(path + '/evidence_MC.txt', 'w')
        OUT_JC = open(path + '/evidence_JC.txt', 'w')
        OUT_UN = open(path + '/evidence_UN.txt', 'w')
        columns_RA = ['type','number', 'misc', 'seq_id', 'position', 'change', \
            'reference', 'sample', 'total_cov', 'new_cov', 'ref_cov', 'major_cov', \
            'minor_cov', 'major_base', 'minor_base', 'prediction', 'frequency', \
            'polymorphism_frequency', 'major_frequency', 'consensus_score', \
            'polymorphism_score', 'fisher_strand_p_value', 'ks_quality_p_value', \
            'bias_e_value', 'reject']

        columns_MC = ['type', 'number', 'misc', 'seq_id', 'start', 'finish', \
                'zero1', 'zero2', 'left_inside_cov', 'left_outside_cov', \
                'right_inside_cov', 'right_outside_cov']

        columns_JC = ['type', 'number', 'misc', 'seq_id_start', 'start', 'strand'\
            'seq_id_finish', 'finish', 'unknown1', 'unknown2' 'side_1_read_count', \
            'max_right_plus', 'coverage_minus', 'prediction', 'frequency', \
            'polymorphism_frequency', 'max_min_left_plus', 'side_2_annotate_key', \
            'alignment_overlap', 'max_left_plus', 'max_min_left', 'flanking_left', \
            'max_left', 'max_min_right_plus', 'total_non_overlap_reads', \
            'coverage_plus', 'side_2_overlap', 'neg_log10_pos_hash_p_value', \
            'max_left_minus', 'side_2_redundant', 'side_1_possible_overlap_registers', \
            'side_2_coverage', 'side_1_continuation', 'max_right_minus', \
            'side_1_redundant', 'side_2_possible_overlap_registers', \
            'unique_read_sequence', 'max_min_left_minus', \
            'new_junction_read_count', 'max_pos_hash_score', 'key', \
            'side_2_continuation', 'pos_hash_score', \
            'junction_possible_overlap_registers', 'side_1_overlap', \
            'max_min_right', 'flanking_right', 'max_right', 'side_1_coverage', \
            'side_2_read_count', 'max_min_right_minus', 'new_junction_coverage', \
            'side_1_annotate_key']

        columns_UN = ['type', 'number', 'misc', 'seq_id', 'start', 'finish']

        print>> OUT_RA, '\t'.join(columns_RA)
        print>> OUT_MC, '\t'.join(columns_MC)
        print>> OUT_JC, '\t'.join(columns_JC)
        print>> OUT_UN, '\t'.join(columns_UN)

        with open(self.path) as f:
            for line in f:
                line = line.split()
                if line[0] == 'RA':
                    line = self.RA_line(line, columns_RA)
                    print>> OUT_RA, '\t'.join(line)
                elif line[0] == 'MC':
                    line = self.MC_line(line, columns_MC)
                    print>> OUT_MC, '\t'.join(line)
                elif line[0] == 'JC':
                    line = self.JC_line(line, columns_JC)
                    print>> OUT_JC, '\t'.join(line)
                elif line[0] == 'UN':
                    print>> OUT_UN, '\t'.join(line)
                else:
                    continue

        OUT_RA.close()
        OUT_MC.close()
        OUT_JC.close()
        OUT_UN.close()


    def clean_evidence(self):

        IN = pd.read_csv(self.path, sep = '\t')


class cleanBreseq_annotated:

    def __init__(self, path):
        self.path = path

    def clean_value(self, row, value_name):
        value_name = value_name + '='
        data_indices = [i for i, s in enumerate(row) if '=' in s]
        gene_name_index = [i for i, s in enumerate(row) if value_name in s]
        if len(gene_name_index) > 0:
            gene_name_index = gene_name_index[0]
            gene_name_index_end = data_indices.index(gene_name_index)
            if gene_name_index_end+1 >= len(data_indices):
                gene_name = row[gene_name_index:]
                gene_name = '_'.join(gene_name)
                gene_name =  re.sub(r'[^\x00-\x7F]+','-', gene_name)
                return row[:gene_name_index] +  [gene_name]

            else:
                gene_name = row[gene_name_index:data_indices[gene_name_index_end+1]]
                gene_name = '_'.join(gene_name)
                gene_name =  re.sub(r'[^\x00-\x7F]+','-', gene_name)
                return row[:gene_name_index] +  [gene_name] + row[data_indices[gene_name_index_end+1]:]
        else:
            return row

    def variant_line(self, row, columns):
        row = self.clean_value(row, 'gene_product')
        row = self.clean_value(row, 'gene_position')
        row = self.clean_value(row, 'gene_name')
        row = self.clean_value(row, 'locus_tag')
        row = self.clean_value(row, 'between')
        if row[0] == 'SNP':
            # 7 b/c added column for insertion length (which is nan for SNPs)
            columns_6_inf = columns[7:]
            row_0_6 = row[:6]
            row_6_inf = row[6:]
            row_list = []
            values_dict = dict(item.split("=") for item in row_6_inf if '=' in item)
            for column in columns_6_inf:
                if column in values_dict:
                    row_list.append(values_dict[column])
                else:
                    row_list.append(float('nan'))
            row_list = [str(x) for x in row_list]
            # add one for point mutations
            new_row = row_0_6 + ['nan'] + row_list

            return new_row
        elif row[0] == 'SUB':
            row[6], row[5] = row[5], row[6]

            columns_7_inf = columns[7:]
            row_0_7 = row[:7]
            row_7_inf = row[7:]
            row_list = []
            values_dict = dict(item.split("=") for item in row_7_inf)
            for column in columns_7_inf:
                if column in values_dict:
                    row_list.append(values_dict[column])
                else:
                    row_list.append(float('nan'))
            row_list = [str(x) for x in row_list]
            new_row = row_0_7 + row_list
            return new_row

        elif row[0] == 'INS':
            # 7 b/c added column for insertion length (which is nan for SNPs)
            columns_6_inf = columns[7:]
            row_0_6 = row[:6]
            row_6_inf = row[6:]
            row_list = []
            values_dict = dict(item.split("=") for item in row_6_inf)
            for column in columns_6_inf:
                if column in values_dict:
                    row_list.append(values_dict[column])
                else:
                    row_list.append(float('nan'))
            row_list = [str(x) for x in row_list]
            # add length of insertion
            new_row = row_0_6 + [str(len(row_0_6[-1]))] + row_list
            return new_row

        elif row[0] == 'DEL':
            # 7 b/c added column for insertion length (which is nan for SNPs)
            columns_6_inf = columns[7:]
            row_0_6 = row[:6]
            row_6_inf = row[6:]
            row_list = []
            values_dict = dict(item.split("=") for item in row_6_inf if '=' in item)
            for column in columns_6_inf:
                if column in values_dict:
                    row_list.append(values_dict[column])
                else:
                    row_list.append(float('nan'))
            row_list = [str(x) for x in row_list]
            # add length of insertion
            new_row = row_0_6 + ['nan'] + row_list
            return new_row

    def RA_line(self, row, columns):
        row = self.clean_value(row, 'gene_product')
        row = self.clean_value(row, 'gene_position')
        row = self.clean_value(row, 'gene_name')
        row = self.clean_value(row, 'locus_tag')
        row = self.clean_value(row, 'between')
        columns = columns[8:]
        row_0_8 = row[:8]
        row_8_inf = row[8:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf if '=' in item)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row

    def MC_line(self, row, columns):
        row = self.clean_value(row, 'gene_product')
        row = self.clean_value(row, 'gene_position')
        row = self.clean_value(row, 'locus_tag')
        row = self.clean_value(row, 'gene_name')
        row = self.clean_value(row, 'gene_list')
        row = self.clean_value(row, 'html_gene_name')

        columns = columns[8:]
        row_0_8 = row[:8]
        row_8_inf = row[8:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row

    def JC_line(self, row, columns):
        columns = columns[10:]
        row_0_8 = row[:10]
        row_8_inf = row[10:]
        row_list = []
        values_dict = dict(item.split("=") for item in row_8_inf)
        for column in columns:
            if column in values_dict:
                row_list.append(values_dict[column])
            else:
                row_list.append(float('nan'))
        row_list = [str(x) for x in row_list]

        new_row = row_0_8 + row_list
        return new_row


    def split_annotated(self):
        path_split = self.path.split('/')
        path = '/'.join(path_split[:-4]) + '/breseq_output_gbk_essentials_split/' + '/'.join(path_split[8:10])
        try:
            os.stat(path)
        except:
            os.mkdir(path)
        OUT_RA = open(path + '/evidence_RA.txt', 'w')
        OUT_MC = open(path + '/evidence_MC.txt', 'w')
        OUT_JC = open(path + '/evidence_JC.txt', 'w')
        OUT_UN = open(path + '/evidence_UN.txt', 'w')
        OUT_variants = open(path + '/evidence_variants.txt', 'w')

        columns_variants = ['type','number', 'file_number', 'seq_id', 'position', \
            'mutation', 'size', 'frequency', 'gene_list', 'gene_name', 'gene_position', \
            'gene_product', 'locus_tag', 'snp_type', 'aa_new_seq', 'aa_position', \
            'aa_ref_seq', 'codon_new_seq', 'codon_number', 'codon_position', \
            'codon_ref_seq', 'gene_strand', 'transl_table', 'insert_position', 'between']

        columns_RA = ['type','number', 'misc', 'seq_id', 'position', 'change', \
            'reference', 'sample', 'total_cov', 'new_cov', 'ref_cov', 'major_cov', \
            'minor_cov', 'major_base', 'minor_base', 'prediction', 'frequency', \
            'polymorphism_frequency', 'major_frequency', 'consensus_score', \
            'polymorphism_score', 'fisher_strand_p_value', 'ks_quality_p_value', \
            'bias_e_value', 'reject', 'snp_type', 'locus_tag', 'gene_product', \
            'gene_position', 'gene_list', 'gene_name', 'bias_p_value']

        columns_MC = ['type', 'number', 'misc', 'seq_id', 'start', 'finish', \
                'zero1', 'zero2', 'left_inside_cov', 'left_outside_cov', \
                'right_inside_cov', 'right_outside_cov', 'gene_list', 'gene_name', \
                'gene_position', 'gene_product', 'locus_tag']

        columns_JC = ['type', 'number', 'misc', 'seq_id_start', 'start', 'strand'\
            'seq_id_finish', 'finish', 'unknown1', 'unknown2' 'side_1_read_count', \
            'max_right_plus', 'coverage_minus', 'prediction', 'frequency', \
            'polymorphism_frequency', 'max_min_left_plus', 'side_2_annotate_key', \
            'alignment_overlap', 'max_left_plus', 'max_min_left', 'flanking_left', \
            'max_left', 'max_min_right_plus', 'total_non_overlap_reads', \
            'coverage_plus', 'side_2_overlap', 'neg_log10_pos_hash_p_value', \
            'max_left_minus', 'side_2_redundant', 'side_1_possible_overlap_registers', \
            'side_2_coverage', 'side_1_continuation', 'max_right_minus', \
            'side_1_redundant', 'side_2_possible_overlap_registers', \
            'unique_read_sequence', 'max_min_left_minus', \
            'new_junction_read_count', 'max_pos_hash_score', 'key', \
            'side_2_continuation', 'pos_hash_score', \
            'junction_possible_overlap_registers', 'side_1_overlap', \
            'max_min_right', 'flanking_right', 'max_right', 'side_1_coverage', \
            'side_2_read_count', 'max_min_right_minus', 'new_junction_coverage', \
            'side_1_annotate_key']

        columns_UN = ['type', 'number', 'misc', 'seq_id', 'start', 'finish']

        print>> OUT_RA, '\t'.join(columns_RA)
        print>> OUT_MC, '\t'.join(columns_MC)
        print>> OUT_JC, '\t'.join(columns_JC)
        print>> OUT_UN, '\t'.join(columns_UN)
        print>> OUT_variants, '\t'.join(columns_variants)


        set_test = []
        with open(self.path) as f:
            for row in f:
                row = row.split()
                if len(row) < 3:
                    continue
                if row[0] == 'DEL' or row[0] == 'SNP' or \
                    row[0] == 'SUB' or row[0] == 'INS':
                    row_clean = self.variant_line(row, columns_variants)
                    print>> OUT_variants, '\t'.join(row_clean)
                elif row[0] == 'RA':
                    row_clean = self.RA_line(row, columns_RA)
                    print>> OUT_RA, '\t'.join(row_clean)
                elif row[0] == 'MC':
                    row_clean = self.MC_line(row, columns_RA)
                elif row[0] == 'JC':
                    row_clean = self.JC_line(row, columns_RA)
                elif row[0] == 'UN':
                    print>> OUT_UN, '\t'.join(row)
                else:
                    continue

        OUT_RA.close()
        OUT_MC.close()
        OUT_JC.close()
        OUT_UN.close()
        OUT_variants.close()

def get_variant_annotated(day, strain, treatments, reps, variant_type):
    for treatment in treatments:
        for rep in reps:
            in_path =  mydir + 'breseq_output_gbk_essentials_split/' + \
                day +'/Sample_L' +  treatment + strain + rep + '-' + day[1:]
            variants_path = in_path + '/evidence_variants.txt'
            RA_path = in_path + '/evidence_RA.txt'
            if os.path.exists(variants_path) == True:
                IN_variants = pd.read_csv(variants_path, sep = '\t', header = 'infer')
                IN_RA = pd.read_csv(RA_path, sep = '\t', header = 'infer')
                IN_variants_SNPs = IN_variants.loc[IN_variants['type'] == variant_type]
                IN_variants_SNPs = IN_variants_SNPs.drop(['size'], axis=1)
                drop_RA = ['gene_product', 'frequency', 'gene_list', 'gene_name', \
                        'gene_position', 'locus_tag', 'snp_type', 'type', 'number']
                IN_RA = IN_RA.drop(drop_RA, axis=1)
                IN_variants_SNPs.position = IN_variants_SNPs.position.astype(str)
                IN_variants_SNPs.seq_id = IN_variants_SNPs.seq_id.astype(str)

                IN_RA.position = IN_RA.position.astype(str)
                IN_RA.seq_id = IN_RA.seq_id.astype(str)
                IN_merged = pd.merge(IN_variants_SNPs, IN_RA, how='inner',
                    on=['seq_id', 'position'])
                out_path = mydir + 'breseq_output_gbk_essentials_split_clean/' \
                        + day +'/Sample_L' +  \
                    treatment + strain + rep
                if variant_type == 'INS':
                    drop_dups_on = ['type', 'seq_id', 'position']
                    IN_merged = IN_merged.drop_duplicates(subset=drop_dups_on, keep="first")
                if not os.path.exists(out_path):
                    os.makedirs(out_path)
                #print IN_variants_SNPs
                IN_merged.to_csv(out_path + '/' + variant_type +'.txt', sep = '\t', index = False)


def merge_variant_annotated(day, strain, variant_type):
    treatments = ['0', '1', '2']
    reps = ['1', '2', '3', '4', '5']
    count = 0
    for treatment in treatments:
        for rep in reps:
            path = mydir + 'breseq_output_gbk_essentials_split_clean/' + day \
                + '/Sample_L' + treatment + strain + rep + '/' + variant_type +'.txt'
            if os.path.exists(path) == True:
                IN = pd.read_csv(path, sep = '\t', header = 'infer')
                renamed = []
                for x in to_rename:
                    x_renamed = x + '_L' + treatment +  strain + rep + '_' + day
                    IN = IN.rename(columns = {x : x_renamed})
                    renamed.append(x_renamed)

                if count == 0:
                    merged = IN
                    frequency = 'frequency_L' + treatment +  strain + rep + '_' + day
                    merged_freq = merged[frequency]
                    merged.drop(labels=[frequency], axis=1,inplace = True)
                    merged.insert(len(merged.columns)-1, frequency, merged_freq)
                else:
                    merged_keep = renamed + merged_on

                    merged = pd.merge(merged, IN[merged_keep], \
                            how='outer', on = merged_on)
                    #print merged
                count += 1
        test = merged.columns.tolist()
        for i, column in enumerate(merged_on):
            test.remove(column)
            test.insert(i, column)
        merged = merged.reindex_axis(test, axis=1)
        OUT_path = mydir + 'breseq_output_gbk_essentials_split_clean_merged/' \
                + day + '/' + strain
        if not os.path.exists(OUT_path):
            os.makedirs(OUT_path)
        OUT_name = OUT_path + '/Strain_' + strain + '_' + day + '_' + variant_type + '.txt'
        #if 'seq_id' in merged.columns:
        #    merged = merged.drop_duplicates(subset = ['position', 'seq_id'])
        #else:
        #    merged = merged.drop_duplicates(subset = ['position'])
        merged.to_csv(OUT_name, sep = '\t', index = False)


def get_unique_mutations(day, strain, variant_type):
    path = mydir + 'breseq_output_gbk_essentials_split_clean_merged/' + day +'/' \
        + strain + '/Strain_' + strain + '_' + day + '_' + variant_type + '.txt'
    IN = pd.read_csv(path, sep = '\t', header = 'infer')
    if not os.path.exists(mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique/' \
            + day + '/'):
        os.makedirs(mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique/' \
            + day + '/')
    sample_freqs = [x for x in IN.columns if 'frequency_' in x]
    IN_freqs = IN[sample_freqs]
    samples = IN_freqs.shape[1]
    NoDups = IN_freqs[IN_freqs.apply(lambda x:  x.isnull().sum() == samples - 1, 1)]
    NoDups_index = NoDups.index.values
    IN_unique = IN.ix[NoDups_index]
    OUTname = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique/' \
        + day + '/Strain_' + strain + '_' + day + '_' + variant_type + '.txt'
    IN_unique.to_csv(OUTname, sep = '\t', index = False)


def get_sample_by_gene_matrix(strain):
    table_out_strain = mydir + 'gene_by_sample/' + strain
    if not os.path.exists(table_out_strain):
        os.makedirs(table_out_strain)
    path = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique_merged/Strain_' \
        + strain + '_SNP.txt'
    IN = pd.read_csv(path, sep = '\t', header = 'infer')
    sample_freqs = [x for x in IN.columns if 'frequency_' in x and '_frequency_' not in x]
    list_slice = ['seq_id'] + ['position'] + ['gene_name'] + ['locus_tag'] + sample_freqs
    IN_slice = IN[list_slice]
    mut_dict = {}
    for index, row in IN_slice.iterrows():
        gene_name = row['locus_tag']
        if gene_name not in mut_dict:
            mut_dict[gene_name] = {}
        for sample_freq in sample_freqs:
            # get the frequency
            sample_freq_row = row[sample_freq]
            # check whether or not there's a mutation
            if np.isnan(sample_freq_row) == True:
                sample_freq_row = 0
            else:
                sample_freq_row = 1
            # check whether the sample is in the dictionary for the given gene
            # if the sample isn't in the dict, add the sample
            if sample_freq not in mut_dict[gene_name]:
                mut_dict[gene_name][sample_freq] = sample_freq_row
            # otherwise add the gene to the existing value
            else:
                mut_dict[gene_name][sample_freq] += sample_freq_row
    mut_df = pd.DataFrame.from_dict(mut_dict, orient='index')
    #'''rename the sample names?'''
    #rename_dict {}
    #for sample_freq in sample_freqs:
    #    sample_freq_new = sample_freq.split('_')
    #    rename_dict[sample_freq_new] = sample_freq_new[1]
    #mut_df = mut_df.rename(columns = {'two':'new_name'})
    mut_df = mut_df.T
    # add columns of zero for genes with no observed mutations
    gene_path = mydir + 'data/reference_assemblies_task2/reference_assemblies_task2_table/B.txt'
    IN_gene = pd.read_csv(gene_path, sep = ' ', header = 'infer')
    print(IN_gene)
    # use set theory to get genes that aren't listed in matrix
    # add those genes as empty columns

    OUTpath = mydir + 'gene_by_sample/' + strain
    if not os.path.exists(OUTpath):
        os.makedirs(OUTpath)
    OUTname1 = OUTpath + '/sample_by_gene.txt'
    mut_df.to_csv(OUTname1, sep = '\t', index = True)
    # remove singletons
    multiple_hits = mut_df.loc[:, mut_df.sum(axis=0) > 1]
    OUTname2 = OUTpath + '/sample_by_gene_multiple_hits.txt'
    multiple_hits.to_csv(OUTname2, sep = '\t', index = True)


def get_sample_by_gene_matrix_gscore(strain):
    if strain == 'S':
        gene_path = mydir + 'reference_assemblies_task2/reference_assemblies_task2_table/B.txt'
    else:
        gene_path = mydir + 'reference_assemblies_task2/reference_assemblies_task2_table/' + strain + '.txt'
    IN_gene = pd.read_csv(gene_path, sep = ' ', header = 'infer')
    gene_by_pop_path = mydir + 'gene_by_sample/' + strain + '/sample_by_gene.txt'
    IN_gbp = pd.read_csv(gene_by_pop_path, sep = '\t', header = 'infer', index_col = 0)
    # remove intergenic regions
    cols = [c for c in IN_gbp.columns if '/' not in c]
    IN_gbp = IN_gbp[cols]
    L_tot = sum(IN_gene.Size.values)
    size_dict = dict(zip(IN_gene.LocusTag, IN_gene.Size))
    rel_sizes = np.asarray([size_dict[x] for x in IN_gbp.columns]) / L_tot
    for index, row in IN_gbp.iterrows():
        E_i = sum(row) * rel_sizes
        E_i_and_G_i = zip(E_i, row)
        G_i = [2*x[1]*np.log(x[1]/x[0]) for x in E_i_and_G_i]
        G_i = [0 if math.isnan(x) else x for x in G_i]
        # make negatives zero for now. Need to worry about this later
        #num_neg = [x for x in G_i if x < 0]
        #if len(num_neg) > 0:
        #    print len(num_neg)
        G_i = [0 if x < 0 else x for x in G_i]
        #value_when_true if condition else value_when_false
        IN_gbp.loc[index,:] = G_i

    OUTname = mydir + 'gene_by_sample/' + strain + '/sample_by_gene_Gscore.txt'
    IN_gbp.to_csv(OUTname, sep = '\t', index = True)



def merge_B_S_sample_by_gene_matrix_gscore():
    B_path = mydir + 'gene_by_sample/B/sample_by_gene_Gscore.txt'
    S_path = mydir + 'gene_by_sample/S/sample_by_gene_Gscore.txt'
    B = pd.read_csv(B_path, sep = '\t', header = 'infer', index_col = 0)
    S = pd.read_csv(S_path, sep = '\t', header = 'infer', index_col = 0)
    #B_S_merged = pd.concat([B, S], axis = 1, join='outer')
    B_S_merged = B.append(S)
    B_S_merged = B_S_merged.fillna(0.0)
    #B_names = B.columns.values
    #S_names = S.columns.values
    #print len(S_names)
    #print len(B_names)
    #print set(B_names + S_names)
    OUTname = mydir + 'gene_by_sample/B_S/sample_by_gene_Gscore.txt'
    B_S_merged.to_csv(OUTname, sep = '\t', index = True)

def get_column_name(row, row_name):
    names = row.index.tolist()
    frequency_samples = [y for y in names if 'frequency_' in y]
    frequency_row = row[frequency_samples]
    sample = frequency_row[frequency_row.notnull()].index.tolist()[0]
    if row_name == 'treatment':
        return sample[-3]
    elif row_name == 'replicate':
        return sample[-1]
    else:
        return float('NaN')


def merge_unique_mutations(days, strain, variant_type):
    days_dfs = []
    count = 0
    for day in days:
        path = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique/' \
            + day + '/Strain_' + strain + '_' + day + '_' + variant_type + '.txt'
        if os.path.exists(path):
            IN = pd.read_csv(path, sep = '\t', header = 'infer')
        else:
            continue
        if count == 0:
            merged = IN
        else:
            merged = pd.merge(merged, IN, \
                    how='outer', on = merged_on)
        count += 1
    test = merged.columns.tolist()
    for i, column in enumerate(merged_on):
        test.remove(column)
        test.insert(i, column)
    merged = merged.reindex_axis(test, axis=1)
    OUT_name = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique_merged/' \
            + 'Strain_' + strain + '_'+ variant_type + '.txt'
    if 'seq_id' in merged.columns:
        merged = merged.drop_duplicates(subset = ['position', 'seq_id'])
    else:
        merged = merged.drop_duplicates(subset = ['position'])
    merged.to_csv(OUT_name, sep = '\t', index = False)



def get_split_unique(day, strain, variant_type):
    path = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique/' + day + '/Strain_' \
            + strain + '_' + day + '_' + variant_type + '.txt'
    IN = pd.read_csv(path, sep = '\t', header = 'infer')
    if not os.path.exists(mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique_split/' + day + '/'):
        os.makedirs(mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique_split/' + day + '/')
    if not os.path.exists(mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique_split/' + day + '/' + strain):
        os.makedirs(mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique_split/' + day + '/' + strain)
    samples = [x.split('_')[1] for x in IN.columns if 'frequency_' in x]
    for sample in samples:
        to_rename_sample = [x + '_' + sample + '_' + day for x in to_rename]
        to_slice = merged_on + to_rename_sample
        IN_sample = IN[to_slice]
        IN_sample = IN_sample[np.isfinite(IN_sample['bias_p_value_' + sample + '_' + day])]
        OUTname = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique_split/' + day + '/' \
                    + strain + '/' + sample + '_' + day + '_' + variant_type + '.txt'
        IN_sample.to_csv(OUTname, sep = '\t', index = False)


def run_everything(day, strain, split = False, get_variants = False, \
    merge_variants = False,  unique_mutations = False, \
    split_unique = False, variant_type = 'SNP'):
    treatments = ['0', '1', '2']
    reps = ['1', '2', '3', '4', '5']
    table_out = mydir + 'gene_by_sample'
    if not os.path.exists(table_out):
        os.makedirs(table_out)
    for treatment in treatments:
        for rep in reps:
            if split == True:
                if not os.path.exists(mydir + 'breseq_output_gbk_essentials_split/' + day + '/'):
                    os.makedirs(mydir + 'breseq_output_gbk_essentials_split/' + day + '/')
                path = mydir + 'breseq_output_gbk_essentials/' + day + '/Sample_L' \
                        + treatment + strain + rep + '-' + day[1:]
                evidence_path = path + '/evidence.gd'
                annotated_path = path + '/annotated.gd'
                if os.path.exists(evidence_path) == True:
                    cleanBreseq_annotated(annotated_path).split_annotated()
            if get_variants == True:
                get_variant_annotated(day, strain, treatments, reps, variant_type)
    if merge_variants == True:
        #print "merge_variants"
        merge_variant_annotated(day, strain, variant_type)
    if unique_mutations == True:
        #print "unique_mutations"
        get_unique_mutations(day, strain, variant_type)

get_sample_by_gene_matrix('B')
