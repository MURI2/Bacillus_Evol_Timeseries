from __future__ import division
import cleanData as cd
import pandas as pd
import os
import scipy.stats as stats
import numpy as np

mydir = os.path.expanduser("~/GitHub/Task2/PoolPopSeq/data/")

ts_tv_dict = {
    ("A", "C"):"V", ("A", "G"):"S", ("A", "T"):"V",
    ("C", "A"):"V", ("C", "G"):"V", ("C", "T"):"S",
    ("G", "A"):"S", ("G", "C"):"V", ("G", "T"):"V",
    ("T", "A"):"V", ("T", "C"):"S", ("T", "G"):"V",
    }

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


bp_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}

mut_bp_dict = {"B" : 3.35 * (10**-10), "C" : 3.35 * (10**-10), \
            "D" : 3.35 * (10**-10), "J" : 3.35 * (10**-10), "S" : 3.35 * (10**-10)}

mut_id_dict = {"B" : 1.20 * (10**-10), "C" : 3.46 * (10**-10), \
            "D" : 3.35 * (10**-10), "J" : 3.35 * (10**-10), "S" : 3.35 * (10**-10)}

mut_bias_dict = {"B" : 1.2739, "C" : 0.5011, "D" : 0.5511, "J" : 2.9555, "S" : 1.2739}



L_samples = {'C': 4042929, 'D': 3284156, 'F': 5836693, 'P':6592875, \
    'B':4215606, 'J':6082545, 'S':4215606,}

species_dict = {'B': 'Bacillus', 'C':'Caulobacter', 'D':'Deinococcus', \
    'F': 'Pedobacter', 'J': 'Janthinobacterium', 'P':'Pseudomonas', \
    'S': 'Bacillus_spoA'}



def getTransferTime(x):
    x = x.split('_')[1][1]
    if x == str(0):
        return str(1)
    elif x == str(1):
        return str(10)
    elif x == str(2):
        return str(100)

def common_entries(*dcts):
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)

def g_b_s_long(strain):
    IN_file_mut = mydir + 'gene_by_sample/' + strain +  '/sample_by_gene.txt'
    if strain == 'S':
        IN_file_genes = mydir + 'reference_assemblies_task2/reference_assemblies_task2_table/B.txt'
    else:
        IN_file_genes = mydir + 'reference_assemblies_task2/reference_assemblies_task2_table/' \
            + strain + '.txt'
    IN_mut = pd.read_csv(IN_file_mut, sep = '\t')
    # remove columns with all zeros
    #print IN_mut.loc[:, (IN_mut != 0).any(axis=0)]
    IN_genes = pd.read_csv(IN_file_genes, sep = ' ')
    IN_mut = IN_mut.T.reset_index()
    IN_mut.columns = IN_mut.iloc[0]
    IN_mut = IN_mut.reindex(IN_mut.index.drop(0))
    IN_mut.rename(columns={'Unnamed: 0':'LocusTag'}, inplace=True)
    IN_merged = pd.merge(IN_mut, IN_genes, how='inner',
        on=['LocusTag'])
    IN_merged['SizeLog10'] = np.log10(IN_merged['Size'])
    # turn the wide form data into long form
    value_vars_strain = [x for x in IN_merged.columns if 'frequency_' in x]
    #print value_vars_strain
    IN_merged_long = pd.melt(IN_merged, id_vars=['LocusTag', 'SizeLog10', \
        'GC', 'Sequence'], value_vars=value_vars_strain, \
        var_name = 'Line', value_name = 'Muts')
    IN_merged_long['TransferTime'] = IN_merged_long['Line'].apply(getTransferTime)
    IN_merged_long['Strain'] = IN_merged_long['Line'].apply(lambda x: x.split('_')[1][-2])
    IN_merged_long['Replicate'] = IN_merged_long['Line'].apply(lambda x: x.split('_')[1][-1])
    IN_merged_long['Day'] = IN_merged_long['Line'].apply(lambda x: x.split('_')[2][1:])
    # make sure mutations are read as integers
    IN_merged_long.Muts = IN_merged_long.Muts.astype(int)
    IN_merged_long.TransferTime = IN_merged_long.TransferTime.astype(int)
    #IN_merged_long = IN_merged_long.sort_values(['TransferTime', 'Replicate'], \
    #    ascending=[True, True])
    #cov_matrix = IN_merged_long[['TransferTime', 'Replicate']]
    #cov_matrix['TransferTime'] = cov_matrix['TransferTime'].apply(int)
    #cov_matrix['Replicate'] = cov_matrix['Replicate'].apply(int)
    #cov_matrix['TransferTime'] = np.log10(cov_matrix['TransferTime'])
    #cov_matrix['Replicate'] = cov_matrix['Replicate'] - 1
    #cov_matrix['TransferTime'] = cov_matrix['TransferTime'].apply(int)
    IN_merged_long.to_csv(mydir + 'gene_by_sample_long/gene_by_sample_long_' + strain + '.txt', sep = '\t')


def p_q(day, strain):
    strain_path = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique_split/D100/' \
                + strain
    if not os.path.exists(mydir + 'p_q/' + day):
        os.makedirs(mydir + 'p_q/' + day)
    if not os.path.exists(mydir + 'p_q/' + day + '/' + strain):
        os.makedirs(mydir + 'p_q/' + day + '/' + strain)
    for pop in os.listdir(strain_path):
        IN = pd.read_csv(strain_path + '/' + pop, sep = '\t')
        sample_list = [x for x in IN.columns if 'snp_type_' in x]
        not_nuc = ['.', 'N']
        nuc_list = ['A', 'C', 'G', 'T']
        # dictionaries containing the number of polymorphis at each gene for different n-fold sites
        q_1 = {}
        q_2 = {}
        q_2V = {}
        q_2S = {}
        q_3 = {}
        q_4 = {}
        p_1 = {}
        p_2 = {}
        p_2V = {}
        p_2S = {}
        p_3 = {}
        p_4 = {}
        # q = transversion
        # p = transition
        dict_list = [p_1, p_2, p_2V, p_2S, p_3, p_4, q_1, q_2, q_2V, q_2S, q_3, q_4]
        for index, row in IN.iterrows():
            # get sample designation
            sample = [x for x in sample_list if pd.isnull(row[x]) == False][0].split('_')[2]
            ref_codon = row['codon_ref_seq_' + sample]
            new_codon = row['codon_new_seq_' + sample]
            if pd.isnull(ref_codon) == False and pd.isnull(new_codon) == False and \
                (row['mutation_' + sample] not in not_nuc and \
                row['reference'] not in not_nuc ):
                if len(set(list(ref_codon)) & set(not_nuc)) != 0:
                    continue
                codon_pos = int(row['codon_position_' + sample])
                if row['gene_strand'] == '<':
                    ref = bp_dict[row['reference']]
                    mut = bp_dict[row['mutation_' + sample]]
                else:
                    ref = row['reference']
                    mut = row['mutation_' + sample]
                pos = int(row['codon_position_' + sample]) - 1
                locus_tag = row['locus_tag']
                s_v = ts_tv_dict[(mut, ref)]
                ref_aa = 0
                for dict_i in dict_list:
                    #if sample not in dict_i:
                    #    dict_i[sample]:
                    if locus_tag not in dict_i:
                        dict_i[locus_tag] = 0
                ref_aa =  row['aa_ref_seq_' + sample]
                mut_aa =  row['aa_new_seq_' + sample]

                fold_count = 0
                fold_2_S_i = 0
                fold_2_V_i = 0
                # is the site 1, 2, 2V, 2S, 3, or 4 fold redundant?
                for nuc in nuc_list:
                    ref_codon_list = list(ref_codon)
                    if ref_codon_list[pos] == nuc:
                        continue
                    ref_codon_list[pos] = nuc
                    codon_nuc_mut = "".join(ref_codon_list)
                    S_V_i = ts_tv_dict[(ref_codon[pos], codon_nuc_mut[pos])]
                    #print ref_codon[pos], codon_nuc_mut[pos], S_V
                    if codon_dict[ref_codon] == codon_dict[codon_nuc_mut]:
                        fold_count += 1
                        if S_V_i == 'S':
                            fold_2_S_i += 1
                        elif S_V_i == 'V':
                            fold_2_V_i += 1
                # check wheter mutation was transition of transversion
                mut_S_V = ts_tv_dict[(ref_codon[pos], new_codon[pos])]
                # transitions = p
                if mut_S_V == 'S':
                    if fold_count == 3:
                        p_4[locus_tag] += 1
                    elif fold_count == 2:
                        p_3[locus_tag] += 1
                    elif fold_count == 1:
                        p_2[locus_tag] += 1
                        if fold_2_S_i == 1 and fold_2_V_i == 0:
                            p_2S[locus_tag] += 1
                        elif fold_2_S_i == 0 and fold_2_V_i == 1:
                            p_2V[locus_tag] += 1
                        else:
                            print fold_S_count, fold_V_count
                    elif fold_count == 0:
                        p_1[locus_tag] += 1
                # transversions = q
                if mut_S_V == 'V':
                    if fold_count == 3:
                        q_4[locus_tag] += 1
                    elif fold_count == 2:
                        q_3[locus_tag] += 1
                    elif fold_count == 1:
                        q_2[locus_tag] += 1
                        if fold_2_S_i == 1 and fold_2_V_i == 0:
                            q_2S[locus_tag] += 1
                        elif fold_2_S_i == 0 and fold_2_V_i == 1:
                            q_2V[locus_tag] += 1
                        else:
                            print fold_S_count, fold_V_count
                    elif fold_count == 0:
                        q_1[locus_tag] += 1
        dict_merged = list(common_entries(p_1, p_2, p_2V, p_2S, p_3, p_4, \
                                            q_1, q_2, q_2V, q_2S, q_3, q_4))
        if len(dict_merged) == 0:
            continue
        OUT =  mydir + 'p_q/' + day + '/' + strain + '/' + pop
        OUT = open(OUT, 'w')
        print>> OUT, 'locus_tag', 'p_1', 'p_2', 'p_2V', 'p_2S', 'p_3', 'p_4', \
                    'q_1', 'q_2', 'q_2V', 'q_2S', 'q_3', 'q_4'
        for x in dict_merged:
            print>> OUT, x[0], x[1], x[2], x[3], x[4], x[5], \
                        x[6], x[7], x[8], x[9], x[10], x[11], x[12]
        OUT.close()

def dN_dS(day):
    strains = ['B', 'C', 'D', 'F', 'J', 'P']
    OUT =  mydir + 'pop_gen_stats/D100/dN_dS.txt'
    OUT = open(OUT, 'w')
    print>> OUT, 'strain', 'treatment', 'replicate', 'dN', 'dS', 'dN/dS'
    for strain in strains:
        print strain
        strain_L_path = mydir + 'reference_assemblies_task2_table/' + strain + '.txt'
        strain_L = pd.read_csv(strain_L_path, sep = ' ')
        strain_path = mydir + 'p_q/' + day + '/' + strain
        if not os.path.exists(mydir + 'dN_dS/' + day):
            os.makedirs(mydir + 'dN_dS/' + day)
        if not os.path.exists(mydir + 'dN_dS/' + day + '/' + strain):
            os.makedirs(mydir + 'dN_dS/' + day + '/' + strain)
        for pop in os.listdir(strain_path):
            if pop == '.DS_Store':
                continue
            IN_p_q_file = mydir + 'p_q/' + day + '/' + strain + '/' + pop
            IN_p_q = pd.read_csv(IN_p_q_file, sep = ' ')
            Q_S = (sum(IN_p_q.q_2V.values) + sum(IN_p_q.q_4.values)) / \
                    (sum(strain_L.Fold_2_V.values) + sum(strain_L.Fold_4.values))
            B_s = -(1/2) * np.log(1 - (2*Q_S))
            Q_A = (sum(IN_p_q.q_1.values) + sum(IN_p_q.q_2S.values)) / \
                    (sum(strain_L.Fold_1.values) + sum(strain_L.Fold_2_S.values))
            B_a = -(1/2) * np.log(1 - (2*Q_A))
            P_2S = sum(IN_p_q.p_2S.values) / sum(strain_L.Fold_2_S.values)
            P_4 = sum(IN_p_q.p_4.values) / sum(strain_L.Fold_4.values)

            A_2S = (-(1/2) * np.log(1 - (2*P_2S) - Q_A)) + ((1/4) * np.log(1 - (2*Q_A)))
            Q_4 = sum(IN_p_q.q_4.values) /  sum(strain_L.Fold_4.values)
            A_4 = (-(1/2) * np.log(1 - (2*P_4) - Q_4)) + ((1/4) * np.log(1 - (2*Q_4)))
            A_s = ((sum(strain_L.Fold_2_S.values) * A_2S) + (sum(strain_L.Fold_4.values) * A_4)) /\
                    (sum(strain_L.Fold_2_S.values) + sum(strain_L.Fold_4.values))

            P_2V = sum(IN_p_q.p_2V.values) / sum(strain_L.Fold_2_V.values)
            P_0 = sum(IN_p_q.p_1.values) / sum(strain_L.Fold_1.values)
            Q_0 = sum(IN_p_q.q_1.values) /  sum(strain_L.Fold_4.values)
            A_2V = (-(1/2) * np.log(1 - (2*P_2V) - Q_S)) + ((1/4) * np.log(1 - (2*Q_S)))
            A_0 = (-(1/2) * np.log(1 - (2*P_0) - Q_0)) + ((1/4) * np.log(1 - (2*Q_0)))
            A_a = ((sum(strain_L.Fold_2_V.values) * A_2V) + (sum(strain_L.Fold_1.values) * A_0)) / \
                    (sum(strain_L.Fold_2_V.values) + sum(strain_L.Fold_1.values))

            dS = A_s + B_s
            dN = A_a + B_a

            print>> OUT, species_dict[pop[2]], pop[1], pop[3], str(dN), str(dS), str(dN/dS)

    OUT.close()


def mut_bias():
    mut_strains = ['B', 'C', 'D', 'J', 'S']
    OUT =  mydir + 'pop_gen_stats/D100/mut_bias.txt'
    OUT = open(OUT, 'w')
    print>> OUT, 'strain', 'treatment', 'replicate', 'm_sample_ma'
    for strain in mut_strains:
        IN_file = mydir + 'breseq_output_gbk_essentials_split_clean_merged_unique/D100/Strain_' + strain +  '_D100_SNP.txt'
        IN = pd.read_csv(IN_file, sep = '\t')
        AT_GC = {}
        GC_AT = {}
        mut_bias_sample_dict = {}
        sample_snp_list = [x for x in IN.columns if 'snp_type_' in x]
        sample_list = [x.split('_')[2] for x in sample_snp_list]
        for sample in sample_list:
            AT_GC[sample] = 0
            GC_AT[sample] = 0
            mut_bias_sample_dict[sample] = [0, 0]
        for index, row in IN.iterrows():
            sample_row = [x for x in sample_snp_list if pd.isnull(row[x]) == False][0].split('_')[2]
            ref = row['reference']
            mut = row['mutation_' + sample_row + '_D100']
            #print sample_row, ref, mut
            if (ref == 'A' and mut == 'C') or \
                (ref == 'A' and mut == 'G') or \
                (ref == 'T' and mut == 'C') or \
                (ref == 'T' and mut == 'G'):
                AT_GC[sample_row] += 1
            elif (ref == 'C' and mut == 'A') or \
                (ref == 'C' and mut == 'T') or \
                (ref == 'G' and mut == 'A') or \
                (ref == 'G' and mut == 'T'):
                GC_AT[sample_row] += 1
            else:
                continue
        AT_GC_list = list(common_entries(GC_AT, AT_GC))
        AT_GC_dict = {}
        for x in AT_GC_list:
            if (x[1] == 0 and x[2] == 0):
                continue
            else:
                print x
                AT_GC_dict[x[0]] = round(((x[1] + 1) / (x[2] + 1)) / mut_bias_dict[strain], 3)
                #AT_GC_dict[x[0]] = round((x[1] + 1) / (x[2] + 1), 3)
        for key, value in AT_GC_dict.iteritems():
            #print key, value
            print>> OUT, species_dict[key[2]], key[1], key[3], value
    OUT.close()


def run_everything():
    strains = ['B', 'C', 'D', 'F', 'J', 'P', 'S']
    #strains = ['B']
    days = ['D100', 'D200', 'D300']
    #days = ['D100']
    #strains = ['P']
    #variant_types = ['SNP', 'INS', 'DEL']
    variant_types = ['SNP']
    for strain in strains:
    #    for variant_type in variant_types:
    #        for day in days:
    #            if strain == 'S' and day != 'D100':
    #                continue
    #            print day, strain, variant_type
    #            cd.run_everything(day, strain, split = False, get_variants = False, \
    #                merge_variants = False,  unique_mutations = False, \
    #                split_unique = False, variant_type = variant_type)
    #        #cd.merge_unique_mutations(days, strain, variant_type)
    #        #if variant_type == 'SNP':
    #        #    cd.get_sample_by_gene_matrix(strain)
    #    #g_b_s_long(strain)
        if strain != 'S':
            cd.cleanGBK(strain)
    #    cd.get_sample_by_gene_matrix_gscore(strain)

#cd.cleanGBK(strain)
#p_q('D100', strain)

#cd.get_sample_by_gene_matrix('D100', strain)
#dN_dS('D100')
mut_bias()
#run_everything()
#cd.merge_B_S_sample_by_gene_matrix_gscore()
