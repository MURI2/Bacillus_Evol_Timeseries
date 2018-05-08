from __future__ import division
import os, re
import bacillus_tools as bt
import numpy as np
import pandas as pd
from skbio.diversity import beta_diversity
import  matplotlib.pyplot as plt

class likelihood_matrix:

    def get_gene_lengths(self):
        gene_df_path = bt.get_path() + '/data/gene_table.txt'
        gene_df = pd.read_csv(gene_df_path, sep = '\t', header = 'infer', index_col = 0)
        length_df = gene_df.loc[:, 'Size']
        length_df.columns = ['Size']
        return length_df
        #return pd.Series(gene_df.Size.values, index=gene_df.index).to_dict()

    def get_likelihood_matrix(self):
        df_in = bt.get_path() + '/data/pool_pop_seq/gene_by_pop.txt'
        df = pd.read_csv(df_in, sep = '\t', header = 'infer', index_col = 0)
        genes = df.columns.tolist()
        genes_lengths = self.get_gene_lengths().loc[genes]
        genes_lengths = genes_lengths.reindex(genes)
        L_mean = np.mean(list(genes_lengths.values))
        L_i = np.asarray(list(genes_lengths.values))
        N_genes = len(genes)
        m_mean = df.sum(axis=1) / N_genes
        for index, row in df.iterrows():
            m_mean_j = m_mean[index]
            delta_j = row * np.log((row * (L_mean / L_i)) / m_mean_j)
            df.loc[index,:] = delta_j
        out_name = bt.get_path() + '/data/pool_pop_seq/gene_by_pop_delta.txt'
        df_new = df.fillna(0)
        # remove colums with all zeros
        df_new.loc[:, (df_new != 0).any(axis=0)]
        # replace negative values with zero
        df_new[df_new < 0] = 0
        df_new.to_csv(out_name, sep = '\t', index = True)


def pcoa():
    print('what')
    # can't use floats, just write your own functions
    data = [[0.5, 64, 14, 0, 0, 3, 1],
         [0, 3, 35, 42, 0, 12, 1],
         [0, 5, 5, 0, 40, 40, 0],
         [44, 35, 9, 0, 1, 0, 0],
         [0, 2, 8, 0, 35, 45, 1],
         [0, 0, 25, 35, 0, 19, 0]]
    ids = list('ABCDEF')
    print(beta_diversity("braycurtis", data, ids))
    #pw_distances(data, ids, "braycurtis")


class mut_bias:

    def get_mut_bias(self):
        out_df = open(bt.get_path() + '/data/mut_bias.txt', 'w')
        out_df.write('\t'.join(['Sample', 'Strain', 'Treatment', 'Replicate', 'Time' ,'m_sample_ma']) + '\n')
        AT_GC = {}
        GC_AT = {}
        to_exclude = bt.mutations_to_exclude()
        gene_pop_matrix = {}
        directory = os.fsencode(bt.get_path() + '/data/pool_pop_seq/rebreseq_annotated')
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            if filename.endswith('-100.gd'):
                in_df = open(os.path.join(str(directory, 'utf-8'), filename), 'r')
                pop = filename.split('.')[0]
                if pop not in AT_GC:
                    AT_GC[pop] = 0
                if pop not in GC_AT:
                    GC_AT[pop] = 0
                to_keep = []
                for line in in_df:
                    line_split = line.strip().split()
                    if line_split[0] == 'SNP':
                        to_keep.append(int(line_split[2]))
                    # all RA occur after SNPs
                    if (line_split[0] == 'RA') and (int(line_split[1]) in to_keep):
                        ref = line_split[6]
                        mut = line_split[7]
                        if (ref == 'A' and mut == 'C') or \
                            (ref == 'A' and mut == 'G') or \
                            (ref == 'T' and mut == 'C') or \
                            (ref == 'T' and mut == 'G'):
                            AT_GC[pop] += 1
                        elif (ref == 'C' and mut == 'A') or \
                            (ref == 'C' and mut == 'T') or \
                            (ref == 'G' and mut == 'A') or \
                            (ref == 'G' and mut == 'T'):
                            GC_AT[pop] += 1
                        else:
                            continue
        AT_GC_list = list(bt.common_entries(GC_AT, AT_GC))
        AT_GC_dict = {}
        for x in AT_GC_list:
            if (x[1] == 0 and x[2] == 0):
                continue
            else:
                AT_GC_dict[x[0]] = round(((x[1] + 1) / (x[2] + 1)) / bt.get_bacillus_mut_bias(), 3)
        for key, value in AT_GC_dict.items():
            print(key, value)
            key_split = re.split(r'[-_]+', key)
            out_df.write('\t'.join([key, key_split[1][2], key_split[1][1], key_split[1][3], key_split[2], str(value)]) + '\n')
        out_df.close()

    def plot_mut_bias(self):
        df = pd.read_csv(bt.get_path() + '/data/mut_bias.txt', sep = '\t', header = 'infer', index_col = 0)
        strains = ['B', 'S']
        B_1 = df[(df.Strain == 'B') & (df.Treatment == 0)]['m_sample_ma'].tolist()
        B_10 = df[(df.Strain == 'B') & (df.Treatment == 1)]['m_sample_ma'].tolist()
        B_100 = df[(df.Strain == 'B') & (df.Treatment == 2)]['m_sample_ma'].tolist()
        S_1 = df[(df.Strain == 'S') & (df.Treatment == 0)]['m_sample_ma'].tolist()
        S_10 = df[(df.Strain == 'S') & (df.Treatment == 1)]['m_sample_ma'].tolist()
        S_100 = df[(df.Strain == 'S') & (df.Treatment == 2)]['m_sample_ma'].tolist()
        fig, ax = plt.subplots()
        ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
        ax.plot([0]* len(B_1), B_1, marker='o', linestyle='', \
            ms=14, color = bt.get_colors()['0'], alpha = 0.9)
        ax.plot([0.5]* len(S_1), S_1, marker='o', linestyle='', \
            ms=14, color = bt.get_colors()['0'], alpha = 0.9, markeredgewidth=2, mfc='none')
        ax.plot([1.5]* len(B_10), B_10, marker='o', linestyle='', \
            ms=14, color = bt.get_colors()['1'], alpha = 0.9)
        ax.plot([2]* len(S_10), S_10, marker='o', linestyle='', \
            ms=14, color = bt.get_colors()['1'], alpha = 0.9, markeredgewidth=2, mfc='none')
        ax.plot([3]* len(B_100), B_100, marker='o', linestyle='', \
            ms=14, color = bt.get_colors()['2'], alpha = 0.9)
        ax.plot([3.5]* len(S_100), S_100, marker='o', linestyle='', \
            ms=14, color = bt.get_colors()['2'], alpha = 0.9, markeredgewidth=2, mfc='none')
        #fig.canvas.draw()
        ax.text(0, 6.8, r'$m=\frac{\mathrm{G + C} \rightarrow \mathrm{A + T} }{\mathrm{A + T} \rightarrow \mathrm{G + C} }$', fontsize=18)
        labels = [item.get_text() for item in ax.get_xticklabels()]
        labels[1] = '        1-Day'
        labels[4] = '       10-Day'
        labels[7] = '      100-Day'
        plt.axhline(y=1, color='grey', linestyle='--', lw = 4)
        plt.ylim([0,8])
        ax.set_xticklabels(labels, fontsize = 18)
        ax.set_ylabel(r'$m_{pop} / m_{Ancestor}$', fontsize = 20)
        #plt.ticklabel_format(style='sci', axis='y')
        #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
        fig.savefig(bt.get_path() + '/figs/mut.png', bbox_inches='tight',  dpi = 600)
        plt.close()




mut_bias().plot_mut_bias()



#def get_hellinger():
#pca()
#likelihood_matrix().get_likelihood_matrix()
#pcoa()
